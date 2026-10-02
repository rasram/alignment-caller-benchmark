#!/usr/bin/env python3
"""
STEP 7b — WHY are the errors where they are? Three mechanism tests.

The sweep says WHICH pipelines make errors; these tests check a specific,
falsifiable explanation for each pattern instead of leaving it as a story.

  1. fp_near_indel.tsv  Bowtie2 + FreeBayes/BCFtools lose F1 as depth RISES.
     Hypothesis: Bowtie2 aligns end-to-end (no soft-clipping), so a read that
     crosses a true indel near its end is forced to absorb it as mismatches;
     with more depth, those mismatches pile up into confident false SNVs.
     Test: distance from every false-positive SNV to the nearest TRUE indel,
     against the background fraction of the genome within that window.

  2. fn_repeats.tsv     Every pipeline plateaus at ~40-50 missed SNVs, however
     deep the data. Hypothesis: those variants sit in repeats, where no read
     can be placed uniquely. Test: fraction of each missed site's reads with
     MAPQ >= 20, against the same statistic for all true variant sites.

  3. phix_errors.tsv    phiX is not perfect on every seed. Every error outside
     5x coverage is listed by position, so the cause can be read off directly.

Inputs: vcfeval outputs (results/vcfeval/), truth VCFs, alignment BAMs (work/).
Needs samtools on PATH. Outputs -> results/analysis/
Usage: python3 scripts/diagnose_errors.py
"""
import bisect
import gzip
import os
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", "analysis")
VCFEVAL = os.path.join(REPO, "results", "vcfeval")
ALIGNERS = ["bwa", "bowtie2", "minimap2"]
CALLERS = ["gatk", "freebayes", "bcftools"]
WINDOW = 150          # one read length
MAPQ_MIN = 20
SAMTOOLS = shutil.which("samtools") or sys.exit("samtools not on PATH")


def vcf_records(path):
    """(pos, ref, alt, qual) for every record of a (b)gzipped VCF."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 6)
            qual = float(f[5]) if f[5] != "." else float("nan")
            for alt in f[4].split(","):
                yield int(f[1]), f[3], alt, qual


def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1


def run_dir(tag, aligner, caller, callset="raw"):
    d = os.path.join(VCFEVAL, f"{tag}__{aligner}__{caller}__{callset}__all")
    if not os.path.exists(os.path.join(d, "summary.txt")):
        sys.exit(f"missing vcfeval output: {d} — run the sweep first")
    return d


# --- 1. false-positive SNVs vs true indels -----------------------------------
def background_fraction(indel_pos, genome_len, w):
    """Fraction of the genome within w bp of a true indel (union of windows)."""
    covered, end = 0, 0
    for p in sorted(indel_pos):
        lo, hi = max(1, p - w), min(genome_len, p + w)
        if hi > end:
            covered += hi - max(lo, end + 1) + 1
            end = hi
    return covered / genome_len


def fp_near_indel():
    truth = os.path.join(REPO, "data", "truth", "ecoli.truth.vcf.gz")
    indels = sorted(p for p, r, a, _ in vcf_records(truth) if not is_snv(r, a))
    glen = int(open(os.path.join(REPO, "data", "refs", "ecoli.fa.fai")).read().split()[1])
    bg = background_fraction(indels, glen, WINDOW)
    rows = []
    for cov in [5, 10, 20, 30, 50, 100]:
        tag = f"ecoli_cov{cov}_len150_err0_seed1"
        for a in ALIGNERS:
            for c in CALLERS:
                fps = [p for p, r, alt, _ in vcf_records(os.path.join(run_dir(tag, a, c), "fp.vcf.gz"))
                       if is_snv(r, alt)]
                dist = []
                for p in fps:
                    i = bisect.bisect_left(indels, p)
                    dist.append(min(abs(p - indels[j]) for j in (i - 1, i) if 0 <= j < len(indels)))
                near = sum(d <= WINDOW for d in dist)
                rows.append(dict(coverage=cov, aligner=a, caller=c, fp_snv=len(fps),
                                 near_indel=near,
                                 frac_near=near / len(fps) if fps else float("nan"),
                                 median_dist=statistics.median(dist) if dist else float("nan"),
                                 background_frac=bg))
    return rows


# --- 2. missed variants vs mappability ---------------------------------------
def low_mapq_fraction(bam, positions):
    """Per site: is <50% of its reads at MAPQ >= MAPQ_MIN?  (samtools depth twice)"""
    if not positions:
        return []
    bed = os.path.join(OUT, ".sites.bed")
    with open(bed, "w") as fh:
        for chrom, p in positions:
            fh.write(f"{chrom}\t{p - 1}\t{p}\n")

    def depth(extra):
        o = subprocess.run([SAMTOOLS, "depth", "-a", "-b", bed, *extra, bam],
                           check=True, capture_output=True, text=True).stdout
        return {(l.split()[0], int(l.split()[1])): int(l.split()[2]) for l in o.splitlines()}
    d_all, d_q = depth([]), depth(["-Q", str(MAPQ_MIN)])
    os.remove(bed)
    return [d_all.get(s, 0) > 0 and d_q.get(s, 0) / d_all[s] < 0.5 for s in positions]


def fn_repeats():
    truth = os.path.join(REPO, "data", "truth", "ecoli.truth.vcf.gz")
    chrom = open(os.path.join(REPO, "data", "refs", "ecoli.fa.fai")).read().split()[0]
    all_sites = sorted({(chrom, p) for p, *_ in vcf_records(truth)})
    rows = []
    for cov in [30, 100]:
        tag = f"ecoli_cov{cov}_len150_err0_seed1"
        for a in ALIGNERS:
            bam = os.path.join(REPO, "work", f"{tag}.{a}.md.bam")
            bg = low_mapq_fraction(bam, all_sites)
            for c in CALLERS:
                fn = sorted({(chrom, p) for p, *_ in
                             vcf_records(os.path.join(run_dir(tag, a, c), "fn.vcf.gz"))})
                lo = low_mapq_fraction(bam, fn)
                rows.append(dict(coverage=cov, aligner=a, caller=c, fn_sites=len(fn),
                                 fn_low_mapq=sum(lo),
                                 frac_fn_low_mapq=sum(lo) / len(lo) if lo else float("nan"),
                                 truth_sites=len(bg), frac_truth_low_mapq=sum(bg) / len(bg)))
    return rows


# --- 3. every phiX error outside 5x ------------------------------------------
def phix_errors():
    ev = defaultdict(lambda: dict(runs=0, pipelines=set(), conditions=set(), seeds=set(),
                                  quals=[]))
    n_runs = 0
    for d in sorted(os.listdir(VCFEVAL)):
        if not (d.startswith("phiX_") and d.endswith("__raw__all")) or "_cov5_" in d:
            continue
        n_runs += 1
        tag, a, c = d.split("__")[:3]
        cond, seed = tag.split("_seed")
        for kind in ("fn", "fp"):
            for p, r, alt, q in vcf_records(os.path.join(VCFEVAL, d, f"{kind}.vcf.gz")):
                e = ev[(kind.upper(), p, f"{r}>{alt}")]
                e["runs"] += 1
                e["pipelines"].add(f"{a}+{c}")
                e["conditions"].add(cond.replace("phiX_", ""))
                e["seeds"].add(int(seed))
                if kind == "fp":
                    e["quals"].append(q)
    rows = []
    for (kind, p, change), e in sorted(ev.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        callers = sorted({x.split("+")[1] for x in e["pipelines"]})
        rows.append(dict(kind=kind, pos=p, change=change, runs=e["runs"], of_runs=n_runs,
                         callers=",".join(callers), n_pipelines=len(e["pipelines"]),
                         conditions=",".join(sorted(e["conditions"])),
                         seeds=",".join(map(str, sorted(e["seeds"]))),
                         max_qual=max(e["quals"]) if e["quals"] else float("nan")))
    return rows


def write(rows, name):
    import csv
    path = os.path.join(OUT, name)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {os.path.relpath(path, REPO)} ({len(rows)} rows)")


def main():
    os.makedirs(OUT, exist_ok=True)
    write(fp_near_indel(), "fp_near_indel.tsv")
    write(fn_repeats(), "fn_repeats.tsv")
    write(phix_errors(), "phix_errors.tsv")


if __name__ == "__main__":
    main()
