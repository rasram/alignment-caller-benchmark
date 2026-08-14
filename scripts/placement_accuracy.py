#!/usr/bin/env python3
"""
PHASE 5 — read placement accuracy: a pure aligner metric.

WHAT IT MEASURES
----------------
ART's `-sam` output records where every simulated read ACTUALLY came from. This
script compares that against where each aligner PUT the read, and reports the
fraction placed within a tolerance (default +/-10 bp) of the truth.

This isolates the aligner completely. Variant-calling F1 confounds two things —
how well reads were placed and how well the caller reasoned about the pileup.
Placement accuracy answers only the first, with no caller involved.

THE COORDINATE TRAP (see NOTES 3.4)
-----------------------------------
The two files are in DIFFERENT coordinate systems:

    ART truth SAM  -> MUTATED genome coordinates   (@SQ LN:5400)
    aligner BAM    -> ORIGINAL reference coords    (@SQ LN:5386)

The contig has the SAME NAME in both and the lengths differ only by the net indel
balance, so comparing them directly raises no error and no warning — it just
silently returns wrong answers. Positions drift by the cumulative indel offset,
reaching 187 bp by the end of E. coli, far beyond a 10 bp tolerance. Reads near
the end of the genome would nearly all be scored as misplaced, making every
aligner look bad in a way that gets worse along the genome.

So this script builds an explicit mutated -> reference coordinate map from the
indels simuG injected, and converts before comparing. The map is self-validated
against simuG's own `sim_start` annotations (see build_coord_map).

Usage:
    python3 scripts/placement_accuracy.py \
        --truth-sam work/<tag>_.sam \
        --indel-vcf data/truth/<gen>.refseq2simseq.INDEL.vcf \
        --bam work/<tag>.bwa.md.bam --aligner bwa \
        [--tolerance 10] [--out results/placement.tsv]
"""

import argparse
import bisect
import os
import subprocess
import sys

SAMTOOLS = os.path.join(
    os.environ.get("CONDA_BASE", os.path.expanduser("~/miniforge3")),
    "envs", "align", "bin", "samtools")

# SAM FLAG bits we need
FLAG_UNMAPPED = 0x4
FLAG_READ1 = 0x40
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def parse_indels(vcf_path):
    """Read simuG's raw INDEL VCF -> [(ref_pos, len_ref, len_alt, sim_start)].

    simuG's *raw* VCF is used rather than the normalised truth set: normalisation
    left-aligns indels, which changes their reported position by a few bases in
    homopolymers. For coordinate mapping we want the edit exactly as simuG
    applied it, and the raw file also carries simuG's own `sim_start`, which lets
    the map validate itself.
    """
    out = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            pos, ref, alt, info = int(f[1]), f[3], f[4], f[7]
            sim_start = None
            for kv in info.split(";"):
                if kv.startswith("sim_start="):
                    try:
                        sim_start = int(kv.split("=", 1)[1])
                    except ValueError:
                        pass
            out.append((pos, len(ref), len(alt), sim_start))
    out.sort(key=lambda r: r[0])
    return out


def build_coord_map(indels):
    """Build a piecewise mutated->reference coordinate map.

    Returns (starts, kinds, values) for bisect lookup, where each segment covers
    mutated positions [starts[i], starts[i+1]) and is either:
        ('shift', cum)  -> ref = mut - cum
        ('ins',  anchor)-> mutated base is INSIDE an inserted stretch and has no
                           reference equivalent; map it to the anchor ref base.

    Handling insertions explicitly matters: a plain global shift would map bases
    inside a 50 bp insertion to a position up to 50 bp away, which is larger than
    the tolerance and would show up as a fake placement error concentrated at
    insertion sites.
    """
    starts, kinds, values = [], [], []
    cum = 0
    prev_ref = 1

    for pos, lref, lalt, _sim in indels:
        delta = lalt - lref
        # Unchanged stretch: ref [prev_ref .. pos] -> mut [prev_ref+cum .. pos+cum]
        starts.append(prev_ref + cum); kinds.append("shift"); values.append(cum)
        if delta > 0:
            # Novel inserted bases occupy mut [pos+cum+1 .. pos+cum+delta]
            starts.append(pos + cum + 1); kinds.append("ins"); values.append(pos)
            prev_ref = pos + 1
        else:
            # Deleted ref bases [pos+1 .. pos+lref-1] have no mutated equivalent
            prev_ref = pos + lref
        cum += delta

    starts.append(prev_ref + cum); kinds.append("shift"); values.append(cum)
    return starts, kinds, values


def make_converter(starts, kinds, values):
    def mut_to_ref(m):
        i = bisect.bisect_right(starts, m) - 1
        if i < 0:
            return m
        if kinds[i] == "ins":
            return values[i]
        return m - values[i]
    return mut_to_ref


def validate_map(indels, mut_to_ref):
    """Self-check: simuG records each indel's own mutated coordinate (sim_start).

    Our independently-derived map must reproduce it. If this fails the map is
    wrong and every placement number would be quietly wrong too, so it is a hard
    error rather than a warning.
    """
    checked = bad = 0
    for pos, lref, lalt, sim_start in indels:
        if sim_start is None:
            continue
        checked += 1
        # The anchor base at ref `pos` sits at mutated position sim_start (for
        # deletions and 1-base-anchored insertions alike, within 1 base).
        back = mut_to_ref(sim_start)
        if abs(back - pos) > 1:
            bad += 1
    return checked, bad


def read_truth_sam(path):
    """(qname, is_read1) -> true position in MUTATED coordinates."""
    truth = {}
    opener = subprocess.Popen([SAMTOOLS, "view", path], stdout=subprocess.PIPE,
                              text=True) if path.endswith(".bam") else None
    fh = opener.stdout if opener else open(path)
    try:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 5)
            flag = int(f[1])
            if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
                continue
            truth[(f[0], bool(flag & FLAG_READ1))] = int(f[3])
    finally:
        if opener:
            fh.close(); opener.wait()
        else:
            fh.close()
    return truth


def score_bam(bam, truth, mut_to_ref, tolerance):
    """Compare aligner placements against converted truth positions."""
    n_total = n_mapped = n_correct = n_missing = 0
    mapq_sum = 0
    correct_by_mapq0 = 0
    n_mapq0 = 0
    dists = []

    p = subprocess.Popen([SAMTOOLS, "view", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
            continue
        n_total += 1
        key = (f[0], bool(flag & FLAG_READ1))
        t_mut = truth.get(key)
        if t_mut is None:
            n_missing += 1
            continue
        if flag & FLAG_UNMAPPED:
            continue                      # counted as placed incorrectly
        n_mapped += 1
        mapq = int(f[4])
        mapq_sum += mapq
        t_ref = mut_to_ref(t_mut)
        d = abs(int(f[3]) - t_ref)
        if len(dists) < 200000:
            dists.append(d)
        if d <= tolerance:
            n_correct += 1
            if mapq == 0:
                correct_by_mapq0 += 1
        if mapq == 0:
            n_mapq0 += 1
    p.wait()

    dists.sort()
    med = dists[len(dists) // 2] if dists else -1
    return dict(
        total_reads=n_total, mapped=n_mapped, missing_from_truth=n_missing,
        correct=n_correct,
        placement_accuracy=(n_correct / n_total) if n_total else 0.0,
        placement_accuracy_of_mapped=(n_correct / n_mapped) if n_mapped else 0.0,
        mapping_rate=(n_mapped / n_total) if n_total else 0.0,
        mean_mapq=(mapq_sum / n_mapped) if n_mapped else 0.0,
        mapq0_reads=n_mapq0, median_abs_offset=med,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-sam", required=True)
    ap.add_argument("--indel-vcf", required=True)
    ap.add_argument("--bam", required=True, action="append")
    ap.add_argument("--aligner", required=True, action="append")
    ap.add_argument("--tolerance", type=int, default=10)
    ap.add_argument("--tag", default="")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if len(args.bam) != len(args.aligner):
        print("--bam and --aligner must be given the same number of times",
              file=sys.stderr)
        return 2

    indels = parse_indels(args.indel_vcf)
    starts, kinds, values = build_coord_map(indels)
    mut_to_ref = make_converter(starts, kinds, values)

    checked, bad = validate_map(indels, mut_to_ref)
    if bad:
        print(f"FATAL: coordinate map disagrees with simuG sim_start for "
              f"{bad}/{checked} indels. Placement numbers would be wrong.",
              file=sys.stderr)
        return 1
    print(f"coordinate map: {len(indels)} indels, validated against simuG "
          f"sim_start on {checked} of them — all consistent")

    truth = read_truth_sam(args.truth_sam)
    print(f"truth SAM: {len(truth)} read records")

    rows = []
    for bam, aligner in zip(args.bam, args.aligner):
        r = score_bam(bam, truth, mut_to_ref, args.tolerance)
        r["aligner"] = aligner
        r["tag"] = args.tag
        r["tolerance"] = args.tolerance
        rows.append(r)

    cols = ["tag", "aligner", "tolerance", "total_reads", "mapped",
            "mapping_rate", "correct", "placement_accuracy",
            "placement_accuracy_of_mapped", "mean_mapq", "mapq0_reads",
            "median_abs_offset", "missing_from_truth"]

    print()
    print(f"{'aligner':<10} {'mapping%':>9} {'placement%':>11} "
          f"{'placed/mapped%':>15} {'meanMAPQ':>9} {'MAPQ0':>8} {'medOff':>7}")
    for r in rows:
        print(f"{r['aligner']:<10} {100*r['mapping_rate']:>8.3f}% "
              f"{100*r['placement_accuracy']:>10.3f}% "
              f"{100*r['placement_accuracy_of_mapped']:>14.3f}% "
              f"{r['mean_mapq']:>9.2f} {r['mapq0_reads']:>8d} "
              f"{r['median_abs_offset']:>7d}")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        new = not os.path.exists(args.out)
        with open(args.out, "a") as fh:
            if new:
                fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(
                    f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c])
                    for c in cols) + "\n")
        print(f"\nAppended {len(rows)} row(s) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
