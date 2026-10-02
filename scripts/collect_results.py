#!/usr/bin/env python3
"""
PHASE 7a — walk every vcfeval output directory and emit results/results.tsv.

Joins four sources into one row per (condition x aligner x caller x variant_type):

  results/vcfeval/<tag>__<aligner>__<caller>__<set>__<type>/summary.txt
                                       -> TP, FP, FN, precision, recall, f1
  results/align_metrics.tsv            -> mapping_rate, mean_mapq, mean_depth,
                                          align_seconds, peak_rss_mb
  results/placement_accuracy.tsv       -> placement_accuracy
  logs/call_timing.tsv                 -> call_seconds

The condition axes (genome, coverage, read_length, qs_shift, seed) are parsed out
of the tag, which is why the tag encodes them: <genome>_cov30_len150_err0_seed1.

Usage: python3 scripts/collect_results.py [--set raw|filt] [--out results/results.tsv]
"""

import argparse
import csv
import gzip
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# vcfeval summary.txt columns, verified against real output:
#   1 Threshold  2 True-pos-baseline  3 True-pos-call  4 False-pos
#   5 False-neg  6 Precision          7 Sensitivity    8 F-measure
#
# TP-baseline (matched TRUTH records) is used as "TP", because precision/recall
# here are defined against the truth set. TP-call can differ when several call
# records match one truth record.
SUMMARY_COLS = 8

TAG_RE = re.compile(
    r"^(?P<genome>[^_]+)_cov(?P<coverage>[^_]+)_len(?P<read_length>[^_]+)"
    r"_err(?P<qs_shift>[^_]+)_seed(?P<seed>[^_]+)$")

# Directory name: <tag>__<aligner>__<caller>__<set>__<type>
DIR_RE = re.compile(r"^(?P<tag>.+?)__(?P<aligner>[^_]+)__(?P<caller>[^_]+)"
                    r"__(?P<set>raw|filt)__(?P<vtype>all|snps|indels)$")

TYPE_LABEL = {"snps": "SNV", "indels": "INDEL", "all": "ALL"}

# RTG's per-type ROC columns (verified against real output):
#   1 score  2 true_positives_baseline  3 false_positives
#   4 true_positives_call  5 false_negatives  6 precision
#   7 sensitivity  8 f_measure
# The LAST row is the lowest score threshold, i.e. the whole call set.
ROC_FILES = {"SNV": "snp_roc.tsv.gz", "INDEL": "non_snp_roc.tsv.gz"}


def parse_roc(path):
    """Per-type TP/FP/FN from a single full-callset vcfeval run.

    PREFERRED over pre-splitting the VCF by type. vcfeval matches variants
    haplotype-aware: it reconstructs local haplotypes and compares sequence. If
    the indels are stripped out before scoring, a SNV adjacent to an indel can no
    longer be reconciled, and vcfeval charges it as BOTH a false negative and a
    false positive. Measured on E. coli bwa/freebayes indels, pre-splitting gave
    TP=983 FP=9 FN=17 where the single run gave TP=994 FP=0 FN=6 — and 70.6% of
    the spurious FNs had another truth variant within 50 bp. See NOTES 7.6.
    """
    if not os.path.exists(path):
        return None
    last = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 8:
                last = f
    if not last:
        return None
    try:
        return f1_from_counts(dict(TP=int(float(last[1])), FP=int(float(last[2])),
                    FN=int(float(last[4])), precision=float(last[5]),
                    recall=float(last[6]), f1=float(last[7])))
    except (ValueError, IndexError):
        return None


def f1_from_counts(d):
    """F1 = 2TP / (2TP + FP + FN).

    vcfeval computes F1 from precision and recall, so when a pipeline makes NO
    calls of a type, precision is 0/0 and vcfeval reports F1 as NaN. The count
    form is the standard definition, is identical to the harmonic mean whenever
    that is defined, and correctly gives 0 when TP = 0. Precision itself is left
    as NaN in that case — it genuinely is undefined.
    """
    den = 2 * d["TP"] + d["FP"] + d["FN"]
    if d["f1"] != d["f1"] and den > 0:          # NaN check
        d["f1"] = 2 * d["TP"] / den
    return d


def parse_summary(path):
    """Return the unthresholded ('None') row — performance of the whole callset."""
    if not os.path.exists(path):
        return None
    best = None
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < SUMMARY_COLS:
                continue
            if f[0] in ("Threshold",) or f[0].startswith("---"):
                continue
            # The 'None' row is the full call set with no score threshold.
            if f[0] == "None":
                best = f
                break
            best = f          # fall back to last row if 'None' absent
    if not best:
        return None
    try:
        return f1_from_counts(dict(TP=int(float(best[1])), FP=int(float(best[3])),
                    FN=int(float(best[4])), precision=float(best[5]),
                    recall=float(best[6]), f1=float(best[7])))
    except (ValueError, IndexError):
        return None


def load_tsv(path, keyfn):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[keyfn(row)] = row
            except KeyError:
                continue
    return out


def load_call_timing(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                out[(f[0], f[1], f[2])] = f[3]
    return out


def num(row, key, default=""):
    if not row:
        return default
    v = row.get(key, default)
    return v if v not in (None, "") else default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", dest="callset", default="raw",
                    choices=["raw", "filt", "both"])
    ap.add_argument("--out", default=os.path.join(REPO, "results", "results.tsv"))
    ap.add_argument("--include-all-type", action="store_true",
                    help="also emit the combined (SNV+INDEL) rows")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any metric column is empty. The workflow "
                         "uses this so a missing join fails loudly instead of "
                         "producing rows with blank feature columns.")
    ap.add_argument("--tags", nargs="*", default=None,
                    help="only collect these run tags. The workflow passes exactly the "
                         "tags it built, so stale vcfeval directories from a different "
                         "run mode are never mixed into the table.")
    args = ap.parse_args()
    want = set(args.tags) if args.tags else None

    ve = os.path.join(REPO, "results", "vcfeval")
    if not os.path.isdir(ve):
        print(f"no vcfeval output at {ve}", file=sys.stderr)
        return 1

    align = load_tsv(os.path.join(REPO, "results", "align_metrics.tsv"),
                     lambda r: (r["tag"], r["aligner"]))
    place = load_tsv(os.path.join(REPO, "results", "placement_accuracy.tsv"),
                     lambda r: (r["tag"], r["aligner"]))
    reads = load_tsv(os.path.join(REPO, "results", "read_metrics.tsv"),
                     lambda r: r["tag"])

    # Runtime: Snakemake benchmark aggregate if present, else the legacy
    # /usr/bin/time logs written by the Phase 5/6 shell scripts.
    rt_path = os.path.join(REPO, "results", "runtime.tsv")
    align_rt, call_rt = {}, {}
    if os.path.exists(rt_path):
        with open(rt_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["stage"] == "align":
                    align_rt[(r["tag"], r["aligner"])] = r
                else:
                    call_rt[(r["tag"], r["aligner"], r["caller"])] = r
    else:
        for (tg, al, ca), sec in load_call_timing(
                os.path.join(REPO, "logs", "call_timing.tsv")).items():
            call_rt[(tg, al, ca)] = {"seconds": sec, "max_rss_mb": ""}

    cols = ["genome", "coverage", "read_length", "qs_shift", "seed",
            "aligner", "caller", "variant_type", "callset", "scoring_method",
            "TP", "FP", "FN", "precision", "recall", "f1",
            "mapping_rate", "mean_mapq", "mean_depth", "placement_accuracy",
            "align_seconds", "call_seconds", "peak_rss_mb",
            # additions beyond the brief's schema — features for the model
            "call_peak_rss_mb", "actual_coverage", "mean_q", "mean_p",
            "timing_exclusive"]

    rows = []
    for name in sorted(os.listdir(ve)):
        d = os.path.join(ve, name)
        if not os.path.isdir(d):
            continue
        m = DIR_RE.match(name)
        if not m:
            continue

        if args.callset != "both" and m["set"] != args.callset:
            continue
        if want is not None and m["tag"] not in want:
            continue

        entries = []
        if m["vtype"] == "all":
            # Primary numbers: RTG's own per-type breakdown of the single run.
            for label, roc in ROC_FILES.items():
                rr = parse_roc(os.path.join(d, roc))
                if rr:
                    entries.append((label, "single_run", rr))
            if args.include_all_type:
                sa = parse_summary(os.path.join(d, "summary.txt"))
                if sa:
                    entries.append(("ALL", "single_run", sa))
        else:
            sp = parse_summary(os.path.join(d, "summary.txt"))
            if sp:
                entries.append((TYPE_LABEL[m["vtype"]], "pre_split", sp))
        if not entries:
            print(f"WARNING: no usable summary in {name}", file=sys.stderr)
            continue

        tm = TAG_RE.match(m["tag"])
        cond = tm.groupdict() if tm else dict(
            genome=m["tag"].split("_")[0], coverage="", read_length="",
            qs_shift="", seed="")

        a = align.get((m["tag"], m["aligner"]))
        p = place.get((m["tag"], m["aligner"]))
        q = reads.get(m["tag"])
        art = align_rt.get((m["tag"], m["aligner"]))
        crt = call_rt.get((m["tag"], m["aligner"], m["caller"]))

        for vlabel, method, s in entries:
          rows.append({
            "genome": cond["genome"], "coverage": cond["coverage"],
            "read_length": cond["read_length"], "qs_shift": cond["qs_shift"],
            "seed": cond["seed"],
            "aligner": m["aligner"], "caller": m["caller"],
            "variant_type": vlabel, "callset": m["set"],
            "scoring_method": method,
            "TP": s["TP"], "FP": s["FP"], "FN": s["FN"],
            "precision": f"{s['precision']:.6f}",
            "recall": f"{s['recall']:.6f}", "f1": f"{s['f1']:.6f}",
            "mapping_rate": num(a, "mapping_rate"),
            "mean_mapq": num(a, "mean_mapq"),
            "mean_depth": num(a, "mean_depth"),
            "placement_accuracy": num(p, "placement_accuracy"),
            "align_seconds": num(art, "seconds") or num(a, "align_seconds"),
            "call_seconds": num(crt, "seconds"),
            "peak_rss_mb": num(art, "max_rss_mb") or num(a, "peak_rss_mb"),
            "call_peak_rss_mb": num(crt, "max_rss_mb"),
            "actual_coverage": num(q, "actual_coverage"),
            "mean_q": num(q, "mean_q"),
            "mean_p": num(q, "mean_p"),
            # yes = timed with the machine to itself (clean); no = timed under
            # contention, recorded for completeness but excluded from runtime
            # analysis. Both align and call must be clean for "yes".
            "timing_exclusive": ("yes" if num(art, "exclusive") == "yes"
                                 and num(crt, "exclusive") == "yes" else "no"),
        })

    order = {"SNV": 0, "INDEL": 1, "ALL": 2}
    rows.sort(key=lambda r: (r["genome"], r["scoring_method"], r["aligner"],
                             r["caller"], order.get(r["variant_type"], 9)))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.out}")

    feature_cols = ["mapping_rate", "mean_mapq", "mean_depth", "placement_accuracy",
                    "align_seconds", "call_seconds", "peak_rss_mb",
                    "actual_coverage", "mean_q", "mean_p"]
    # "nan"/"NA" count as missing too: the first pass of the timing refactor wrote
    # "nan" for peak RSS and an empty-string-only check let it through silently.
    def missing(v):
        return str(v).strip().lower() in ("", "nan", "na", "none")
    blank = {c: sum(1 for r in rows if missing(r[c])) for c in feature_cols}
    blank = {c: n for c, n in blank.items() if n}
    if blank:
        msg = ", ".join(f"{c}={n}" for c, n in blank.items())
        print(f"{'ERROR' if args.strict else 'WARNING'}: empty metric cells: {msg}",
              file=sys.stderr)
        if args.strict:
            return 2

    # 3x3 F1 matrices, printed per genome and variant type.
    genomes = sorted({r["genome"] for r in rows})
    types = [t for t in ("SNV", "INDEL", "ALL") if any(r["variant_type"] == t for r in rows)]
    aligners = ["bwa", "bowtie2", "minimap2"]
    callers = ["gatk", "freebayes", "bcftools"]

    for method in ("single_run", "pre_split"):
      if not any(r["scoring_method"] == method for r in rows):
          continue
      tail = "PRIMARY" if method == "single_run" else "comparison only (see NOTES 7.6)"
      print(f"\n########## scoring method: {method} — {tail} ##########")
      for g in genomes:
        for t in types:
            hits_any = [r for r in rows if r["genome"] == g
                        and r["variant_type"] == t and r["scoring_method"] == method]
            if not hits_any:
                continue
            print(f"\n=== {g} — {t} — F1 (aligner x caller) ===")
            print(f"{'':<10}" + "".join(f"{c:>12}" for c in callers))
            for a in aligners:
                cells = []
                for c in callers:
                    hit = [r for r in hits_any if r["aligner"] == a and r["caller"] == c]
                    cells.append(f"{float(hit[0]['f1']):.4f}" if hit else "--")
                print(f"{a:<10}" + "".join(f"{v:>12}" for v in cells))
    return 0


if __name__ == "__main__":
    sys.exit(main())
