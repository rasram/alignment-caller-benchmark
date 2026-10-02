#!/usr/bin/env python3
"""
Per-condition read metrics — ONE output row per run tag (Snakemake rule `read_metrics`).

Reports, for the paired FASTQs of one simulated run:
    actual_coverage  sequenced bases / MUTATED genome length (the template the
                     reads were drawn from — see NOTES 3.3)
    mean_q           arithmetic mean Phred score, both mates pooled
    mean_p           mean per-base error PROBABILITY, both mates pooled
    q_effective      -10*log10(mean_p)
    frac_q30         fraction of bases >= Q30

`mean_p` — not ART's arbitrary `-qs` knob — is the physically meaningful
error-rate feature for the model. Averaging Phred scores directly understates
the error rate because Q is logarithmic (NOTES 4.4).

The histogram/statistics code is imported from extract_mean_q.py rather than
re-implemented, so the two can never disagree.

Usage: read_metrics.py --tag T --r1 R1.fq --r2 R2.fq --genome-fa MUT.fa --out OUT.tsv
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_mean_q import scan_fastq, stats_from_hist  # noqa: E402

TAG_RE = re.compile(r"^(?P<genome>[^_]+)_cov(?P<coverage>[^_]+)_len(?P<read_length>[^_]+)"
                    r"_err(?P<qs_shift>[^_]+)_seed(?P<seed>[^_]+)$")

COLS = ["tag", "genome", "coverage", "read_length", "qs_shift", "seed",
        "reads", "total_bases", "mutated_genome_len", "actual_coverage",
        "mean_q", "mean_p", "q_effective", "frac_q30", "r1_mean_q", "r2_mean_q"]


def genome_len(path):
    n = 0
    with open(path) as fh:            # awk-equivalent; never grep (NOTES 3.8)
        for line in fh:
            if not line.startswith(">"):
                n += len(line.rstrip("\n"))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--genome-fa", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    m = TAG_RE.match(a.tag)
    if not m:
        sys.exit(f"tag does not encode the condition: {a.tag}")

    h1, _, c1 = scan_fastq(a.r1)
    h2, _, c2 = scan_fastq(a.r2)
    s1, s2, both = stats_from_hist(h1), stats_from_hist(h2), stats_from_hist(h1 + h2)
    if not (s1 and s2 and both):
        sys.exit(f"no quality data in {a.r1} / {a.r2}")

    glen = genome_len(a.genome_fa)
    reads = int(c1.max()) + int(c2.max())
    row = dict(m.groupdict(), tag=a.tag, reads=reads, total_bases=both["total_bases"],
               mutated_genome_len=glen,
               actual_coverage=f"{both['total_bases'] / glen:.4f}",
               mean_q=f"{both['mean_q']:.4f}", mean_p=f"{both['mean_p']:.6e}",
               q_effective=f"{both['q_effective']:.4f}", frac_q30=f"{both['frac_q30']:.6f}",
               r1_mean_q=f"{s1['mean_q']:.4f}", r2_mean_q=f"{s2['mean_q']:.4f}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.write("\t".join(str(row[c]) for c in COLS) + "\n")


if __name__ == "__main__":
    main()
