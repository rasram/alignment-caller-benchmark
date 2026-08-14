#!/usr/bin/env python3
"""
PHASE 4 — extract per-dataset base-quality statistics into a table.

WHY THIS EXISTS
---------------
The read-simulation sweep varies ART's `-qs` flag over 0, -2, -5, -10. That number
is a *quality-score shift*: an arbitrary knob, not a physical quantity. A model
fitted on "qs_shift = -5" learns nothing transferable, because -5 means nothing
outside ART.

This script converts that knob into the physical quantity it actually controls:
the **sequencing error rate**. That IS transferable — it can be measured on any
real dataset and compared against these results.

THE SUBTLETY THAT MATTERS
-------------------------
A Phred score Q encodes an error probability:  P(error) = 10 ** (-Q / 10).
The relationship is logarithmic, so **the mean of Q is not the Q of the mean error
rate**. Averaging Phred scores directly overstates data quality, because a few
terrible bases (Q2, P=0.63) are hidden by many good ones (Q40, P=0.0001).

Worked example — two bases, Q40 and Q10:
    mean Q            = (40 + 10) / 2 = 25.0        -> implies P = 0.0032
    mean P            = (0.0001 + 0.1) / 2 = 0.05005
    Q of the mean P   = -10 * log10(0.05005) = 13.0

Those differ by 12 Phred points, i.e. a ~16x difference in implied error rate.

So this script reports BOTH:
  * mean_q      arithmetic mean of Phred scores  (what FastQC shows; comparable
                to what people quote in papers)
  * mean_p      mean per-base error probability  (the physically correct average)
  * q_effective -10*log10(mean_p)                (mean_p expressed back as Phred)

`mean_p` is the feature the Phase 9/10 model should use.

Usage:
    python3 scripts/extract_mean_q.py work/*.fq
    python3 scripts/extract_mean_q.py            # defaults to baseline FASTQs
Writes: results/qc/mean_q.tsv
"""

import glob
import math
import os
import sys

import numpy as np

PHRED_OFFSET = 33  # Illumina 1.8+ / Sanger encoding
CHUNK_READS = 200_000


def scan_fastq(path):
    """Return (qual_histogram, per_position_sum, per_position_count).

    Reads only the quality lines (every 4th). Accumulates a histogram over all
    Phred values and, separately, position-wise sums so the 3' decline can be
    shown. Chunked so a multi-GB FASTQ never has to fit in memory.
    """
    hist = np.zeros(256, dtype=np.int64)
    pos_sum = np.zeros(0, dtype=np.int64)
    pos_cnt = np.zeros(0, dtype=np.int64)

    with open(path, "rb") as fh:
        buf = []
        for i, line in enumerate(fh):
            if i % 4 != 3:            # quality line is the 4th of each record
                continue
            buf.append(line.rstrip(b"\n"))
            if len(buf) >= CHUNK_READS:
                hist, pos_sum, pos_cnt = _absorb(buf, hist, pos_sum, pos_cnt)
                buf = []
        if buf:
            hist, pos_sum, pos_cnt = _absorb(buf, hist, pos_sum, pos_cnt)

    return hist, pos_sum, pos_cnt


def _absorb(buf, hist, pos_sum, pos_cnt):
    lengths = {len(b) for b in buf}
    maxlen = max(lengths)

    if len(pos_sum) < maxlen:
        pos_sum = np.pad(pos_sum, (0, maxlen - len(pos_sum)))
        pos_cnt = np.pad(pos_cnt, (0, maxlen - len(pos_cnt)))

    if len(lengths) == 1:
        # Uniform read length (the normal case): one clean 2-D view.
        arr = np.frombuffer(b"".join(buf), dtype=np.uint8).reshape(len(buf), maxlen)
        hist += np.bincount(arr.ravel(), minlength=256)
        pos_sum[:maxlen] += arr.sum(axis=0, dtype=np.int64)
        pos_cnt[:maxlen] += len(buf)
    else:
        # Ragged (e.g. trimmed data) — fall back to per-read accumulation.
        for b in buf:
            a = np.frombuffer(b, dtype=np.uint8)
            hist += np.bincount(a, minlength=256)
            pos_sum[: len(a)] += a
            pos_cnt[: len(a)] += 1

    return hist, pos_sum, pos_cnt


def stats_from_hist(hist):
    """Collapse the character histogram into Q and error-rate statistics."""
    codes = np.arange(256)
    q = codes - PHRED_OFFSET
    keep = (q >= 0) & (hist > 0)
    qv = q[keep].astype(np.float64)
    n = hist[keep].astype(np.float64)
    total = n.sum()
    if total == 0:
        return None

    mean_q = float((qv * n).sum() / total)

    p = 10.0 ** (-qv / 10.0)                     # per-base error probability
    mean_p = float((p * n).sum() / total)
    q_eff = -10.0 * math.log10(mean_p) if mean_p > 0 else float("inf")

    # Q30 rate: the standard industry summary ("% bases >= Q30").
    q30 = float(n[qv >= 30].sum() / total)
    return dict(total_bases=int(total), mean_q=mean_q, mean_p=mean_p,
                q_effective=q_eff, frac_q30=q30,
                min_q=int(qv.min()), max_q=int(qv.max()))


def main(argv):
    files = argv[1:]
    if not files:
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        files = sorted(glob.glob(os.path.join(repo, "work", "*_seed*_[12].fq")))
    if not files:
        print("no FASTQ files given or found", file=sys.stderr)
        return 1

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outdir = os.path.join(repo, "results", "qc")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "mean_q.tsv")

    rows = []
    for path in files:
        base = os.path.basename(path)
        stem = base[:-3] if base.endswith(".fq") else base
        # tag looks like: <genome>_cov30_len150_err0_seed1_1
        parts = stem.split("_")
        genome = parts[0] if parts else "?"

        def field(prefix, cast=str, default=""):
            for p in parts:
                if p.startswith(prefix):
                    try:
                        return cast(p[len(prefix):])
                    except ValueError:
                        return default
            return default

        hist, pos_sum, pos_cnt = scan_fastq(path)
        s = stats_from_hist(hist)
        if s is None:
            print(f"WARNING: no quality data in {path}", file=sys.stderr)
            continue

        with np.errstate(invalid="ignore", divide="ignore"):
            pos_mean = np.where(pos_cnt > 0, pos_sum / np.maximum(pos_cnt, 1), np.nan) - PHRED_OFFSET

        first10 = float(np.nanmean(pos_mean[:10]))
        last10 = float(np.nanmean(pos_mean[-10:]))

        # The "3' decline" must be measured from the PEAK, not from cycle 1.
        # Real Illumina reads (and ART's empirical profile) start LOW at cycle 1,
        # rise over the first ~20 cycles as the cluster signal stabilises, then
        # decay. Comparing first-10 against last-10 straddles that rise and
        # cancels most of the decline out, making good data look flat.
        peak_cycle = int(np.nanargmax(pos_mean)) + 1
        peak_q = float(np.nanmax(pos_mean))
        decline_from_peak = peak_q - last10

        rows.append(dict(
            dataset=stem, genome=genome,
            coverage=field("cov"), read_length=field("len"),
            qs_shift=field("err"), seed=field("seed"),
            mate=parts[-1],
            reads=int(pos_cnt.max()) if len(pos_cnt) else 0,
            total_bases=s["total_bases"],
            mean_q=s["mean_q"], mean_p=s["mean_p"],
            q_effective=s["q_effective"], frac_q30=s["frac_q30"],
            min_q=s["min_q"], max_q=s["max_q"],
            cycle1_q=float(pos_mean[0]),
            peak_cycle=peak_cycle, peak_q=peak_q,
            mean_q_first10=first10, mean_q_last10=last10,
            decline_from_peak=decline_from_peak,
        ))
        rows[-1]["_pos_mean"] = pos_mean

    cols = ["dataset", "genome", "coverage", "read_length", "qs_shift", "seed",
            "mate", "reads", "total_bases", "mean_q", "mean_p", "q_effective",
            "frac_q30", "min_q", "max_q", "cycle1_q", "peak_cycle", "peak_q",
            "mean_q_first10", "mean_q_last10", "decline_from_peak"]

    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    # Console summary
    print(f"{'dataset':<34} {'meanQ':>7} {'meanP':>10} {'Qeff':>6} {'%Q30':>6} "
          f"{'cyc1':>6} {'peak':>13} {'last10':>7} {'drop':>6}")
    for r in rows:
        print(f"{r['dataset']:<34} {r['mean_q']:>7.2f} {r['mean_p']:>10.3e} "
              f"{r['q_effective']:>6.2f} {100*r['frac_q30']:>5.1f}% "
              f"{r['cycle1_q']:>6.1f} "
              f"{r['peak_q']:>6.1f}@c{r['peak_cycle']:<5d} "
              f"{r['mean_q_last10']:>7.2f} {r['decline_from_peak']:>6.2f}")

    print(f"\nWrote {out}")

    # Per-position profile, printed coarsely so the 3' decline is visible in text.
    print("\nPer-position mean Q (every 10th cycle):")
    for r in rows:
        pm = r["_pos_mean"]
        idx = list(range(0, len(pm), 10))
        print(f"  {r['dataset']}")
        print("    cycle : " + " ".join(f"{i+1:>5d}" for i in idx))
        print("    meanQ : " + " ".join(f"{pm[i]:>5.1f}" for i in idx))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
