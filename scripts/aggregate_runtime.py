#!/usr/bin/env python3
"""
Aggregate Snakemake `benchmark:` files into results/runtime.tsv.

Replaces the hand-rolled `/usr/bin/time -l` timing. Snakemake's benchmark
directive records wall-clock seconds, peak RSS and CPU time for each job's whole
process tree, and works identically on macOS and Linux — which /usr/bin/time does
not (BSD `-l` reports RSS in bytes, GNU `-v` in kilobytes; NOTES 5.4).

FAIRNESS (R8): every timed rule claims the whole machine via a `machine` resource
(see profiles/default/config.yaml), so no timed job ever runs concurrently with
anything else. Without that, a parallel sweep would time each tool under whatever
contention happened to be present, and the runtime comparison would be noise.

Input files (written by scripts/lib/measure.sh):
    benchmarks/measure/align/<tag>.<aligner>.tsv
    benchmarks/measure/call/<tag>.<aligner>.<caller>.tsv

Usage: aggregate_runtime.py --out results/runtime.tsv FILE [FILE ...]
"""
import argparse
import csv
import os
import sys

COLS = ["tag", "stage", "aligner", "caller", "seconds", "max_rss_mb", "cpu_seconds",
        "job_seconds", "exclusive"]


def _read_one(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty measurement file: {path}")
    return rows[-1]


def _num(r, key):
    try:
        return float(r.get(key, ""))
    except (TypeError, ValueError):
        return float("nan")


def parse(path):
    """path = benchmarks/measure/<stage>/<name>.tsv, written by scripts/lib/measure.sh.

    seconds / max_rss_mb / cpu_seconds come from /usr/bin/time around the tool
    alone (kernel-exact peak RSS). job_seconds is Snakemake's whole-job wall time
    from benchmarks/<stage>/<name>.tsv, which additionally includes activating the
    conda environment — reported so that overhead is visible, never used as the
    primary number.
    """
    stage = os.path.basename(os.path.dirname(path))          # align | call
    name = os.path.basename(path)[:-len(".tsv")]
    parts = name.split(".")
    if stage == "align" and len(parts) == 2:
        tag, aligner, caller = parts[0], parts[1], ""
    elif stage == "call" and len(parts) == 3:
        tag, aligner, caller = parts
    else:
        raise ValueError(f"unexpected measurement path: {path}")

    m = _read_one(path)
    if m.get("exit", "0") != "0":
        raise ValueError(f"measured command failed (exit {m['exit']}): {path}")

    job = float("nan")
    snk = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(path))),
                       stage, name + ".tsv")
    if os.path.exists(snk):
        job = _num(_read_one(snk), "s")
    return dict(tag=tag, stage=stage, aligner=aligner, caller=caller,
                seconds=f"{_num(m, 's'):.2f}", max_rss_mb=f"{_num(m, 'max_rss_mb'):.1f}",
                cpu_seconds=f"{_num(m, 'cpu_s'):.2f}", job_seconds=f"{job:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--timing-seeds", default="1",
                    help="comma list of seeds whose timed jobs ran machine-exclusive")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    tseeds = set(a.timing_seeds.split(","))
    rows = [parse(f) for f in a.files]
    for r in rows:
        r["exclusive"] = "yes" if r["tag"].rsplit("_seed", 1)[-1] in tseeds else "no"
    rows.sort(key=lambda r: (r["tag"], r["stage"], r["aligner"], r["caller"]))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"aggregated {len(rows)} benchmark files -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
