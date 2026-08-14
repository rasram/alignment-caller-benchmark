#!/usr/bin/env python3
"""
PHASE 7b — generate config/conditions.tsv, the sweep design table.

DESIGN: one-factor-at-a-time (OFAT). Start from the baseline and vary exactly one
axis at a time, holding the others at baseline:

    coverage     5, 10, 20, 30*, 50, 100     -> 6 conditions
    read_length  75, 100, 150*               -> 3 conditions
    qs_shift     0*, -2, -5, -10             -> 4 conditions
                                                = 13, minus the baseline counted
                                                  three times = 11 unique

(* = baseline: coverage 30, read length 150, qs shift 0)

WHY OFAT AND NOT A FULL FACTORIAL
A full grid would be 6 x 3 x 4 = 72 conditions per genome, and with 5 seeds x 2
genomes x 9 pipelines that is 6,480 pipeline runs — far beyond the compute budget.
OFAT answers "how does each factor affect each pipeline, on its own" with 11.

The cost is real and should be stated in the report: OFAT cannot detect
INTERACTIONS. If a pipeline only degrades when coverage is low AND reads are
short, this design will miss it, because it never varies two axes at once. That
is a deliberate trade-off, and a candidate for Phase 8+ to revisit by adding a
small factorial patch around any interesting region the OFAT sweep reveals.

Usage: python3 scripts/make_conditions.py [--seeds 5] [--out config/conditions.tsv]
"""

import argparse
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GENOMES = ["phiX", "ecoli"]
BASE = dict(coverage=30, read_length=150, qs_shift=0)
AXES = dict(
    coverage=[5, 10, 20, 30, 50, 100],
    read_length=[75, 100, 150],
    qs_shift=[0, -2, -5, -10],
)


def conditions():
    """Baseline first, then each off-baseline level of each axis. De-duplicated."""
    seen = set()
    out = []
    base = tuple(BASE[k] for k in ("coverage", "read_length", "qs_shift"))
    seen.add(base)
    out.append(dict(BASE, axis="baseline"))
    for axis, levels in AXES.items():
        for lv in levels:
            cond = dict(BASE)
            cond[axis] = lv
            key = tuple(cond[k] for k in ("coverage", "read_length", "qs_shift"))
            if key in seen:
                continue
            seen.add(key)
            out.append(dict(cond, axis=axis))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(REPO, "config", "conditions.tsv"))
    args = ap.parse_args()

    conds = conditions()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    cols = ["condition_id", "genome", "coverage", "read_length", "qs_shift",
            "seed", "axis", "is_baseline", "tag"]
    n = 0
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in GENOMES:
            for i, c in enumerate(conds):
                for seed in range(1, args.seeds + 1):
                    is_base = (c["axis"] == "baseline")
                    tag = (f'{g}_cov{c["coverage"]}_len{c["read_length"]}'
                           f'_err{c["qs_shift"]}_seed{seed}')
                    fh.write("\t".join(str(x) for x in [
                        f"{g}_c{i:02d}", g, c["coverage"], c["read_length"],
                        c["qs_shift"], seed, c["axis"],
                        "yes" if is_base else "no", tag]) + "\n")
                    n += 1

    print(f"Wrote {args.out}")
    print(f"  unique conditions per genome : {len(conds)}")
    print(f"  genomes                      : {len(GENOMES)}")
    print(f"  seeds                        : {args.seeds}")
    print(f"  total rows (condition x seed): {n}")
    print(f"  pipeline runs (x 9)          : {n * 9}")
    print()
    for i, c in enumerate(conds):
        star = " *baseline" if c["axis"] == "baseline" else ""
        print(f"  c{i:02d}  cov={c['coverage']:<4} len={c['read_length']:<4} "
              f"qs={c['qs_shift']:<4} (varying: {c['axis']}){star}")


if __name__ == "__main__":
    main()
