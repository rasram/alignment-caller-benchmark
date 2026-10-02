#!/usr/bin/env python3
"""
STEP 9 — figures for the final report.

Design rules (dataviz skill; palette validated with validate_palette.js):
  * colour = ALIGNER, using the reference palette's first three categorical
    slots — the only cap that passes the all-pairs CVD and normal-vision checks.
    The project deck's teal/coral/amber was tested and FAILED (teal below the
    chroma floor; coral vs amber dE 13.9 < 15), so it is not reused here.
  * panel = CALLER (small multiples), so no panel carries more than 3 series and
    identity never rests on colour alone.
  * 2px lines, >= 8px markers, solid hairline grid, recessive axes, one y-axis
    per panel, legend AND direct labels, text in ink — never in series colour.
  * the aqua slot sits below 3:1 contrast on the surface; the relief rule is met
    by direct labels plus a data table for every figure in the report.

Inputs:  results/analysis/*.tsv, results/runtime.tsv, results/model/tree.pkl
Outputs: results/figures/F*.png
"""
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(REPO, "results", "analysis")
OUT = os.path.join(REPO, "results", "figures")

# --- tokens ----------------------------------------------------------------
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
ALN_COLOR = {"bwa": "#2a78d6", "bowtie2": "#eb6834", "minimap2": "#1baf7a"}
ALN_NAME = {"bwa": "BWA-MEM", "bowtie2": "Bowtie2", "minimap2": "minimap2"}
CAL_NAME = {"gatk": "GATK", "freebayes": "FreeBayes", "bcftools": "BCFtools"}
CAL_MARK = {"gatk": "o", "freebayes": "s", "bcftools": "^"}
ALIGNERS, CALLERS = ["bwa", "bowtie2", "minimap2"], ["gatk", "freebayes", "bcftools"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9.5, "text.color": INK, "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "xtick.major.size": 0, "ytick.major.size": 0,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "legend.frameon": False, "legend.fontsize": 9,
})
LW, MS = 1.6, 6.5          # ~2px lines, ~8.5px markers at 100 dpi


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", os.path.relpath(p, REPO))


def legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], color=ALN_COLOR[a], lw=LW, marker="o", ms=MS,
                   markeredgecolor=SURFACE, label=ALN_NAME[a]) for a in ALIGNERS]


def direct_labels(ax, items, min_gap_frac=0.055):
    """Labels at the lines' right ends, in INK (the line carries the colour).

    Lines that converge would put their labels on top of each other, so labels
    are sorted by y and pushed apart to a minimum gap (a fraction of the axis
    height); a thin leader joins each moved label to its line end.
    """
    if not items:
        return
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * min_gap_frac
    items = sorted(items, key=lambda it: it[1])
    ys = [it[1] for it in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    over = ys[-1] - (hi - gap * 0.3)
    if over > 0:
        ys = [y - over for y in ys]
    # Every label is an offset in POINTS from its line end: same 10 pt rightward
    # gap whether or not it was moved, so a moved label never lands on a marker.
    fig = ax.figure
    ax_h_pt = ax.get_window_extent().height * 72.0 / fig.dpi
    for (x, y, text), ly in zip(items, ys):
        dy = (ly - y) / (hi - lo) * ax_h_pt
        moved = abs(dy) > 0.5
        ax.annotate(text, (x, y), xytext=(12, dy), textcoords="offset points",
                    va="center", fontsize=8, color=INK2, annotation_clip=False,
                    arrowprops=(dict(arrowstyle="-", color=AXIS, lw=0.6, shrinkA=0,
                                     shrinkB=4) if moved else None))


def small_multiples(summ, axis, xcol, xlabel, title, fname, logx=False, xticks=None):
    s = summ[(summ.genome == "ecoli") &
             ((summ.axis == axis) | (summ.axis == "baseline"))].copy()
    if s.empty:
        return
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True)
    for i, vt in enumerate(["SNV", "INDEL"]):
        for j, ca in enumerate(CALLERS):
            ax = axes[i, j]
            ends = []
            for al in ALIGNERS:
                d = s[(s.variant_type == vt) & (s.caller == ca) & (s.aligner == al)] \
                    .sort_values(xcol)
                if d.empty:
                    continue
                ax.errorbar(d[xcol], d.f1_mean, yerr=d.f1_sd.fillna(0), color=ALN_COLOR[al],
                            lw=LW, marker="o", ms=MS, markeredgecolor=SURFACE,
                            markeredgewidth=1.2, elinewidth=0.9, capsize=0, zorder=3)
                ends.append((d[xcol].values[-1], d.f1_mean.values[-1], ALN_NAME[al]))
            if logx:
                ax.set_xscale("log")
            if xticks is not None:
                ax.set_xticks(xticks)
                ax.set_xticklabels([str(x) for x in xticks])
                ax.minorticks_off()
            if j == 2:
                direct_labels(ax, ends)
            if i == 0:
                ax.set_title(CAL_NAME[ca])
            if j == 0:
                ax.set_ylabel(f"{vt} F1")
            if i == 1:
                ax.set_xlabel(xlabel)
    fig.legend(handles=legend_handles(), loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.03), title=None)
    fig.suptitle(title, y=1.08, fontsize=11.5, fontweight="bold", color=INK)
    fig.text(0.5, -0.02, "E. coli · mean of 5 seeds · bars = ±1 sd across seeds",
             ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    save(fig, fname)


def heatmaps(summ):
    covs = [5, 30, 100]
    s = summ[(summ.genome == "ecoli") & (summ.read_length == 150) & (summ.qs_shift == 0)
             & (summ.coverage.isin(covs))]
    if s.empty:
        return
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
    for i, vt in enumerate(["SNV", "INDEL"]):
        sv = s[s.variant_type == vt]
        lo, hi = sv.f1_mean.min(), sv.f1_mean.max()   # shared scale across coverages
        for j, cov in enumerate(covs):
            ax = axes[i, j]
            ax.grid(False)
            m = np.full((3, 3), np.nan)
            for a, al in enumerate(ALIGNERS):
                for c, ca in enumerate(CALLERS):
                    d = sv[(sv.coverage == cov) & (sv.aligner == al) & (sv.caller == ca)]
                    if len(d):
                        m[a, c] = d.f1_mean.iloc[0]
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list("seq", SEQ)
            ax.imshow(m, cmap=cmap, vmin=lo, vmax=hi, aspect="equal")
            for a in range(3):
                for c in range(3):
                    if np.isnan(m[a, c]):
                        continue
                    t = (m[a, c] - lo) / (hi - lo) if hi > lo else 1
                    ax.text(c, a, f"{m[a, c]:.4f}", ha="center", va="center", fontsize=8.5,
                            color="#ffffff" if t > 0.55 else INK)
            ax.set_xticks(range(3))
            ax.set_xticklabels([CAL_NAME[c] for c in CALLERS], fontsize=8.5)
            ax.set_yticks(range(3))
            ax.set_yticklabels([ALN_NAME[a] for a in ALIGNERS] if j == 0 else [], fontsize=8.5)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if i == 0:
                ax.set_title(f"{cov}× coverage")
            if j == 0:
                ax.set_ylabel(vt, fontsize=10, fontweight="bold", color=INK)
            # white gaps between cells (2px surface spacer)
            for k in range(4):
                ax.axhline(k - 0.5, color=SURFACE, lw=2.5)
                ax.axvline(k - 0.5, color=SURFACE, lw=2.5)
    fig.suptitle("F1 by pipeline at low, baseline and high coverage", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.01, "E. coli, 150 bp, qs 0 · mean of 5 seeds · darker = higher F1 · "
             "one colour scale per row", ha="center", fontsize=8.5, color=MUTED)
    fig.tight_layout()
    save(fig, "F4_heatmaps_by_coverage.png")


ORDER = ["coverage 5x", "coverage 10x", "coverage 20x", "baseline (30x, 150bp, qs0)",
         "coverage 50x", "coverage 100x", "read length 75bp", "read length 100bp",
         "qs shift -2", "qs shift -5", "qs shift -10"]
SHORT = {"baseline (30x, 150bp, qs0)": "30× (base)", "read length 75bp": "75 bp",
         "read length 100bp": "100 bp", "qs shift -2": "qs −2", "qs shift -5": "qs −5",
         "qs shift -10": "qs −10"}


def variance_shares(eff):
    e = eff[(eff.genome == "ecoli")].copy()
    if "aligner_share" not in e or e.aligner_share.isna().all():
        return
    e["order"] = e.condition.map({c: i for i, c in enumerate(ORDER)})
    e = e.sort_values("order")
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.4), sharex=True)
    parts = [("aligner_share", "Aligner", ALN_COLOR["bwa"]),
             ("caller_share", "Caller", ALN_COLOR["bowtie2"]),
             ("interaction_share", "Aligner × caller", "#c3c2b7")]
    for i, vt in enumerate(["SNV", "INDEL"]):
        ax = axes[i]
        d = e[e.variant_type == vt].reset_index(drop=True)
        x = np.arange(len(d))
        bottom = np.zeros(len(d))
        for col, lab, colr in parts:
            v = d[col].fillna(0).values
            ax.bar(x, v, bottom=bottom, color=colr, width=0.72, label=lab,
                   edgecolor=SURFACE, linewidth=1.5, zorder=3)
            bottom += v
        # Significance marker. A condition with fewer than 2 seeds was NOT TESTED —
        # labelling it "n.s." would misreport an untested condition as a negative.
        for k, r in d.iterrows():
            if pd.isna(r.get("aligner_q")):
                lab = "not tested"
            else:
                tag = [x for x, q in (("A", r.aligner_q), ("C", r.caller_q)) if q < 0.05]
                lab = "·".join(tag) if tag else "n.s."
            ax.text(k, 1.03, lab, ha="center", fontsize=7.5, color=INK2)
        ax.set_ylim(0, 1.12)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylabel(f"{vt}\nshare of between-\npipeline variance")
        # separators between the coverage / read-length / error groups, from the data
        grp = d.condition.map(lambda c: "len" if c.startswith("read") else
                              "qs" if c.startswith("qs") else "cov").values
        for k in range(1, len(grp)):
            if grp[k] != grp[k - 1]:
                ax.axvline(k - 0.5, color=AXIS, lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT.get(c, c.replace("coverage ", "").replace("x", "×"))
                            for c in d.condition],
                           fontsize=8.5)
    axes[0].legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.08))
    fig.suptitle("Which choice drives the difference — aligner or caller?", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.0)
    fig.text(0.5, -0.02, "Type-II sums of squares from F1 ~ aligner × caller + seed, per "
             "condition (E. coli, 5 seeds). A / C = aligner / caller effect significant at "
             "BH q < 0.05; n.s. = tested, neither significant.", ha="center", fontsize=8,
             color=MUTED, wrap=True)
    fig.tight_layout()
    save(fig, "F7_variance_shares.png")


def runtime_vs_accuracy(summ):
    rp = os.path.join(A, "runtime_by_condition.tsv")
    if not os.path.exists(rp):
        return
    rt = pd.read_csv(rp, sep="\t")
    rt = rt[(rt.genome == "ecoli") & (rt.coverage == 30) & (rt.read_length == 150)
            & (rt.qs_shift == 0) & (rt.seed == 1)]
    if rt.empty:
        return
    al = rt[rt.stage == "align"].set_index("aligner").seconds
    ca = rt[rt.stage == "call"].set_index(["aligner", "caller"]).seconds
    base = summ[(summ.genome == "ecoli") & (summ.axis == "baseline")]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharex=True)
    for i, vt in enumerate(["SNV", "INDEL"]):
        ax = axes[i]
        pts = []
        for a in ALIGNERS:
            for c in CALLERS:
                d = base[(base.variant_type == vt) & (base.aligner == a) & (base.caller == c)]
                if d.empty or (a, c) not in ca.index:
                    continue
                t = al[a] + ca[(a, c)]
                pts.append((t, d.f1_mean.iloc[0], f"{ALN_NAME[a]}+{CAL_NAME[c]}"))
                ax.scatter(t, d.f1_mean.iloc[0], s=70, color=ALN_COLOR[a], marker=CAL_MARK[c],
                           edgecolor=SURFACE, linewidth=1.5, zorder=3)
        # Label only the Pareto front — pipelines no other pipeline beats on BOTH
        # speed and accuracy. Labelling all nine piled names on top of each other.
        front = [p for p in pts if not any(q[0] <= p[0] and q[1] >= p[1] and q != p
                                           and (q[0] < p[0] or q[1] > p[1]) for q in pts)]
        # Higher point labels upward, lower point downward: alternating by x-order
        # instead made the leader lines cross and the labels land on each other.
        front.sort(key=lambda p: -p[1])
        for k, (x, y, name) in enumerate(front):
            ax.annotate(name, (x, y), xytext=(16, 16 if k % 2 == 0 else -22),
                        textcoords="offset points", fontsize=7.5, color=INK2,
                        arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.6))
        lo, hi = ax.get_ylim()                   # headroom so upward labels never clip
        ax.set_ylim(lo, hi + (hi - lo) * 0.12)
        ax.set_title(vt)
        ax.set_xlabel("align + call wall time (s), timed with the machine to itself")
        if i == 0:
            ax.set_ylabel("F1 (mean of 5 seeds)")
    from matplotlib.lines import Line2D
    h = [Line2D([0], [0], color=ALN_COLOR[a], lw=0, marker="o", ms=MS, label=ALN_NAME[a])
         for a in ALIGNERS] + \
        [Line2D([0], [0], color=INK2, lw=0, marker=CAL_MARK[c], ms=MS, label=CAL_NAME[c])
         for c in CALLERS]
    fig.legend(handles=h, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.07))
    fig.suptitle("Speed vs accuracy at baseline (E. coli, 30×)", fontsize=11.5,
                 fontweight="bold", color=INK, y=1.14)
    fig.text(0.5, -0.04, "Labelled points are Pareto-optimal: no other pipeline is both faster "
             "and more accurate. Timing from seed 1, each tool run with the machine to itself; "
             "F1 = mean of 5 seeds.", ha="center", fontsize=8, color=MUTED)
    fig.tight_layout()
    save(fig, "F6_runtime_vs_accuracy.png")


def _question(feature, thr):
    """Turn a sklearn split into a plain-language yes/no question.
    Returns (question, label for the <= branch, label for the > branch)."""
    if feature == "coverage":
        return f"coverage ≤ {thr:g}×?", "yes", "no"
    if feature == "read_length":
        return f"reads ≤ {thr:g} bp?", "yes", "no"
    if feature == "mean_p":
        return f"error rate ≤ {100 * thr:.2f}%?", "yes", "no"
    if feature == "is_indel":
        return "indel?", "no (SNV)", "yes"
    if feature.startswith("aln_"):
        return f"aligner is {ALN_NAME[feature[4:]]}?", "no", "yes"
    if feature.startswith("cal_"):
        return f"caller is {CAL_NAME[feature[4:]]}?", "no", "yes"
    return f"{feature} ≤ {thr:.3g}?", "yes", "no"


def _cond_label(cov, rl, qs):
    if rl != 150:
        return f"the {rl:g} bp reads"
    if qs != 0:
        return f"qs {qs:g}"
    return "baseline" if cov == 30 else f"{cov:g}×"


def _split_caveats(model, feats):
    """Say what a measured-error split ACTUALLY separates.

    The sweep is one-factor-at-a-time, so measured error is not an independent
    axis: shorter reads also have a lower mean error (ART's error rises along the
    read), and within one condition it differs slightly from seed to seed. A tree
    split on it can therefore mean "75 bp reads" or merely "seed 3 vs the rest".
    Both are worked out from the training rows, not asserted."""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from fit_model import load
    df = load()
    paths = model.decision_path(df[feats]).tocsc()
    tr, out = model.tree_, {}
    for n in range(tr.node_count):
        if tr.children_left[n] == -1 or feats[tr.feature[n]] != "mean_p":
            continue
        rows = df.iloc[paths[:, n].nonzero()[0]]
        left = rows.mean_p <= tr.threshold[n]
        key = lambda r: tuple(r[["coverage", "read_length", "qs_shift"]])   # noqa: E731
        lc = {key(r) for _, r in rows[left].iterrows()}
        rc = {key(r) for _, r in rows[~left].iterrows()}
        if lc & rc:
            out[n] = "splits one condition's seeds: noise"
        elif len(lc) == 1:
            out[n] = "i.e. " + _cond_label(*next(iter(lc)))
        elif len(rc) == 1:
            out[n] = "no = " + _cond_label(*next(iter(rc)))
    return out


def tree_figure(max_depth=3):
    """Readable rendering of the fitted tree's top levels.

    sklearn's plot_tree printed machine-speak ("cal_gatk <= 0.5") and let boxes
    collide; the point of choosing a tree is that its rules can be read and
    defended, so it is drawn here with plain-language questions, yes/no edges and
    a layout derived from the tree itself (leaves spaced evenly by in-order index).
    """
    p = os.path.join(REPO, "results", "model", "tree.pkl")
    if not os.path.exists(p):
        return
    obj = pickle.load(open(p, "rb"))
    tr, feats = obj["tree"].tree_, obj["features"]
    caveat = _split_caveats(obj["tree"], feats)

    pos, order = {}, []
    def walk(n, d):
        leaf = tr.children_left[n] == -1 or d == max_depth
        if leaf:
            order.append(n); pos[n] = (len(order) - 1, d); return
        walk(tr.children_left[n], d + 1); walk(tr.children_right[n], d + 1)
        pos[n] = ((pos[tr.children_left[n]][0] + pos[tr.children_right[n]][0]) / 2, d)
    walk(0, 0)

    nleaf = len(order)
    fig, ax = plt.subplots(figsize=(max(11, 1.55 * nleaf), 6.2))
    ax.set_axis_off()
    ax.grid(False)
    X = lambda i: i * 1.0                       # noqa: E731
    Y = lambda d: -d * 1.35                     # noqa: E731
    vals = [tr.value[n][0][0] for n in order]
    vlo, vhi = min(vals), max(vals)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ[:5])

    def draw(n):
        x, d = pos[n]
        terminal = n in order
        v = tr.value[n][0][0]
        if not terminal:
            q, lab_l, lab_r = _question(feats[tr.feature[n]], tr.threshold[n])
            if n in caveat:
                q += f"\n({caveat[n]})"
            for child, lab in ((tr.children_left[n], lab_l), (tr.children_right[n], lab_r)):
                cx, cd = pos[child]
                ax.plot([X(x), X(cx)], [Y(d) - 0.22, Y(cd) + 0.24], color=AXIS, lw=1.0, zorder=1)
                ax.text((X(x) + X(cx)) / 2, (Y(d) + Y(cd)) / 2 - 0.02, lab, ha="center",
                        va="center", fontsize=8, color=INK2,
                        bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none"), zorder=2)
                draw(child)
            ax.text(X(x), Y(d), f"{q}\n{tr.n_node_samples[n]} results · mean F1 {v:.4f}",
                    ha="center", va="center", fontsize=8.5, color=INK, zorder=3,
                    bbox=dict(boxstyle="round,pad=0.45", fc="#ffffff", ec=AXIS, lw=0.8))
        else:
            t = (v - vlo) / (vhi - vlo) if vhi > vlo else 1
            more = "" if tr.children_left[n] == -1 else "\n(split further)"
            ax.text(X(x), Y(d), f"F1 {v:.4f}\n{tr.n_node_samples[n]} results{more}",
                    ha="center", va="center", fontsize=8.5, zorder=3,
                    color="#ffffff" if t > 0.6 else INK,
                    bbox=dict(boxstyle="round,pad=0.45", fc=cmap(t), ec="none"))

    draw(0)
    ax.set_xlim(-0.7, nleaf - 0.3)
    ax.set_ylim(Y(max_depth) - 0.6, 0.6)
    fig.suptitle("What predicts F1? — top of the fitted decision tree", fontsize=11.5,
                 fontweight="bold", color=INK, y=0.99)
    fig.text(0.5, 0.01, f"Regression tree (depth 4) on 990 E. coli results (495 pipeline runs × SNV/indel); "
             f"top {max_depth} levels shown. Terminal boxes: mean predicted F1, darker = higher.", ha="center",
             fontsize=8, color=MUTED)
    save(fig, "F5_decision_tree.png")


def main():
    summ = pd.read_csv(os.path.join(A, "summary_by_condition.tsv"), sep="\t")
    eff = pd.read_csv(os.path.join(A, "effects.tsv"), sep="\t")
    small_multiples(summ, "coverage", "coverage", "coverage (×, log scale)",
                    "F1 vs sequencing depth", "F1_f1_vs_coverage.png", logx=True,
                    xticks=[5, 10, 20, 30, 50, 100])
    small_multiples(summ, "read_length", "read_length", "read length (bp)",
                    "F1 vs read length", "F2_f1_vs_read_length.png", xticks=[75, 100, 150])
    # error axis: use the MEASURED error rate, joined from read metrics
    rm = pd.read_csv(os.path.join(REPO, "results", "read_metrics.tsv"), sep="\t")
    pe = (rm[rm.genome == "ecoli"].groupby(["coverage", "read_length", "qs_shift"])
          .mean_p.mean().rename("mean_p").reset_index())
    s2 = summ.merge(pe, on=["coverage", "read_length", "qs_shift"], how="left")
    s2["mean_p_pct"] = s2.mean_p * 100
    small_multiples(s2, "qs_shift", "mean_p_pct",
                    "measured error rate (%)",
                    "F1 vs sequencing error rate (ART quality shift 0, −2, −5, −10)",
                    "F3_f1_vs_error_rate.png")
    heatmaps(summ)
    tree_figure()
    runtime_vs_accuracy(summ)
    variance_shares(eff)


if __name__ == "__main__":
    main()
