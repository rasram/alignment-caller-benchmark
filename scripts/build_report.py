#!/usr/bin/env python3
"""
STEP 10 — assemble docs/FINAL_REPORT.md from a template plus GENERATED tables.

Every number table in the report is computed here from results/ — none is typed
by hand — so the report cannot drift from the data, and re-running the workflow
and then this script regenerates a consistent report.

The narrative lives in docs/report/FINAL_REPORT.template.md and refers to tables
as {{NAME}} placeholders. A placeholder with no generator, or a generated table
the template never uses, is reported, so nothing is silently dropped.

Usage: python3 scripts/build_report.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda *p: os.path.join(REPO, "results", *p)          # noqa: E731
ALN = {"bwa": "BWA-MEM", "bowtie2": "Bowtie2", "minimap2": "minimap2"}
CAL = {"gatk": "GATK", "freebayes": "FreeBayes", "bcftools": "BCFtools"}
ALIGNERS, CALLERS = list(ALN), list(CAL)
ORDER = ["coverage 5x", "coverage 10x", "coverage 20x", "baseline (30x, 150bp, qs0)",
         "coverage 50x", "coverage 100x", "read length 75bp", "read length 100bp",
         "qs shift -2", "qs shift -5", "qs shift -10"]


def md(df, floatfmt=None):
    """Plain GitHub-flavoured markdown table (no tabulate dependency)."""
    cols = list(df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float) and floatfmt and c in floatfmt:
                cells.append("—" if pd.isna(v) else floatfmt[c].format(v))
            elif isinstance(v, float) and pd.isna(v):
                cells.append("—")
            else:
                cells.append(str(v))
        rows.append(cells)
    # Pandoc sizes a wrapped pipe table's columns by the RELATIVE number of dashes
    # in the separator row, so "|---|---|" gives every column the same width and a
    # short "Condition" column steals space from long ones. Size each column by
    # its longest cell, but never below its longest unbreakable word (words cannot
    # wrap, so a narrower column overprints its neighbour) and capped, so one long
    # wrapping column cannot squeeze the rest.
    def word(t):
        return max((len(x) for x in str(t).split()), default=0)
    w = []
    for i, c in enumerate(cols):
        floor = max([word(c)] + [word(r[i]) for r in rows]) + 2
        w.append(min(max([len(r[i]) for r in rows] + [floor]), max(floor, 30)))
    out = ["| " + " | ".join(str(c) for c in cols) + " |",
           "|" + "|".join("-" * n for n in w) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def short_cond(c):
    """'coverage 5x' -> '5×', 'read length 75bp' -> '75 bp', 'qs shift -2' -> 'qs −2'."""
    if c.startswith("baseline"):
        return "baseline (30×)"
    c = re.sub(r"^coverage (\d+)x$", r"\1×", c)
    c = re.sub(r"^read length (\d+)bp$", r"\1 bp reads", c)
    return re.sub(r"^qs shift -(\d+)$", r"qs −\1", c)


def pname(p):
    a, c = p.split("+")
    return f"{ALN.get(a, a)} + {CAL.get(c, c)}"


def cond_sort(df, col="condition"):
    return df.assign(_o=df[col].map({c: i for i, c in enumerate(ORDER)})).sort_values("_o") \
             .drop(columns="_o")


# ---------------------------------------------------------------------------
def t_baseline(summ, vt):
    s = summ[(summ.genome == "ecoli") & (summ.axis == "baseline") & (summ.variant_type == vt)]
    rows = []
    for a in ALIGNERS:
        row = {"": f"**{ALN[a]}**"}
        for c in CALLERS:
            d = s[(s.aligner == a) & (s.caller == c)]
            row[CAL[c]] = (f"{d.f1_mean.iloc[0]:.4f} ± {d.f1_sd.iloc[0]:.4f}"
                           if len(d) else "—")
        rows.append(row)
    return md(pd.DataFrame(rows))


def t_effects(eff, vt):
    e = cond_sort(eff[(eff.genome == "ecoli") & (eff.variant_type == vt)].copy())
    out = pd.DataFrame({
        "Condition": e.condition.map(short_cond),
        "Best pipeline": e.best_pipeline.map(pname),
        "Aligner spread": e.aligner_spread, "Caller spread": e.caller_spread,
        "Seed sd": e.seed_sd,
        "Aligner share": e.aligner_share, "Caller share": e.caller_share,
        "q (aligner)": e.aligner_q, "q (caller)": e.caller_q,
    })
    f = {k: "{:.4f}" for k in ["Aligner spread", "Caller spread"]}
    f.update({"Seed sd": "{:.5f}", "Aligner share": "{:.2f}", "Caller share": "{:.2f}",
              "q (aligner)": "{:.1e}", "q (caller)": "{:.1e}"})
    return md(out, f)


def t_effect_summary(eff):
    rows = []
    for vt in ["SNV", "INDEL"]:
        e = eff[(eff.genome == "ecoli") & (eff.variant_type == vt) & eff.aligner_q.notna()]
        rows.append({
            "Variant type": vt, "Conditions tested": len(e),
            "Aligner effect significant": f"{int((e.aligner_q < 0.05).sum())}/{len(e)}",
            "Caller effect significant": f"{int((e.caller_q < 0.05).sum())}/{len(e)}",
            "Interaction significant": f"{int((e.interaction_q < 0.05).sum())}/{len(e)}",
            "Aligner share > caller share": f"{int((e.aligner_share > e.caller_share).sum())}/{len(e)}",
            "Median aligner share": f"{e.aligner_share.median():.2f}",
            "Median caller share": f"{e.caller_share.median():.2f}",
            "Friedman significant": f"{int((e.friedman_q < 0.05).sum())}/{len(e)}",
        })
    return md(pd.DataFrame(rows))


def t_align_metrics():
    am = pd.read_csv(R("align_metrics.tsv"), sep="\t")
    pa = pd.read_csv(R("placement_accuracy.tsv"), sep="\t")
    x = am.merge(pa[["tag", "aligner", "placement_accuracy"]], on=["tag", "aligner"])
    # NB: never name these columns "cov": `x.cov` is DataFrame.cov (the covariance
    # METHOD), so `x.cov == 5` silently compares a function to 5 and matches nothing
    # — which is exactly how this table first came out empty.
    p = x.tag.str.extract(r"^(?P<genome>[^_]+)_cov(?P<c_cov>\d+)_len(?P<c_len>\d+)_err(?P<c_qs>-?\d+)_seed")
    x = pd.concat([x, p], axis=1)
    x = x[x["genome"] == "ecoli"]
    for c in ["c_cov", "c_len", "c_qs"]:
        x[c] = x[c].astype(int)
    keep = [(5, 150, 0), (30, 150, 0), (100, 150, 0), (30, 75, 0), (30, 150, -10)]
    lab = {(5, 150, 0): "5×", (30, 150, 0): "30× (baseline)", (100, 150, 0): "100×",
           (30, 75, 0): "75 bp reads", (30, 150, -10): "qs −10"}
    rows = []
    for k in keep:
        for a in ALIGNERS:
            d = x[(x["c_cov"] == k[0]) & (x["c_len"] == k[1]) & (x["c_qs"] == k[2])
                  & (x["aligner"] == a)]
            if d.empty:
                continue
            rows.append({"Condition": lab[k], "Aligner": ALN[a],
                         "Placement ±10 bp": f"{100*d.placement_accuracy.mean():.3f}%",
                         "Mapping rate": f"{100*d.mapping_rate.mean():.3f}%",
                         "Mean MAPQ": f"{d.mean_mapq.mean():.1f}",
                         "MAPQ0 reads": f"{100*d.mapq0_frac.mean():.3f}%",
                         "Mean depth": f"{d.mean_depth.mean():.1f}"})
    return md(pd.DataFrame(rows))


def t_runtime():
    rt = pd.read_csv(R("analysis", "runtime_by_condition.tsv"), sep="\t")
    rt = rt[(rt.genome == "ecoli") & (rt.read_length == 150) & (rt.qs_shift == 0)]
    rows = []
    for cov in sorted(rt.coverage.unique()):
        d = rt[rt.coverage == cov]
        row = {"Coverage": f"{int(cov)}×"}
        for a in ALIGNERS:
            s = d[(d.stage == "align") & (d.aligner == a)].seconds
            row[f"{ALN[a]} (s)"] = f"{s.mean():.1f}" if len(s) else "—"
        for c in CALLERS:
            s = d[(d.stage == "call") & (d.caller == c)].seconds
            row[f"{CAL[c]} (s)"] = f"{s.mean():.1f}" if len(s) else "—"
        rows.append(row)
    mem = rt[rt.coverage == 30]
    note = ("\n\nPeak memory at 30× — "
            + ", ".join(f"{ALN[a]} {mem[(mem.stage == 'align') & (mem.aligner == a)].max_rss_mb.mean():.0f} MB"
                        for a in ALIGNERS)
            + "; " + ", ".join(f"{CAL[c]} {mem[(mem.stage == 'call') & (mem.caller == c)].max_rss_mb.mean():.0f} MB"
                              for c in CALLERS) + ".")
    return md(pd.DataFrame(rows)) + note


def t_model_validation():
    v = pd.read_csv(R("model", "validation.tsv"), sep="\t")
    name = {"pipeline_mean": "Pipeline mean (condition-blind)", "tree": "Decision tree (depth 4)",
            "forest": "Random forest", "always_train_best": "Always pick training-set best",
            "random_pick": "Random pick (context)"}
    v["Model"] = v.model.map(name).fillna(v.model)
    def seeds(m):
        n = [int(x) for x in re.findall(r"\d+", m.group(1))]
        return f"{n[0]}–{n[-1]}" if n == list(range(n[0], n[-1] + 1)) else ", ".join(map(str, n))
    v["Validation"] = v.validation.str.replace(r"\[(.*)\]", seeds, regex=True)
    # same row order in both validation blocks
    v = v.assign(_o=v.model.map({k: i for i, k in enumerate(name)})) \
         .sort_values(["validation", "_o"])    # "held-out" sorts before "leave-one"
    out = v[["Validation", "Model", "mae", "r2", "regret", "decisions"]].rename(
        columns={"mae": "MAE (F1)", "r2": "R²", "regret": "Mean regret (F1)",
                 "decisions": "Decisions"})
    return md(out, {"MAE (F1)": "{:.5f}", "R²": "{:.3f}", "Mean regret (F1)": "{:.5f}"})


def t_importance():
    imp = pd.read_csv(R("model", "importances.tsv"), sep="\t")
    g = imp.groupby("factor").importance_mean.sum().sort_values(ascending=False)
    lab = {"coverage": "Coverage", "is_indel": "Variant type (SNV/indel)", "caller": "Caller",
           "aligner": "Aligner", "mean_p": "Measured error rate", "read_length": "Read length"}
    return md(pd.DataFrame({"Factor": [lab.get(k, k) for k in g.index],
                            "Permutation importance (drop in R²)": g.values}),
              {"Permutation importance (drop in R²)": "{:.4f}"})


def describe_set(pipes):
    """Name a set of pipelines compactly. A tree leaf usually covers a full
    aligners x callers block, so describe it that way when it is one."""
    pipes = set(pipes.split(","))
    al = [a for a in ALIGNERS if any(p.startswith(a + "+") for p in pipes)]
    ca = [c for c in CALLERS if any(p.endswith("+" + c) for p in pipes)]
    if len(pipes) == 1:
        return pname(next(iter(pipes)))
    if len(pipes) == 9:
        return "all 9 tied (no preference)"
    if pipes == {f"{a}+{c}" for a in al for c in ca}:
        a_txt = "any aligner" if len(al) == 3 else "/".join(ALN[a] for a in al)
        c_txt = "any caller" if len(ca) == 3 else "/".join(CAL[c] for c in ca)
        return f"{a_txt} + {c_txt} ({len(pipes)} tied)"
    return "; ".join(pname(p) for p in sorted(pipes)) + f" ({len(pipes)} tied)"


def t_recommendations():
    r = pd.read_csv(R("model", "recommendations.tsv"), sep="\t")

    def lab(x):
        if x.coverage != 30:
            return f"{int(x.coverage)}×"
        if x.read_length != 150:
            return f"{int(x.read_length)} bp reads"
        return f"qs −{-int(x.qs_shift)}" if x.qs_shift else "baseline (30×)"
    out = pd.DataFrame({
        "Condition": r.apply(lab, axis=1), "Type": r.variant_type,
        "Actual best (mean F1)": r.actual_best.map(pname) + r.actual_best_f1.map(lambda v: f" ({v:.4f})"),
        "Tree's top-ranked set (mean F1)": r.model_top.map(describe_set)
            + r.model_top_mean_f1.map(lambda v: f" ({v:.4f})"),
        "Worst (mean F1)": r.worst.map(pname) + r.worst_f1.map(lambda v: f" ({v:.4f})"),
    })
    return md(out)


def t_titv():
    t = pd.read_csv(R("analysis", "titv_experiment.tsv"), sep="\t")
    t = t[t.genome == "ecoli"]
    rows = []
    for vt in ["SNV", "INDEL"]:
        for a in ALIGNERS:
            d = t[(t.variant_type == vt) & (t.aligner == a)]
            rows.append({"Type": vt, "Aligner": ALN[a],
                         **{CAL[c]: f"{d[d.caller == c].f1_titv05.iloc[0]:.4f} → "
                                    f"{d[d.caller == c].f1_titv2.iloc[0]:.4f}" for c in CALLERS}})
    return md(pd.DataFrame(rows))


def t_scoring():
    s = pd.read_csv(R("analysis", "scoring_method_effect.tsv"), sep="\t")
    s = s[s.genome == "ecoli"]
    rows = []
    for vt in ["SNV", "INDEL"]:
        for c in CALLERS:
            d = s[(s.variant_type == vt) & (s.caller == c)]
            rows.append({"Type": vt, "Caller": CAL[c], "Mean understatement": d["mean"].mean(),
                         "Largest (any condition)": d["max"].max()})
    return md(pd.DataFrame(rows), {"Mean understatement": "{:+.4f}", "Largest (any condition)": "{:+.4f}"})


def t_hardfilter():
    """Shown separately at 5x: there the DP >= 5 threshold sits AT mean depth and
    dominates any average, hiding how small the effect is everywhere else."""
    h = pd.read_csv(R("analysis", "hard_filter_effect.tsv"), sep="\t")
    h = h[h.genome == "ecoli"]
    rows = []
    for vt in ["SNV", "INDEL"]:
        for c in CALLERS:
            d = h[(h.variant_type == vt) & (h.caller == c)]
            low = d[d.condition == "coverage 5x"]
            rest = d[d.condition != "coverage 5x"]
            rows.append({"Type": vt, "Caller": CAL[c],
                         "ΔF1 at 5×": low.d_f1.mean(), "Δrecall at 5×": low.d_recall.mean(),
                         "ΔF1, other 10 conditions": rest.d_f1.mean(),
                         "Δprecision, other 10": rest.d_precision.mean(),
                         "Δrecall, other 10": rest.d_recall.mean()})
    f = {k: "{:+.4f}" for k in ["ΔF1 at 5×", "Δrecall at 5×", "ΔF1, other 10 conditions",
                                "Δprecision, other 10", "Δrecall, other 10"]}
    return md(pd.DataFrame(rows), f)


def t_phix(summ):
    s = summ[summ.genome == "phiX"]
    g = s.groupby("condition").agg(pipelines=("f1_mean", "size"),
                                   perfect=("f1_mean", lambda v: int((v >= 0.99999).sum())),
                                   min_f1=("f1_min", "min")).reset_index()
    g = cond_sort(g)
    g["Pipeline × type cells perfect on all 5 seeds"] = g.apply(
        lambda r: f"{r.perfect}/{r.pipelines}", axis=1)
    g["condition"] = g.condition.map(short_cond)
    return md(g[["condition", "Pipeline × type cells perfect on all 5 seeds", "min_f1"]].rename(
        columns={"condition": "Condition", "min_f1": "Lowest F1 (any seed)"}),
        {"Lowest F1 (any seed)": "{:.4f}"})


def t_fp_near_indel():
    f = pd.read_csv(R("analysis", "fp_near_indel.tsv"), sep="\t")
    # NB: bracket access — `f.pipe` is DataFrame.pipe (a METHOD), the same trap as `.cov`
    f["pipe"] = f.aligner + "+" + f.caller
    show = ["bowtie2+freebayes", "bowtie2+bcftools", "bwa+bcftools", "minimap2+bcftools"]
    rows = []
    for cov, g in f.groupby("coverage"):
        row = {"Coverage": f"{cov}×"}
        for p in show:
            d = g[g["pipe"] == p].iloc[0]
            row[pname(p)] = (f"{int(d.fp_snv)} ({100 * d.frac_near:.0f}%, {d.median_dist:g} bp)"
                             if d.fp_snv else "0")
        rest = g[~g["pipe"].isin(show)]
        row["Other 5 pipelines (total)"] = str(int(rest.fp_snv.sum()))
        rows.append(row)
    bg = f.background_frac.iloc[0]
    return (md(pd.DataFrame(rows)) + f"\n\nCells: false-positive SNVs (share within 150 bp of a "
            f"true indel, median distance to it). E. coli, seed 1. For comparison, "
            f"{100 * bg:.1f}% of the genome lies within 150 bp of a true indel.")


def t_fn_repeats():
    f = pd.read_csv(R("analysis", "fn_repeats.tsv"), sep="\t")
    rows = []
    for (a, c), g in f.groupby(["aligner", "caller"], sort=False):
        row = {"Pipeline": pname(f"{a}+{c}")}
        for cov in [30, 100]:
            d = g[g.coverage == cov].iloc[0]
            row[f"Missed at {cov}×"] = int(d.fn_sites)
            row[f"…of which in low-MAPQ sites ({cov}×)"] = f"{100 * d.frac_fn_low_mapq:.0f}%"
        rows.append(row)
    bg = f[f.coverage == 30].groupby("aligner").frac_truth_low_mapq.first()
    note = ", ".join(f"{ALN[a]} {100 * bg[a]:.1f}%" for a in ALIGNERS)
    return (md(pd.DataFrame(rows)) + "\n\nLow-MAPQ site: fewer than half of the reads covering "
            f"it have mapping quality ≥ 20. Share of ALL true variant sites that are low-MAPQ "
            f"(30×): {note}. E. coli, seed 1, all variant types.")


def t_phix_errors():
    f = pd.read_csv(R("analysis", "phix_errors.tsv"), sep="\t")
    out = pd.DataFrame({
        "Error": f.kind.map({"FN": "missed (FN)", "FP": "false call (FP)"}),
        "Position": f.pos, "Change": f.change.str.replace(">", " → ", regex=False),
        "Runs affected": f.runs.astype(str) + " / " + f.of_runs.astype(str),
        "Callers": f.callers.str.split(",").map(lambda cs: ", ".join(CAL[c] for c in cs)),
        "Seeds": f.seeds.astype(str),
        "Highest QUAL": f.max_qual.map(lambda v: "—" if pd.isna(v) else f"{v:.0f}"),
    })
    return md(out)


def values(summ, eff):
    """Scalar numbers quoted in the prose, generated like the tables so the text
    cannot drift from the data. Used as {{V_NAME}} in the template."""
    e = eff[eff.genome == "ecoli"]
    v = {}
    for vt in ["SNV", "INDEL"]:
        x = e[e.variant_type == vt]
        v[f"V_{vt}_ALN_SHARE"] = f"{x.aligner_share.median():.2f}"
        v[f"V_{vt}_CAL_SHARE"] = f"{x.caller_share.median():.2f}"
        v[f"V_{vt}_INT_SHARE"] = f"{x.interaction_share.median():.2f}"
        v[f"V_{vt}_SEED_SD"] = f"{x.seed_sd.median():.5f}"
        v[f"V_{vt}_ALN_SPREAD"] = f"{x.aligner_spread.median():.4f}"
        b = summ[(summ.genome == "ecoli") & (summ.axis == "baseline") & (summ.variant_type == vt)]
        v[f"V_{vt}_BASE_MIN"] = f"{b.f1_mean.min():.4f}"
        v[f"V_{vt}_BASE_MAX"] = f"{b.f1_mean.max():.4f}"
        c5 = summ[(summ.genome == "ecoli") & (summ.condition == "coverage 5x") & (summ.variant_type == vt)]
        v[f"V_{vt}_5X_MIN"] = f"{c5.f1_mean.min():.3f}"
        v[f"V_{vt}_5X_MAX"] = f"{c5.f1_mean.max():.3f}"
    val = pd.read_csv(R("model", "validation.tsv"), sep="\t")
    for lab, key in [("held-out", "HO"), ("leave-one", "LOCO")]:
        x = val[val.validation.str.startswith(lab)].set_index("model")
        for m, k in [("always_train_best", "ALWAYS"), ("forest", "FOREST"), ("tree", "TREE"),
                     ("random_pick", "RANDOM")]:
            v[f"V_{key}_{k}"] = f"{x.loc[m, 'regret']:.5f}"
    lo = pd.read_csv(R("model", "loco_by_condition.tsv"), sep="\t")
    l5 = lo[lo.held_out == "cov5 len150 qs0"].set_index("model")
    v["V_LOCO5_ALWAYS"] = f"{l5.loc['always_train_best', 'regret']:.4f}"
    v["V_LOCO5_TREE"] = f"{l5.loc['tree', 'regret']:.4f}"
    rec = pd.read_csv(R("model", "recommendations.tsv"), sep="\t")
    v["V_REC_BEST_IN_SET"] = f"{int(rec.best_in_model_top.sum())}/{len(rec)}"
    rm = pd.read_csv(R("read_metrics.tsv"), sep="\t")
    rm = rm[rm.tag.str.startswith("ecoli")]
    base = rm[rm.tag.str.contains("_cov30_len150_err0_")].mean_p.mean()
    worst = rm[rm.tag.str.contains("_err-10_")].mean_p.mean()
    v["V_ERR_BASE"] = f"{100 * base:.2f}%"
    v["V_ERR_WORST"] = f"{100 * worst:.2f}%"
    v["V_ERR_FOLD"] = f"{worst / base:.0f}"
    pe = pd.read_csv(R("analysis", "phix_errors.tsv"), sep="\t")
    v["V_PHIX_RUNS"] = str(int(pe.of_runs.iloc[0]))
    return v


def main():
    summ = pd.read_csv(R("analysis", "summary_by_condition.tsv"), sep="\t")
    eff = pd.read_csv(R("analysis", "effects.tsv"), sep="\t")
    tables = {
        "T_BASELINE_SNV": lambda: t_baseline(summ, "SNV"),
        "T_BASELINE_INDEL": lambda: t_baseline(summ, "INDEL"),
        "T_EFFECT_SUMMARY": lambda: t_effect_summary(eff),
        "T_EFFECTS_SNV": lambda: t_effects(eff, "SNV"),
        "T_EFFECTS_INDEL": lambda: t_effects(eff, "INDEL"),
        "T_ALIGN": t_align_metrics, "T_RUNTIME": t_runtime,
        "T_MODEL": t_model_validation, "T_IMPORTANCE": t_importance,
        "T_RECOMMEND": t_recommendations, "T_TITV": t_titv, "T_SCORING": t_scoring,
        "T_HARDFILTER": t_hardfilter, "T_PHIX": lambda: t_phix(summ),
        "T_FP_NEAR_INDEL": t_fp_near_indel, "T_FN_REPEATS": t_fn_repeats,
        "T_PHIX_ERRORS": t_phix_errors,
    }
    vals = values(summ, eff)
    tpl_path = os.path.join(REPO, "docs", "report", "FINAL_REPORT.template.md")
    tpl = open(tpl_path).read()
    used = set(re.findall(r"\{\{(T_[A-Z_]+)\}\}", tpl))
    used_v = set(re.findall(r"\{\{(V_[A-Z0-9_]+)\}\}", tpl))
    missing = (used - set(tables)) | (used_v - set(vals))
    if missing:
        sys.exit(f"template uses placeholders with no generator: {sorted(missing)}")
    for k in sorted(used):
        tpl = tpl.replace("{{" + k + "}}", tables[k]())
    for k in sorted(used_v):
        tpl = tpl.replace("{{" + k + "}}", vals[k])
    left = re.findall(r"\{\{[^}]*\}\}", tpl)
    if left:
        sys.exit(f"unresolved placeholders: {sorted(set(left))}")
    unused = set(tables) - used
    if unused:
        print(f"note: generated but unused tables: {sorted(unused)}", file=sys.stderr)
    out = os.path.join(REPO, "docs", "FINAL_REPORT.md")
    with open(out, "w") as fh:
        fh.write(tpl)
    print(f"wrote {os.path.relpath(out, REPO)} ({len(used)} generated tables, "
          f"{len(used_v)} generated values)")


if __name__ == "__main__":
    main()
