#!/usr/bin/env python3
"""
STEP 7 — Are the differences between pipelines real?

The mid-semester headline ("the aligner matters more than the caller") rested on
ONE seed at ONE condition, with a spread of 0.0047 in SNV F1. This script asks
whether that survives contact with seed-to-seed noise, at every condition.

STATISTICAL DESIGN
------------------
Within a condition, all nine pipelines analyse the SAME reads for a given seed:
seed 3's FASTQs go to all three aligners, and each aligner's BAM to all three
callers. So seed is a BLOCK, not independent noise — this is a randomised
complete block design with a 3x3 factorial treatment (aligner x caller).

  Primary   two-way ANOVA   F1 ~ aligner * caller + seed      (type II SS)
            -> F-tests for aligner, caller, interaction against the residual,
               i.e. against seed-to-seed noise
            -> partial eta^2 per effect, and each effect's SHARE of the
               between-pipeline sum of squares (the cleanest answer to
               "does the aligner matter more than the caller?")
  Check     Friedman test   9 pipelines, seeds as blocks — rank-based, so it does
                            not lean on normality. F1 near 1.0 is bounded and
                            skewed, so the parametric result is cross-checked.
  Multiple  Benjamini-Hochberg FDR within each (genome, variant type) family,
  testing   across its 11 conditions, separately per effect.

Outputs -> results/analysis/
  summary_by_condition.tsv   mean / sd over seeds, every pipeline x condition
  effects.tsv                per condition: spreads, seed noise, ANOVA, Friedman
  titv_experiment.tsv        Ti/Tv 0.5 vs 2.0 at baseline (NOTES Phase 8)
  scoring_method_effect.tsv  single_run minus pre_split F1, by condition & caller
  runtime_by_condition.tsv   CLEAN timings only (timing_exclusive == yes)
  hard_filter_effect.tsv     filtered minus raw F1/precision/recall, by caller
  headline.txt               the numbers the report quotes

Usage: python3 scripts/analyse_sweep.py
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", "analysis")
ALIGNERS = ["bwa", "bowtie2", "minimap2"]
CALLERS = ["gatk", "freebayes", "bcftools"]
KEYS = ["genome", "variant_type", "axis", "coverage", "read_length", "qs_shift"]


def cond_label(r):
    if r["axis"] == "baseline":
        return "baseline (30x, 150bp, qs0)"
    if r["axis"] == "coverage":
        return f"coverage {r['coverage']}x"
    if r["axis"] == "read_length":
        return f"read length {r['read_length']}bp"
    return f"qs shift {r['qs_shift']}"


def load():
    res = pd.read_csv(os.path.join(REPO, "results", "results.tsv"), sep="\t")
    cond = pd.read_csv(os.path.join(REPO, "config", "conditions.tsv"), sep="\t")
    axis = (cond[["genome", "coverage", "read_length", "qs_shift", "axis"]]
            .drop_duplicates())
    res = res.merge(axis, on=["genome", "coverage", "read_length", "qs_shift"], how="left")
    if res["axis"].isna().any():
        sys.exit("results.tsv contains a condition not in config/conditions.tsv")
    for c in ["TP", "FP", "FN", "precision", "recall", "f1", "seed"]:
        res[c] = pd.to_numeric(res[c], errors="coerce")
    res["condition"] = res.apply(cond_label, axis=1)
    return res


def summarise(single):
    g = single.groupby(KEYS + ["condition", "aligner", "caller"], sort=False)
    out = g.agg(n_seeds=("f1", "size"),
                f1_mean=("f1", "mean"), f1_sd=("f1", "std"),
                f1_min=("f1", "min"), f1_max=("f1", "max"),
                recall_mean=("recall", "mean"), precision_mean=("precision", "mean"),
                TP_mean=("TP", "mean"), FP_mean=("FP", "mean"), FN_mean=("FN", "mean"))
    return out.reset_index()


def anova(sub):
    """Two-way blocked ANOVA. Returns {} if there is nothing to test."""
    if sub["seed"].nunique() < 2 or sub["f1"].std() == 0:
        return {}
    m = smf.ols("f1 ~ C(aligner) * C(caller) + C(seed)", data=sub).fit()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tab = sm.stats.anova_lm(m, typ=2)
    ss_res = tab.loc["Residual", "sum_sq"]
    rows = {"aligner": "C(aligner)", "caller": "C(caller)",
            "interaction": "C(aligner):C(caller)"}
    ss_between = sum(tab.loc[v, "sum_sq"] for v in rows.values())
    out = {"resid_df": tab.loc["Residual", "df"]}
    for k, v in rows.items():
        ss = tab.loc[v, "sum_sq"]
        out[f"{k}_F"] = tab.loc[v, "F"]
        out[f"{k}_p"] = tab.loc[v, "PR(>F)"]
        out[f"{k}_eta2p"] = ss / (ss + ss_res) if (ss + ss_res) > 0 else np.nan
        out[f"{k}_share"] = ss / ss_between if ss_between > 0 else np.nan
    return out


def friedman(sub):
    piv = sub.pivot_table(index="seed", columns=["aligner", "caller"], values="f1")
    if piv.shape[0] < 2 or piv.isna().any().any():
        return {}
    try:
        stat, p = stats.friedmanchisquare(*[piv[c].values for c in piv.columns])
    except ValueError:              # e.g. every value tied
        return {}
    return {"friedman_chi2": stat, "friedman_p": p}


def effects(single):
    out = []
    for key, sub in single.groupby(KEYS + ["condition"], sort=False):
        row = dict(zip(KEYS + ["condition"], key))
        row["n_seeds"] = sub["seed"].nunique()
        al = sub.groupby("aligner")["f1"].mean()
        ca = sub.groupby("caller")["f1"].mean()
        row["aligner_spread"] = al.max() - al.min()
        row["caller_spread"] = ca.max() - ca.min()
        row["best_aligner"], row["worst_aligner"] = al.idxmax(), al.idxmin()
        row["best_caller"] = ca.idxmax()
        pipe = sub.groupby(["aligner", "caller"])["f1"]
        var = pipe.var(ddof=1)
        row["seed_sd"] = float(np.sqrt(var.mean())) if row["n_seeds"] > 1 else np.nan
        means = pipe.mean().sort_values(ascending=False)
        (a1, c1), (a2, c2) = means.index[0], means.index[1]
        row["best_pipeline"] = f"{a1}+{c1}"
        row["runner_up"] = f"{a2}+{c2}"
        row["best_minus_runner_up"] = means.iloc[0] - means.iloc[1]
        row["f1_range"] = means.iloc[0] - means.iloc[-1]
        row["worst_pipeline"] = "+".join(means.index[-1])
        row.update(anova(sub))
        row.update(friedman(sub))
        out.append(row)
    eff = pd.DataFrame(out)
    # BH-FDR within each (genome, variant type) family, per effect
    for col in ["aligner_p", "caller_p", "interaction_p", "friedman_p"]:
        if col not in eff:
            continue
        eff[col.replace("_p", "_q")] = np.nan
        for _, idx in eff.groupby(["genome", "variant_type"]).groups.items():
            ps = eff.loc[idx, col]
            ok = ps.notna()
            if ok.sum():
                eff.loc[ps[ok].index, col.replace("_p", "_q")] = \
                    multipletests(ps[ok], method="fdr_bh")[1]
    return eff


def titv(single):
    arc = os.path.join(REPO, "results", "archive", "titv0.5_baseline", "results.tsv")
    if not os.path.exists(arc):
        return pd.DataFrame()
    old = pd.read_csv(arc, sep="\t")
    old = old[old.scoring_method == "single_run"]
    new = single[(single.axis == "baseline") & (single.seed == 1)]
    k = ["genome", "aligner", "caller", "variant_type"]
    m = new[k + ["f1", "TP", "FP", "FN"]].merge(
        old[k + ["f1", "TP", "FP", "FN"]], on=k, suffixes=("_titv2", "_titv05"))
    m["delta_f1"] = m["f1_titv2"] - m["f1_titv05"]
    return m


def scoring_method(res):
    k = ["genome", "variant_type", "condition", "seed", "aligner", "caller"]
    raw = res[res.callset == "raw"]
    a = raw[raw.scoring_method == "single_run"][k + ["f1"]]
    b = raw[raw.scoring_method == "pre_split"][k + ["f1"]]
    m = a.merge(b, on=k, suffixes=("_single", "_presplit"))
    m["understatement"] = m["f1_single"] - m["f1_presplit"]
    return (m.groupby(["genome", "variant_type", "condition", "caller"])["understatement"]
            .agg(["mean", "max", "size"]).reset_index())


def hard_filter_effect(res):
    """Does the shared hard filter (QUAL>=20 && DP>=5) help or hurt each caller?

    Same pipeline, same run, single-run scoring: filtered F1 minus raw F1. Because
    QUAL is not calibrated identically across callers (NOTES 6.7), one threshold can
    be lenient for one caller and harsh for another — this measures that directly.
    """
    k = ["genome", "variant_type", "condition", "seed", "aligner", "caller"]
    s = res[res.scoring_method == "single_run"]
    a = s[s.callset == "raw"][k + ["f1", "precision", "recall"]]
    b = s[s.callset == "filt"][k + ["f1", "precision", "recall"]]
    m = a.merge(b, on=k, suffixes=("_raw", "_filt"))
    if m.empty:
        return m
    for c in ["f1", "precision", "recall"]:
        m[f"d_{c}"] = m[f"{c}_filt"] - m[f"{c}_raw"]
    return (m.groupby(["genome", "variant_type", "condition", "caller"])
            [["d_f1", "d_precision", "d_recall"]].mean().reset_index())


def runtime(res):
    rt_path = os.path.join(REPO, "results", "runtime.tsv")
    if not os.path.exists(rt_path):
        return pd.DataFrame()
    rt = pd.read_csv(rt_path, sep="\t")
    rt = rt[rt.exclusive == "yes"].copy()
    parts = rt["tag"].str.extract(
        r"^(?P<genome>[^_]+)_cov(?P<coverage>\d+)_len(?P<read_length>\d+)_err(?P<qs_shift>-?\d+)_seed(?P<seed>\d+)$")
    for c in ["coverage", "read_length", "qs_shift", "seed"]:   # genome stays a string
        parts[c] = pd.to_numeric(parts[c])
    rt = pd.concat([rt, parts], axis=1)
    rt["caller"] = rt["caller"].fillna("")
    return rt


def headline(eff, titv_df, sm_df, rt):
    L = []
    w = L.append
    ec = eff[eff.genome == "ecoli"]
    for vt in ["SNV", "INDEL"]:
        e = ec[ec.variant_type == vt]
        if e.empty:
            continue
        w(f"== E. coli {vt}: {len(e)} conditions ==")
        if "aligner_q" in e:
            sig_a = int((e.aligner_q < 0.05).sum())
            sig_c = int((e.caller_q < 0.05).sum())
            tested = int(e.aligner_q.notna().sum())
            w(f"  aligner effect significant (BH q<0.05): {sig_a}/{tested} conditions")
            w(f"  caller  effect significant (BH q<0.05): {sig_c}/{tested} conditions")
            w(f"  median share of between-pipeline SS — aligner {e.aligner_share.median():.2f}, "
              f"caller {e.caller_share.median():.2f}, interaction {e.interaction_share.median():.2f}")
            bigger = int((e.aligner_share > e.caller_share).sum())
            w(f"  aligner explains more than caller in {bigger}/{tested} conditions")
        w(f"  median seed sd of F1: {e.seed_sd.median():.5f}")
        for _, r in e.iterrows():
            q = f"q_aln={r.get('aligner_q', np.nan):.2g} q_cal={r.get('caller_q', np.nan):.2g}"
            w(f"    {r['condition']:<28} best {r['best_pipeline']:<20} range {r['f1_range']:.4f}  "
              f"aln {r['aligner_spread']:.4f} cal {r['caller_spread']:.4f} seed_sd {r['seed_sd']:.5f}  {q}")
    if not titv_df.empty:
        e = titv_df[titv_df.genome == "ecoli"]
        w("== Ti/Tv 0.5 -> 2.0 at baseline (E. coli) ==")
        for al in ALIGNERS:
            d = e[e.aligner == al].delta_f1.abs().max()
            w(f"  max |dF1| {al:9s} {d:.4f}")
    if not sm_df.empty:
        e = sm_df[sm_df.genome == "ecoli"]
        w("== Pre-split understatement of F1 (single_run minus pre_split), E. coli ==")
        for vt in ["SNV", "INDEL"]:
            for ca in CALLERS:
                s = e[(e.variant_type == vt) & (e.caller == ca)]["mean"]
                if len(s):
                    w(f"  {vt:5s} {ca:10s} mean {s.mean():+.4f}  max {s.max():+.4f}")
    if not rt.empty:
        w(f"== clean (exclusive) timings: {len(rt)} jobs ==")
    return "\n".join(L) + "\n"


def main():
    os.makedirs(OUT, exist_ok=True)
    res = load()
    # Primary analysis = RAW call sets, single-run scoring. Hard-filtered sets
    # are one operating point on the raw ROC; they get their own comparison.
    single = res[(res.scoring_method == "single_run") & (res.callset == "raw")].copy()

    summ = summarise(single)
    eff = effects(single)
    t = titv(single)
    s = scoring_method(res)
    rt = runtime(res)
    hf = hard_filter_effect(res)
    if not hf.empty:
        hf.to_csv(os.path.join(OUT, "hard_filter_effect.tsv"), sep="\t", index=False,
                  float_format="%.6g")

    summ.to_csv(os.path.join(OUT, "summary_by_condition.tsv"), sep="\t", index=False,
                float_format="%.6g")
    eff.to_csv(os.path.join(OUT, "effects.tsv"), sep="\t", index=False, float_format="%.6g")
    if not t.empty:
        t.to_csv(os.path.join(OUT, "titv_experiment.tsv"), sep="\t", index=False,
                 float_format="%.6g")
    s.to_csv(os.path.join(OUT, "scoring_method_effect.tsv"), sep="\t", index=False,
             float_format="%.6g")
    if not rt.empty:
        rt.to_csv(os.path.join(OUT, "runtime_by_condition.tsv"), sep="\t", index=False)

    txt = headline(eff, t, s, rt)
    with open(os.path.join(OUT, "headline.txt"), "w") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
