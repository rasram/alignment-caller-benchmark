#!/usr/bin/env python3
"""
STEP 8 — Predict which pipeline wins under which conditions.

DATA. E. coli only. phiX saturates — every pipeline scores F1 = 1.0 at almost
every condition — so it carries no information about which pipeline is better
and would only teach the model that "phiX => 1.0".

FEATURES. coverage, read_length, mean_p (MEASURED per-base error rate, not
ART's arbitrary -qs knob — NOTES 4.4), variant type, aligner, caller.
TARGET.   F1 of that pipeline on that run.

MODELS
  pipeline_mean   condition-blind: each pipeline's training-set mean F1
  tree            DecisionTreeRegressor, depth 4 — the PRIMARY model, because
                  its rules can be read, defended, and put on a slide
  forest          RandomForestRegressor — fitted only for permutation importance

VALIDATION — two, because they answer different questions
  held-out seeds    train seeds 1-3, test seeds 4-5. Splitting by SEED, never
                    by row: rows sharing a seed share the same reads, so a random
                    row split trains and tests on near-duplicates and reports a
                    flattering, meaningless score.
  leave-one-        train on 10 conditions, predict the 11th. This is the real
  condition-out     use case (a condition you have not benchmarked) and it is
                    harsh: holding out 5x or 100x coverage forces EXTRAPOLATION,
                    which trees cannot do. Reported honestly, not hidden.

DECISION METRIC — selection regret
  For each test (condition, variant type, seed): the F1 lost by picking the
  model's top-ranked pipeline instead of the true best one. Compared against
  the trivial rule "always pick the training-set's best-on-average pipeline".
  If the model cannot beat that rule, the honest conclusion is that one
  pipeline dominates — itself a finding.

Outputs -> results/model/
Usage: python3 scripts/fit_model.py
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.tree import DecisionTreeRegressor, export_text

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", "model")
SEED = 20260814
NUM = ["coverage", "read_length", "mean_p", "is_indel"]
ALIGNERS = ["bwa", "bowtie2", "minimap2"]
CALLERS = ["gatk", "freebayes", "bcftools"]
GROUP = ["coverage", "read_length", "qs_shift", "variant_type", "seed"]


def load():
    df = pd.read_csv(os.path.join(REPO, "results", "results.tsv"), sep="\t")
    df = df[(df.scoring_method == "single_run") & (df.callset == "raw")
            & (df.genome == "ecoli")].copy()
    if df.empty:
        sys.exit("no E. coli single_run rows in results.tsv")
    for c in ["coverage", "read_length", "qs_shift", "seed", "f1", "mean_p"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["is_indel"] = (df.variant_type == "INDEL").astype(int)
    df["pipeline"] = df.aligner + "+" + df.caller
    for a in ALIGNERS:
        df[f"aln_{a}"] = (df.aligner == a).astype(int)
    for c in CALLERS:
        df[f"cal_{c}"] = (df.caller == c).astype(int)
    return df


FEATS = NUM + [f"aln_{a}" for a in ALIGNERS] + [f"cal_{c}" for c in CALLERS]


TIE = 1e-12


def top_set(g, pred_col):
    """Rows the model ranks first. A regression tree predicts one value per leaf,
    so several pipelines often TIE for the top prediction — the model is then
    indifferent between them. idxmax() would silently break the tie by row order
    (alphabetical pipeline name), crediting or blaming the model for a choice it
    never made. Every tied row is returned instead."""
    return g[g[pred_col] >= g[pred_col].max() - TIE]


def regret(test, pred_col):
    """Mean F1 lost by choosing the model's top pipeline instead of the true best.
    Ties are scored as a uniform random choice among the tied pipelines, i.e.
    best minus the MEAN F1 of the tied set (expected regret)."""
    r = []
    for _, g in test.groupby(GROUP):
        r.append(g.f1.max() - top_set(g, pred_col).f1.mean())
    return float(np.mean(r)), len(r)


def fit_eval(train, test, label):
    rows = []
    # trivial selection rule + condition-blind predictor
    pm = train.groupby("pipeline").f1.mean()
    always = pm.idxmax()
    test = test.copy()
    test["pred_pipeline_mean"] = test.pipeline.map(pm)
    test["pred_always"] = (test.pipeline == always).astype(float)

    tree = DecisionTreeRegressor(max_depth=4, min_samples_leaf=5, random_state=SEED)
    tree.fit(train[FEATS], train.f1)
    test["pred_tree"] = tree.predict(test[FEATS])

    forest = RandomForestRegressor(n_estimators=400, min_samples_leaf=3,
                                   random_state=SEED, n_jobs=-1)
    forest.fit(train[FEATS], train.f1)
    test["pred_forest"] = forest.predict(test[FEATS])

    for name, col in [("pipeline_mean", "pred_pipeline_mean"), ("tree", "pred_tree"),
                      ("forest", "pred_forest")]:
        rg, n = regret(test, col)
        rows.append(dict(validation=label, model=name,
                         mae=mean_absolute_error(test.f1, test[col]),
                         r2=r2_score(test.f1, test[col]), regret=rg, decisions=n))
    rg, n = regret(test, "pred_always")
    # One stable model name across folds; the pipeline it chose is kept separately.
    # (Naming the row after the winner split leave-one-condition-out into several
    # rows, because different folds have different training-set winners.)
    rows.append(dict(validation=label, model="always_train_best", mae=np.nan, r2=np.nan,
                     regret=rg, decisions=n, picked=always))
    # context: expected regret of picking a pipeline uniformly at random
    rnd = float(np.mean([g.f1.max() - g.f1.mean() for _, g in test.groupby(GROUP)]))
    rows.append(dict(validation=label, model="random_pick", mae=np.nan, r2=np.nan,
                     regret=rnd, decisions=n))
    return rows, tree, forest, test


def main():
    os.makedirs(OUT, exist_ok=True)
    df = load()
    seeds = sorted(df.seed.unique())
    conds = df[["coverage", "read_length", "qs_shift"]].drop_duplicates()
    print(f"E. coli rows {len(df)} | seeds {seeds} | conditions {len(conds)}")

    results = []

    # --- A. held-out seeds ---------------------------------------------------
    if len(seeds) >= 3:
        tr_seeds = seeds[:max(1, len(seeds) - 2)]
        te_seeds = [s for s in seeds if s not in tr_seeds]
        rows, _, _, _ = fit_eval(df[df.seed.isin(tr_seeds)], df[df.seed.isin(te_seeds)],
                                 f"held-out seeds {[int(s) for s in te_seeds]}")
        results += rows

    # --- B. leave-one-condition-out ------------------------------------------
    loco = []
    per_cond = []
    for _, c in conds.iterrows():
        m = ((df.coverage == c.coverage) & (df.read_length == c.read_length)
             & (df.qs_shift == c.qs_shift))
        if m.sum() == 0 or (~m).sum() == 0:
            continue
        rows, _, _, te = fit_eval(df[~m], df[m], "leave-one-condition-out")
        loco += rows
        lab = f"cov{int(c.coverage)} len{int(c.read_length)} qs{int(c.qs_shift)}"
        for r in rows:
            per_cond.append(dict(r, held_out=lab))
    if loco:
        L = pd.DataFrame(loco)
        # Regret is averaged per DECISION (weighted by fold size), not per fold.
        L["regret_x_n"] = L.regret * L.decisions
        agg = (L.groupby(["validation", "model"])
               .agg(mae=("mae", "mean"), r2=("r2", "mean"), regret_x_n=("regret_x_n", "sum"),
                    decisions=("decisions", "sum")).reset_index())
        agg["regret"] = agg.regret_x_n / agg.decisions
        agg = agg.drop(columns="regret_x_n")
        results += agg.to_dict("records")
        pd.DataFrame(per_cond).to_csv(os.path.join(OUT, "loco_by_condition.tsv"),
                                      sep="\t", index=False, float_format="%.6g")

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(OUT, "validation.tsv"), sep="\t", index=False, float_format="%.6g")

    # --- final models on ALL data: readable rules + importances --------------
    tree = DecisionTreeRegressor(max_depth=4, min_samples_leaf=5, random_state=SEED)
    tree.fit(df[FEATS], df.f1)
    with open(os.path.join(OUT, "tree_rules.txt"), "w") as fh:
        fh.write("# Decision tree (depth 4) predicting F1 — fitted on all E. coli runs\n")
        fh.write("# aln_*/cal_* are 0/1 indicators; mean_p is the measured error rate\n\n")
        fh.write(export_text(tree, feature_names=FEATS, decimals=4))
    import pickle
    with open(os.path.join(OUT, "tree.pkl"), "wb") as fh:
        pickle.dump({"tree": tree, "features": FEATS}, fh)

    forest = RandomForestRegressor(n_estimators=400, min_samples_leaf=3,
                                   random_state=SEED, n_jobs=-1).fit(df[FEATS], df.f1)
    pi = permutation_importance(forest, df[FEATS], df.f1, n_repeats=20,
                                random_state=SEED, n_jobs=-1)
    imp = (pd.DataFrame({"feature": FEATS, "importance_mean": pi.importances_mean,
                         "importance_sd": pi.importances_std})
           .sort_values("importance_mean", ascending=False))
    # group one-hot indicators back into their factor
    imp["factor"] = imp.feature.str.replace(r"^aln_.*", "aligner", regex=True) \
                               .str.replace(r"^cal_.*", "caller", regex=True)
    imp.to_csv(os.path.join(OUT, "importances.tsv"), sep="\t", index=False, float_format="%.6g")

    # --- recommendation table: per condition x variant type -----------------
    df["pred"] = tree.predict(df[FEATS])
    rec = []
    for (cov, ln, qs, vt), g in df.groupby(["coverage", "read_length", "qs_shift", "variant_type"]):
        mean_by = g.groupby("pipeline")[["f1", "pred"]].mean()
        tied = top_set(mean_by, "pred")          # the tree's top-ranked set (may be >1)
        rec.append(dict(coverage=cov, read_length=ln, qs_shift=qs, variant_type=vt,
                        actual_best=mean_by.f1.idxmax(), actual_best_f1=mean_by.f1.max(),
                        model_top=",".join(sorted(tied.index)), model_top_n=len(tied),
                        model_top_mean_f1=tied.f1.mean(),
                        best_in_model_top=mean_by.f1.idxmax() in tied.index,
                        worst=mean_by.f1.idxmin(), worst_f1=mean_by.f1.min()))
    pd.DataFrame(rec).to_csv(os.path.join(OUT, "recommendations.tsv"), sep="\t",
                             index=False, float_format="%.6g")

    with pd.option_context("display.width", 140, "display.max_columns", 20):
        print(res.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
        print("\nPermutation importance (by factor):")
        print(imp.groupby("factor").importance_mean.sum().sort_values(ascending=False)
              .to_string(float_format=lambda x: f"{x:.6f}"))


if __name__ == "__main__":
    main()
