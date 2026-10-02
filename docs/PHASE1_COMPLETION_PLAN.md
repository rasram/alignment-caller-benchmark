# Phase 1 completion plan — what is done, what is left, and exactly how to finish

Status report and step-by-step instructions to take the project from its current state
to a complete, defensible semester deliverable.

**Verified state at time of writing:** `bash scripts/verify_all.sh` → **83 PASS / 0 FAIL / 0 WARN**.
Only the baseline condition has been executed: `results/results.tsv` contains **72 rows**
covering **2 of 110** condition×seed combinations.

---

## PART A — What is implemented

All of this is built, executed and verified. Evidence is reproducible via
`bash scripts/verify_all.sh`.

| # | Component | Evidence |
|---|---|---|
| 1 | **Environment** — 5 conda envs, pinned + locked; RTG Tools 3.13; simuG | `envs/*.yaml`, `logs/versions_latest.txt` |
| 2 | **References** — phiX (5,386 bp) + *E. coli* (4,641,652 bp), renamed, 5 index types each | `data/refs/`, lengths asserted |
| 3 | **Truth sets** — 50+10 and 5,000+1,000 variants, normalised, haploid `GT=1` | `data/truth/*.truth.vcf.gz` |
| 4 | **Read simulation** — ART HS25, 150 bp PE, exactly 30.00× both genomes, truth SAM kept | `work/*_1.fq`, `work/*_.sam` |
| 5 | **QC** — FastQC + MultiQC; mean-Q extraction with error-rate conversion | `results/qc/` |
| 6 | **Alignment** — 3 aligners × 2 genomes, indexed, duplicate-marked, read groups verified | `work/*.md.bam` |
| 7 | **Placement accuracy** — with mutated→reference coordinate conversion, self-validated | `results/placement_accuracy.tsv` |
| 8 | **Variant calling** — 9 pipelines × 2 genomes = 18 call sets, raw + hard-filtered | `work/*.raw.vcf.gz` |
| 9 | **Ploidy verification** — all 18 confirmed haploid, re-checked live by the audit | `logs/ploidy_verification.txt` |
| 10 | **GA4GH scoring** — `rtg vcfeval`, SNV and indel separately, two scoring methods | `results/vcfeval/` |
| 11 | **Results table + ROC curves** | `results/results.tsv`, `results/roc_*.svg` |
| 12 | **Snakemake workflow** — reproduces the baseline byte-for-byte; 990-run DAG resolves | `Snakefile`, 13,320 jobs |
| 13 | **Verification harness** — 83 checks, re-derived from data not file presence | `scripts/verify_all.sh` |
| 14 | **Documentation** — NOTES (1,551 lines), README, HANDOFF, setup.sh | repo root |
| 15 | **Mid-sem deck** — 14 slides with speaker notes | `docs/midsem_presentation.pptx` |

### Three findings that changed the method

These are the project's real intellectual content and belong in the final report:

1. **Ploidy cannot be verified from the score.** A caller emitting `1/1` scores a *perfect* F1
   against haploid truth. Only grepping the `GT` field detects the misconfiguration.
2. **Coordinate systems differ silently.** The simulator's truth SAM uses mutated-genome
   coordinates; aligners report reference coordinates. Same contig name, no error.
   Uncorrected: placement accuracy reads **8.9%** instead of **99.0%**.
3. **Pre-splitting VCFs by type before `vcfeval` distorts results.** It degrades a
   haplotype-aware comparison into a context-free one, penalising FreeBayes hardest.

---

## PART B — What is pending

| # | Item | Status | Blocking? |
|---|---|---|---|
| B1 | Metric collection is **not** workflow-integrated | 6 scripts are manual | **YES — blocks the sweep** |
| B2 | Disk management for the sweep | no `temp()` on reads/truth SAM | **YES — 50 GB vs 15 GB** |
| B3 | Ti/Tv realism decision | unresolved | **YES — must precede the sweep** |
| B4 | Parameter sweep execution | 2 of 110 conditions run | core deliverable |
| B5 | Seed-variance estimate | 1 seed, no error bars | core deliverable |
| B6 | Predictive model | not started | core deliverable |
| B7 | Final report / write-up | not started | core deliverable |
| B8 | hap.py cross-validation | not started | optional |
| B9 | Real read data (ERX008638) | not started | optional |

### B1 in detail — the one genuinely blocking defect

Six scripts produce the metric columns, and **none is a Snakemake rule**:

```
align_metrics.sh   placement_accuracy.py   call timing
extract_mean_q.py  run_qc.sh               make_roc.sh
```

`collect_results.py` joins them by `(tag, aligner)`. Those tables currently contain **only the
two baseline tags**, so for all 108 unrun condition×seed combinations the join finds nothing
and the metric columns come out **empty** — not stale, *empty*:

```
mapping_rate  mean_mapq  mean_depth  placement_accuracy  align_seconds  call_seconds  peak_rss_mb
```

Those are precisely the features the Phase-B6 model needs. A sweep run today produces 3,960
rows with the scientific columns correct and the feature columns blank.

**The fix is not simply "wrap the scripts in rules."** Each script loops over all aligners and
**appends** to one shared TSV. Under Snakemake's parallel execution, concurrent appends to one
file interleave and corrupt rows. The work must be restructured into one output file per unit,
then aggregated — the standard scatter-gather pattern.

---

## PART C — Critical path

```
C0 housekeeping
      │
C1 metrics → Snakemake rules  ──┐
C2 disk / temp()              ──┤
C3 Ti/Tv decision             ──┘
      │
C4 pilot sweep (1 condition × 5 seeds)   ← GATE: metric columns must be populated
      │
C5 full sweep (990 runs)
      │
C6 validate sweep output                 ← GATE: no empty columns, no F1 = 0
      │
C7 variance analysis  →  C8 model  →  C9 figures  →  C10 report
```

C1–C3 are independent of each other and can be done in any order, but **all three must precede
C4**. Nothing after C4 can start until the sweep data exists.

---

## PART D — Step-by-step instructions

Time estimates assume working sessions, not elapsed days.

---

### STEP 0 — Housekeeping (15 min)

PowerPoint writes a lock file (`~$midsem_presentation.pptx`) whenever the deck is open, and one
was committed by accident.

```bash
cd "/Users/rash/Documents/Files/Self/College/Semester 7/DNA/project/itr1"
printf '\n# Office lock files, created while a document is open\n~$*\n' >> .gitignore
git rm --cached -f 'docs/~$midsem_presentation.pptx' 2>/dev/null || true
git add -A && git commit -m "Ignore Office lock files"
```

**Verify:** `git status --short` is clean.

---

### STEP 1 — Make metric collection workflow-integrated (3–4 h) · **BLOCKING**

The goal: after the sweep, every row of `results.tsv` has its metric columns filled.

#### 1a. Use Snakemake's built-in benchmarking for timing

Delete the hand-rolled `/usr/bin/time -l` plumbing. Add to **every** `align_*`, `call_*` and
`mark_duplicates` rule:

```python
benchmark: "benchmarks/{tag}.{aligner}.bwa.tsv"      # adjust per rule
```

Snakemake writes `s`, `h:m:s`, `max_rss`, `max_vms`, `mean_load` and more. This is portable —
it sidesteps the BSD-vs-GNU `time` difference (`-l` vs `-v`, bytes vs kilobytes) that caused
trouble earlier.

#### 1b. Convert `align_metrics.sh` to one file per (tag, aligner)

```python
rule align_metrics:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
    output: "work/metrics/{tag}.{aligner}.align.tsv"   # ONE unit, no appending
    conda: "envs/align.yaml"
    shell: "bash scripts/align_metrics_one.sh {wildcards.tag} {wildcards.aligner} > {output}"
```

Write `scripts/align_metrics_one.sh` as a single-unit version of the existing script: same
`flagstat` parsing, same **primary**-count fairness fix (R8), but emitting exactly one row to
stdout and taking no loop over aligners.

#### 1c. Same treatment for placement accuracy

`placement_accuracy.py` already accepts `--bam` / `--aligner`; add a mode that writes a single
row to a given path instead of appending:

```python
rule placement_accuracy:
    input:
        truth_sam = "work/{tag}_.sam",
        bam       = "work/{tag}.{aligner}.md.bam",
        indels    = lambda w: f"data/truth/{genome_of(w.tag)}.refseq2simseq.INDEL.vcf",
    output: "work/metrics/{tag}.{aligner}.placement.tsv"
    conda: "envs/ml.yaml"
    shell:
        "python3 scripts/placement_accuracy.py --truth-sam {input.truth_sam} "
        "--indel-vcf {input.indels} --bam {input.bam} --aligner {wildcards.aligner} "
        "--tag {wildcards.tag} --out {output} --single"
```

#### 1d. Add mean-Q as a rule

```python
rule mean_q:
    input: "work/{tag}_1.fq", "work/{tag}_2.fq"
    output: "work/metrics/{tag}.meanq.tsv"
    conda: "envs/ml.yaml"
    shell: "python3 scripts/extract_mean_q.py --out {output} {input}"
```

`mean_p` is the physically meaningful error-rate feature — the model should use it rather than
`qs_shift`, which is an arbitrary ART knob.

#### 1e. Aggregate, then collect

```python
rule aggregate_metrics:
    input:
        align = expand("work/metrics/{tag}.{aligner}.align.tsv", tag=TAGS, aligner=ALIGNERS),
        place = expand("work/metrics/{tag}.{aligner}.placement.tsv", tag=TAGS, aligner=ALIGNERS),
        meanq = expand("work/metrics/{tag}.meanq.tsv", tag=TAGS),
    output:
        align = "results/align_metrics.tsv",
        place = "results/placement_accuracy.tsv",
        meanq = "results/qc/mean_q.tsv",
    shell:
        r"""
        head -1 {input.align[0]} > {output.align}; tail -q -n +2 {input.align} >> {output.align}
        head -1 {input.place[0]} > {output.place}; tail -q -n +2 {input.place} >> {output.place}
        head -1 {input.meanq[0]} > {output.meanq}; tail -q -n +2 {input.meanq} >> {output.meanq}
        """
```

Then make `collect_results` depend on `aggregate_metrics`, and extend it to read Snakemake's
`benchmarks/` files for `align_seconds` / `call_seconds` / `peak_rss_mb`.

**Verify step 1:**

```bash
rm -f results/align_metrics.tsv results/placement_accuracy.tsv
snakemake --use-conda --cores 8 --forceall
"$HOME/miniforge3/envs/ml/bin/python" - <<'PY'
import csv
rows=list(csv.DictReader(open('results/results.tsv'),delimiter='\t'))
cols=["mapping_rate","mean_mapq","mean_depth","placement_accuracy",
      "align_seconds","call_seconds","peak_rss_mb"]
blank=[c for c in cols if any(r[c]=="" for r in rows)]
print("EMPTY COLUMNS:", blank or "none — PASS")
PY
```

---

### STEP 2 — Disk management (45 min) · **BLOCKING**

Measured footprint per 30× *E. coli* condition, retained:

| Artefact | Size |
|---|---|
| 2 × FASTQ | 284 MB |
| ART truth SAM | 315 MB |
| 3 × duplicate-marked BAM | 249 MB |
| **Total** | **≈ 0.83 GB** |

Across the sweep (coverage scale factor 12.17 × 5 seeds):

| Strategy | *E. coli* total |
|---|---|
| Keep everything | **≈ 50 GB** |
| Reads + truth SAM as `temp()` | **≈ 15 GB** |

You currently have **216 GB free**, so 50 GB fits. But 15 GB is better and costs one line each.

Mark the FASTQ outputs `temp()` — Snakemake deletes them once the last aligner has consumed
them. The truth SAM is the subtle one: **it cannot simply be `temp()`'d**, because
`placement_accuracy` needs it. Mark it `temp()` *and* ensure the placement rule lists it as an
input; Snakemake then deletes it after the last placement job for that tag.

```python
rule simulate_reads:
    output:
        r1  = temp("work/{tag}_1.fq"),
        r2  = temp("work/{tag}_2.fq"),
        sam = temp("work/{tag}_.sam"),
```

**Caution:** with `temp()`, re-running a downstream rule later forces re-simulation. That is
correct behaviour and cheap (ART is fast), but do not be surprised by it.

**Verify:** `snakemake -n --config run=all` still resolves; then after step 4 confirm
`du -sh work` stays bounded.

---

### STEP 3 — Resolve the Ti/Tv question (15 min, or 30 min if regenerating) · **BLOCKING**

The truth sets use simuG's default transition/transversion ratio ≈ **0.5** — what uniformly
random substitution produces. Real bacterial genomes run **1–2**.

- It does **not** bias the comparison: all nine pipelines score against the same truth and none
  is Ti/Tv-tuned, so the **ranking holds**.
- It **does** mean absolute recall is not directly transferable to a real resequencing project.

**Decide now, because changing it afterwards invalidates every result.**

If you choose to regenerate:

```bash
rm -f data/truth/{phiX,ecoli}.refseq2simseq.* data/truth/{phiX,ecoli}.simseq.genome.fa
for g in phiX ecoli; do
  snp=$([ $g = phiX ] && echo 50 || echo 5000)
  ind=$([ $g = phiX ] && echo 10 || echo 1000)
  "$HOME/miniforge3/envs/sim/bin/perl" tools/simuG/simuG.pl \
    -refseq data/refs/${g}.fa -snp_count $snp -indel_count $ind \
    -titv_ratio 2.0 -seed 20260814 -prefix data/truth/${g}
  bash scripts/build_truth.sh $g && bash scripts/verify_truth.sh $g
done
snakemake --use-conda --cores 8 --forceall     # rebuild the baseline
```

**Recommendation:** regenerate with `-titv_ratio 2.0`. It costs ~30 minutes total and removes
the single biggest "but is this realistic?" objection from your final report.

---

### STEP 4 — Pilot the sweep (30 min) · **GATE**

Do **not** launch 990 runs against untested plumbing. Run one off-baseline condition at all
5 seeds first — this exercises a coverage level never used before *and* the new metric rules.

```bash
# temporarily restrict conditions.tsv to one non-baseline condition, e.g. cov=5
cp config/conditions.tsv config/conditions.full.tsv
head -1 config/conditions.full.tsv > config/conditions.tsv
grep -P '\tecoli\t5\t150\t0\t' config/conditions.full.tsv >> config/conditions.tsv
snakemake --use-conda --cores 8 --config run=all
```

**Gate criteria — all must hold before proceeding:**

1. All jobs complete, exit 0
2. `results.tsv` gains 5 × 18 = **90** new `single_run` rows
3. **No empty metric columns** (the check from step 1)
4. Ploidy verification passes for every new call set
5. Coverage actually measures ~5× — confirm ART honoured the parameter
6. `du -sh work` bounded as expected

Then restore: `mv config/conditions.full.tsv config/conditions.tsv`

---

### STEP 5 — Run the full sweep (~3 h wall-clock)

```bash
cp results/results.tsv results/results.baseline-backup.tsv
nohup snakemake --use-conda --cores 8 --config run=all \
      --rerun-incomplete --keep-going > logs/sweep.log 2>&1 &
```

- `--keep-going` finishes independent branches if one job fails, rather than halting everything
- `--rerun-incomplete` recovers cleanly if the machine sleeps mid-run

**Measured budget:** 120 s per full 9-pipeline *E. coli* condition at 8 cores; coverage scale
factor 12.17; 5 seeds → **≈ 2.25 h** plus overhead. Budget 3 h and **disable sleep** for the
duration.

Monitor: `tail -f logs/sweep.log` or `grep -c "steps.*done" logs/sweep.log`

---

### STEP 6 — Validate the sweep (45 min) · **GATE**

```bash
bash scripts/verify_all.sh
"$HOME/miniforge3/envs/ml/bin/python" - <<'PY'
import csv, collections
rows=[r for r in csv.DictReader(open('results/results.tsv'),delimiter='\t')
      if r['scoring_method']=='single_run']
print(f"single_run rows: {len(rows)}  (expected 110 x 18 = 1980)")
conds={(r['genome'],r['coverage'],r['read_length'],r['qs_shift'],r['seed']) for r in rows}
print(f"conditions covered: {len(conds)}  (expected 110)")
cols=["mapping_rate","mean_mapq","mean_depth","placement_accuracy",
      "align_seconds","call_seconds","peak_rss_mb"]
print("empty columns:", [c for c in cols if any(r[c]=="" for r in rows)] or "none")
zero=[r for r in rows if r['genome']=='ecoli' and float(r['f1'])==0]
print(f"E. coli rows with F1=0 (bug signature): {len(zero)}")
PY
```

**Expected:** 1,980 `single_run` rows (3,960 total including `pre_split`), 110 conditions,
no empty columns, no *E. coli* F1 = 0.

Also confirm `logs/ploidy_verification.txt` contains **990 PASS** entries and zero FAIL.

---

### STEP 7 — Variance analysis (2–3 h) — *the scientific core*

This answers the question the mid-sem talk could not: **are the pipeline differences real?**

At baseline the nine *E. coli* SNV F1 values span only **0.0047**. Nothing yet shows that
exceeds seed-to-seed noise.

Write `scripts/analyse_sweep.py` to produce:

1. **Per-condition mean ± sd across the 5 seeds**, for every pipeline and variant type.
2. **Effect size vs noise:** is the aligner spread larger than the within-pipeline sd?
3. **A repeated-measures comparison** — seeds are paired (all pipelines see the same seed), so
   compare pipelines *within* seed rather than treating runs as independent. A Friedman test
   (non-parametric, repeated measures) suits 9 related groups; follow with pairwise Wilcoxon
   and a multiple-comparison correction.
4. **The headline plot:** F1 vs coverage, one line per aligner, error bars = sd over seeds.

**Report the result honestly either way.** "The aligner effect exceeds seed noise above 10×
coverage but not below" is a finding. So is "the differences are within noise at 30×" — it
would mean the mid-sem headline needs qualifying, which is exactly what the sweep is for.

---

### STEP 8 — Predictive model (4–6 h)

Goal from the original brief: *predict which pipeline wins under which conditions.*

**Features** (per condition, not per pipeline): `coverage`, `read_length`, **`mean_p`**
(measured error rate, not the arbitrary `qs_shift`), `genome_size`, `variant_type`.

**Target:** two framings, both worth reporting —
- *Regression:* predict F1 for a given (condition, pipeline).
- *Classification:* which of the 9 pipelines has the highest F1 for this condition?

**Model choice:** a shallow **decision tree** (depth 3–4) is the right primary model. Not
because it is the most accurate, but because it is **readable** — it produces rules such as
"below 10× coverage, prefer BWA + GATK" which you can defend in a viva and put on a slide. Fit
a random forest alongside only to report feature importances.

**Validation:** split by **seed**, not randomly. Random splits leak — rows from the same seed
and condition are correlated, so a random split trains and tests on near-identical rows and
reports a flattering, meaningless score. Hold out seeds 4–5 entirely.

**Sanity floor:** compare against the trivial baseline "always pick BWA + GATK". If the model
cannot beat that, the honest conclusion is that one pipeline dominates everywhere — which is
itself a publishable result.

---

### STEP 9 — Figures (2–3 h)

| Figure | Content |
|---|---|
| F1 | F1 vs coverage, per aligner, error bars over seeds — **the headline** |
| F2 | F1 vs read length, same structure |
| F3 | F1 vs measured error rate (`mean_p`) |
| F4 | 3×3 heatmaps at low (5×) and high (100×) coverage, side by side |
| F5 | Decision tree diagram |
| F6 | Runtime vs accuracy scatter — the practical trade-off |
| F7 | ROC curves at a low-coverage condition, where they actually separate |

Reuse the deck's palette for consistency. Note `make_roc.sh` must `rm -f` its target first —
`rtg rocplot` refuses to overwrite.

---

### STEP 10 — Final report (8–12 h)

Suggested structure, with the material already written in `NOTES.md` mapped to each section:

| Section | Source |
|---|---|
| Introduction / motivation | deck slides 2–3 |
| Methods — design, tools, rules | `NOTES.md` phases 0–2, deck slide 4 |
| Methods — correctness controls | `NOTES.md` §2.3, §3.4, §6.6, §7.6 |
| Results — alignment layer | `results/align_metrics.tsv`, slide 7 |
| Results — baseline 3×3 | slide 8 |
| Results — sweep | STEP 7 output |
| Results — model | STEP 8 output |
| Discussion — limitations | `HANDOFF.md` §2 |
| Reproducibility | `README.md`, `verify_all.sh` |

**The three silent failures (deck slide 11) deserve their own methods subsection.** They are
the strongest evidence of understanding in the whole project — most student benchmarks report
numbers without ever demonstrating the numbers could have been wrong.

---

## PART E — Optional extensions

Only if time allows. **None is required for a complete Phase 1.**

| Item | Cost | Value |
|---|---|---|
| **hap.py cross-validation** | ~4 h | Confirms `vcfeval` scores with an independent GA4GH tool. Highest-value optional item — a second tool agreeing is strong evidence. |
| **Real read data (ERX008638)** | ~6 h | Tests whether simulated rankings hold on real *E. coli* reads. No truth set, so compare pipelines to each other, not to truth. |
| **More aligners/callers** | ~8 h | Breadth, not depth. Lowest value — the brief fixed the 3×3 deliberately. |

---

## PART F — Risks and decision points

| Risk | Mitigation |
|---|---|
| Sweep produces blank feature columns | STEP 4 pilot gates exactly this |
| Machine sleeps mid-sweep | `caffeinate -i`; `--rerun-incomplete` recovers |
| Differences turn out to be within noise | Report it — a null result is a result, and STEP 7 is designed to detect it |
| Disk fills | STEP 2 `temp()` cuts 50 GB → 15 GB |
| Ti/Tv objection in the viva | STEP 3 removes it for ~30 min of work |
| Model does not beat "always BWA+GATK" | Report the trivial baseline as the finding |

### Three questions to settle before STEP 5

1. **Regenerate the truth sets with realistic Ti/Tv?** (recommended: yes)
2. **Is the single-run scoring method acceptable** in place of the brief's pre-split method?
   Both are in `results.tsv`; `single_run` is primary and defensible, but if the marking scheme
   expects the brief's method literally, the primary column must change.
3. **Is 5 seeds enough?** If STEP 7 shows the effect is marginal, more seeds are cheaper than
   more conditions — each additional seed costs ~27 min.

---

## PART G — Suggested schedule

| Session | Steps | Hours | Output |
|---|---|---|---|
| 1 | 0, 1 | 4 | Metrics workflow-integrated |
| 2 | 2, 3, 4 | 3 | Pilot passes the gate |
| 3 | 5, 6 | 4 | Full sweep data validated |
| 4 | 7 | 3 | Variance analysis — are differences real? |
| 5 | 8 | 5 | Model fitted and validated |
| 6 | 9 | 3 | Figures |
| 7–8 | 10 | 10 | Final report |
| | **Total** | **≈ 32 h** | |

Plus ~4 h if you take the hap.py cross-validation.

---

## Definition of done for Phase 1

- [ ] `verify_all.sh` passes with the sweep data present
- [ ] `results/results.tsv` has 1,980 `single_run` rows, no empty columns
- [ ] 990 ploidy checks logged, zero failures
- [ ] Variance analysis states whether pipeline differences exceed seed noise
- [ ] Model fitted, validated by held-out **seed**, compared against the trivial baseline
- [ ] Figures F1–F6 generated
- [ ] Final report written, limitations section explicit
- [ ] `HANDOFF.md` updated to describe the finished state
