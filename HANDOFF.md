# HANDOFF

State of the project at the **end of Phase 1** (October 2026): the full benchmark has been run,
analysed and written up. This file records *what is true right now*, what was assumed, what is
not verified, and what Phase 2 should do. **[NOTES.md](NOTES.md)** explains *why* each decision
was made; **[docs/FINAL_REPORT.pdf](docs/FINAL_REPORT.pdf)** is the result.

---

## 1. What is done

`bash scripts/verify_all.sh` — **all checks pass** at full-sweep scope (110 run tags). Output:
[logs/verification_report.txt](logs/verification_report.txt).

| Area | Status |
|---|---|
| Truth sets regenerated with Ti/Tv 2.0 (same positions as the Ti/Tv 0.5 originals) | ✅ |
| Full sweep: 11 conditions × 5 seeds × 2 genomes × 9 pipelines = **990 runs** | ✅ |
| Every call set haploid — 990 VCFs re-grepped live, zero diploid genotypes | ✅ |
| Per-run metrics as Snakemake rules (read, alignment, placement, runtime, ploidy) | ✅ |
| Clean (machine-exclusive) timings for seed 1 of every condition | ✅ |
| Raw **and** hard-filtered call sets scored; single-run and pre-split scoring | ✅ |
| `results/results.tsv` — 5,940 rows, no empty metric cells (`--strict`) | ✅ |
| Variance analysis (blocked ANOVA + Friedman, BH-corrected) | ✅ |
| Error-mechanism tests (`scripts/diagnose_errors.py`) | ✅ |
| Predictive model (tree + forest; held-out-seed and leave-one-condition-out) | ✅ |
| Figures F1–F7, report (MD/PDF/DOCX) — every table and quoted number generated | ✅ |
| One command reproduces everything: `snakemake --config run=all` | ✅ |

**Compute actually used:** about 2.5 h wall-clock for the sweep on an Apple M4 (10 cores, 24 GB),
14,865 Snakemake jobs. Retained on disk: `work/` 17 GB (BAMs, VCFs), `results/vcfeval/` 0.9 GB,
workflow conda envs 2.1 GB. All of it is gitignored and regenerable.

---

## 2. Headline results

Full detail, tables and figures are in the report. In one paragraph: all nine pipelines are
accurate on simulated haploid data (*E. coli* 30× SNV F1 0.989–0.995). Their differences are
real (seed sd ≈ 0.0005, about a tenth of the aligner spread) and **driven mainly by the
aligner**. The aligner accounts for a median 0.78 of between-pipeline variation for SNVs, the
caller 0.02. Bowtie2's end-to-end alignment creates false SNVs beside true indels that grow with
depth; GATK's reassembly absorbs them. About 50 variants sit in repeats that no pipeline recovers.
Depth is the only data property that changes the ranking: at 5× BCFtools wins for SNVs. Error
rate (7× range) barely matters. A model predicts the best pipeline only marginally better than
"always BWA-MEM + GATK", except at low depth.

---

## 3. Assumptions and decisions you should know about

### 3.1 Primary scoring is one vcfeval run, not the brief's pre-split method — **confirm this**

The brief specified splitting call sets by type before scoring. Measured, that understates F1 by
up to 0.0105, and most for FreeBayes, so it changes conclusions rather than shifting all numbers
(report §3.8, [NOTES.md §7.6](NOTES.md)). Both methods are in `results.tsv` under
`scoring_method`; **`single_run` is primary**. If the marker expects the brief's method, the
primary column must change.

### 3.2 Ti/Tv regenerated to 2.0 after the mid-semester review

The original truth sets used simuG's default Ti/Tv 0.5. They were regenerated with 2.0 (measured
2.04 for *E. coli*), keeping the identical 5,000 SNV positions and identical indels. At baseline,
F1 changed by ≤ 0.0002 for BWA-MEM/minimap2 pipelines and ≤ 0.0024 for Bowtie2. The Ti/Tv 0.5
baseline is archived in `results/archive/titv0.5_baseline/`.

### 3.3 Other choices

| Choice | Why | Where |
|---|---|---|
| simuG seed 20260814 fixed; ART seeds 1–5 | the variants are the subject; seeds are replicate *sequencing runs* | NOTES Phase 3 |
| Hard filter `QUAL ≥ 20 && DP ≥ 5` | only fields all three callers emit comparably; QUAL is not equally calibrated | NOTES Phase 6 |
| 4 threads for aligners and GATK | identical thread budget (R8); FreeBayes/BCFtools are single-threaded | config.yaml |
| Clean timing for seed 1 only | runtime depends on condition, not seed; exclusive timing of all seeds would idle 7 cores | NOTES 8.3 |
| Default tool parameters except ploidy | an out-of-the-box comparison; tuning would favour whichever tool was tuned | report §4.3 |
| No BQSR/VQSR for GATK | needs a known-variant database these organisms lack | report §2.4 |
| Tree depth 4, min leaf 5, fixed before fitting | interpretability; not re-tuned after seeing results | report §2.10 |
| Model regret treats tied predictions as ties | a tree ranks whole leaves equal; breaking ties by name was a bug | NOTES 8.12 |

---

## 4. Corrections to earlier statements

- **Mid-semester deck numbers** (`docs/midsem_presentation.pptx`) are from the Ti/Tv 0.5 truth
  sets, seed 1 only. They are archived for traceability; the final numbers are in the report.
- **"All nine pipelines score F1 = 1.0 on phiX"** was true for seed 1 at baseline only. Across five
  seeds, 16 of 18 baseline cells are perfect. Every phiX error outside 5× is explained in report
  §3.7: a variant at position 51, near the end of the linearised circular genome, and low-QUAL
  BCFtools calls.
- **The earlier open question "is the 0.0047 SNV spread more than noise?"** is answered: the seed
  sd is about 0.0005, roughly ten times smaller than the aligner spread.

---

## 5. What is NOT verified

1. **Real data.** Everything is simulated. ART has no PCR duplicates, GC bias or contamination;
   simuG places variants uniformly. Real-data F1 will be lower.
2. **The Bowtie2 mechanism is tested by its prediction, not by intervention.** The false SNVs
   cluster beside true indels exactly as predicted (99% within 150 bp vs 6.3% background), but
   Bowtie2 was never rerun with `--local` to show they disappear. That is the decisive experiment.
3. **GATK's MAPQ-20 read filter** is inferred, not tested, as the reason Bowtie2 + GATK misses
   twice as many variants. Rerunning with `--minimum-mapping-quality` lowered would test it.
4. **`vcfeval` is the only scorer.** A `hap.py` cross-check was out of scope.
5. **Portability.** The workflow ran only on macOS arm64. Linux is untested.
6. **Interactions between axes** (e.g. 5× with 75 bp reads) were never sampled; the design is
   one-factor-at-a-time.

---

## 6. What Phase 2 should do first

1. **Rerun Bowtie2 with `--local`** at 30× and 100× (a small, targeted run). It confirms or
   kills the report's main mechanistic claim.
2. **Real data**: an *E. coli* run with an independent closed assembly as truth, to test whether
   the simulated ranking holds.
3. **`hap.py` cross-check** of `vcfeval` on a subset.
4. **A factorial design at low depth**, where the ranking actually changes.

---

## 7. Known traps, so you do not rediscover them

| Trap | Symptom |
|---|---|
| Running Snakemake without `--use-conda` | every rule fails with **exit 127** (the default profile sets it) |
| Snakemake cannot find `conda` | `Error running conda info` — add `~/miniforge3/bin` to `PATH` |
| `gatk` called by absolute path | `env: python: No such file or directory` — it is a Python launcher |
| `grep` on simuG's mutated FASTA | silently truncates the single 4.6 Mb line; use `awk` |
| `cmd \| head` / `echo \| awk '{…; exit}'` under `pipefail` | intermittent exit 141 (SIGPIPE); use here-strings |
| ART truth SAM vs aligner BAM | different coordinate systems, same contig name — no error, wrong answers |
| RTG output directory already exists | refuses to overwrite; `rm -rf` the target first |
| macOS `bash` 3.2, BSD `/usr/bin/time` | no associative arrays; `-l` reports bytes, not KB |
| `snakemake benchmark:` on macOS | `max_rss` is NA — use `scripts/lib/measure.sh` |
| `kill -INT` on a backgrounded Snakemake | ignored (POSIX: `&` from a script ignores SIGINT); use SIGTERM |
| GATK `--native-pair-hmm-threads 4` | reserves 4 cores, uses ~1; run the sweep with `--cores 12` |
| pandas `df.cov`, `df.pipe` | DataFrame **methods**, not your columns — use `df["cov"]` |
| `rule all` lists the report before its template exists | whole DAG fails to validate; give explicit targets |
| pandoc pipe tables with `\|---\|---\|` | every column the same width; size dashes by content |
| XeLaTeX + Helvetica | no `→` glyph, silently dropped — check the log for "Missing character" |
