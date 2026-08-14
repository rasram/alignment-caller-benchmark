# HANDOFF

State of the project at the end of the Phase 0–7b session, for whoever picks up Phase 8.

Read this before running anything. **[NOTES.md](NOTES.md)** explains *why* each decision was
made; this file records *what is true right now*, what was assumed, and what is not yet
verified.

---

## 1. What works

Verified by `bash scripts/verify_all.sh` — **81 of 83 checks pass**, and the two failures were
this file and `README.md` not yet existing. Full output in
[logs/verification_report.txt](logs/verification_report.txt).

| Area | Status |
|---|---|
| Five conda environments, pinned + locked | ✅ |
| RTG Tools 3.13 (Java 23), simuG `0289e58` | ✅ |
| Both references downloaded, renamed, 5 index types each | ✅ |
| Truth sets: phiX 50 SNV + 10 indel; *E. coli* 5,000 + 1,000 | ✅ |
| Baseline reads, both genomes, exactly 30.00× | ✅ |
| FastQC/MultiQC + mean-Q extraction | ✅ |
| 3 aligners → indexed, duplicate-marked BAMs with read groups | ✅ |
| Placement accuracy (with mutated→reference coordinate conversion) | ✅ |
| 9 aligner×caller VCFs per genome (18 total), raw + hard-filtered | ✅ |
| **Ploidy haploid in all 18, logged and re-verified live** | ✅ |
| `rtg vcfeval` scoring, SNV and indel separately | ✅ |
| `results/results.tsv` — 36 baseline rows + 36 comparison rows | ✅ |
| ROC curves, both genomes, both variant types | ✅ |
| Snakemake workflow reproducing the baseline **byte-for-byte** | ✅ |
| Full 990-run sweep DAG resolves (13,320 jobs) | ✅ |

**The single most reassuring result:** the baseline was rebuilt from scratch through Snakemake
in freshly created conda environments, using command lines written independently of the shell
scripts, and `results.tsv` came out byte-for-byte identical. Two independent implementations
agree.

---

## 2. Assumptions you should verify

These are choices made without explicit instruction. None is hidden in the code; all are
listed here because they affect how the results should be read.

### 2.1 Ti/Tv ratio is ~0.5, not biologically realistic — **most important**

simuG's default `-titv_ratio 0.5` was kept. That is what uniformly random substitution
produces (each base has one transition partner and two transversion partners). **Real bacterial
genomes run Ti/Tv ≈ 1–2.** Measured: phiX 0.667, *E. coli* 0.480.

*Impact:* does **not** bias the comparison — all nine pipelines score against the same truth
and none is Ti/Tv-tuned, so the **ranking holds**. But absolute recall is not directly
transferable to a real resequencing project. Fixing it means regenerating the truth sets with
`-titv_ratio 2.0` and re-running everything (~2.5 min for the baseline).

### 2.2 ART seed 1 for the baseline; simuG seed 20260814 fixed project-wide

The mutated genome is the *experimental subject* and must not change between conditions, so
simuG's seed is fixed. ART's seed varies 1–5 across the sweep as replicate sequencing runs.
**The baseline reported here is a single seed** — no error bars. The sweep's 5 seeds will give
the first estimate of run-to-run variance, and that variance is needed before claiming any
pipeline difference is real (see §5, Q2).

### 2.3 Hard filter is `QUAL>=20 && INFO/DP>=5`

Chosen because QUAL and DP are the only fields all three callers emit comparably. GATK Best
Practices would use `QD`/`FS`/`MQRankSum`, but FreeBayes and BCFtools do not produce them, so
using them would tune GATK's filter and not its competitors'. **QUAL is not calibrated
identically across the three tools**, so one threshold is *procedurally* identical but not
*statistically* equivalent. This is why the ROC curves are the honest comparison and the
hard-filtered numbers are one operating point on them.

### 2.4 Primary scoring uses one vcfeval run, not the brief's pre-split method

The brief specified splitting the VCF by type before scoring. Measured, that **distorts
results** — it degrades vcfeval's haplotype-aware comparison into a context-free one. Both are
in `results.tsv` under `scoring_method`; **`single_run` is primary**. Detail and evidence in
[NOTES.md §7.6](NOTES.md). *This is a deviation from the brief and should be confirmed.*

### 2.5 Threads fixed at 4 for aligners; callers differ in threading

R8 requires an identical thread count across aligners, and that is enforced. But the callers
are not equally parallel: GATK takes `--native-pair-hmm-threads 4`, BCFtools is effectively
single-threaded for the pileup, and FreeBayes is single-threaded. **Caller runtimes are
therefore not a like-for-like speed comparison** and should be reported as "as typically run",
not as a controlled benchmark.

### 2.6 MarkDuplicates is inert here

~0% duplicates, as expected — ART simulates no PCR. It is run for pipeline realism and cannot
be a source of difference between pipelines. Say so in the report rather than presenting it as
a result.

---

## 3. What is NOT verified

Be careful with these.

1. **`align_metrics.tsv`, `placement_accuracy.tsv` and `call_timing.tsv` are not Snakemake
   rules.** They are produced by the Phase 5/6 shell scripts. When the baseline was rebuilt
   through Snakemake, those columns in `results.tsv` were **carried over, not recomputed**. The
   scientific columns (TP/FP/FN/precision/recall/F1) *were* fully regenerated. **Fix this
   before the sweep** — otherwise every condition inherits the baseline's metrics. This is
   Phase 8 item #1.

2. **`--use-conda` was exercised only on this machine (macOS arm64).** The environments solve
   and the workflow runs, but portability to Linux is untested.

3. **Only the baseline condition has been executed.** Nothing at 5×, 100×, 75 bp, or any qs
   shift has ever run. The DAG resolves for them; the *commands* are unexercised at those
   parameter values.

4. **phiX cannot discriminate between pipelines.** All nine score F1 = 1.0000. Verified real
   (not a bug) with a negative control, but phiX is a smoke test, not a result. Do not put a
   9-way phiX comparison in a report as if it means something.

5. **The `filt` call sets exist but were never scored.** Only `raw` was pushed through vcfeval.
   Scoring them is a one-line change (`score_variants.sh ... filt`).

6. **`.snakemake/` was briefly committed to git and then removed by amending the commit.** The
   history was rewritten before any push. If you cloned this repo *very* early, re-clone.

---

## 4. Compute budget for the sweep

**Measured**, not estimated: one full *E. coli* 9-pipeline baseline condition = **120 s
wall-clock on 8 cores** (67 Snakemake jobs). phiX adds ~13 s.

Scaling by coverage (compute tracks total bases sequenced; read length and qs shift do not
change data volume), the 11 conditions sum to **12.17×** a single 30× condition:

| | Estimate |
|---|---|
| *E. coli* sweep (11 conditions × 5 seeds) | **~2.0 h** |
| phiX sweep | ~0.2 h |
| **Total, 8 cores** | **~2.25 h** |

### This is far cheaper than the brief's ~50 h estimate — but disk is the real constraint

| | |
|---|---|
| Intermediates per 30× *E. coli* condition | 1.32 GB (2 FASTQ + ART truth SAM + 3 BAMs) |
| **All *E. coli* intermediates if retained** | **~80 GB** |
| Free space on this machine now | ~121 GB |

It fits, but not comfortably, and the 100× conditions alone account for ~22 GB. **Recommended
before running the sweep:** mark the FASTQ and truth-SAM outputs `temp()` in the Snakefile so
Snakemake deletes them once the BAMs exist. The BAMs are what the callers need; the reads are
regenerable from the recorded seed. That cuts peak disk roughly in half.

The 2.25 h figure assumes linear scaling in coverage and this machine's 8 cores. Treat it as
the right order of magnitude, not a promise.

---

## 5. What Phase 8 should do first

**In this order.**

1. **Make the metrics collection Snakemake rules** (§3.1). Without this the sweep produces
   36 × 11 × 5 rows whose mapping-rate, MAPQ, depth, placement-accuracy and timing columns are
   all copies of the baseline. This is the one blocking defect.

2. **Add `temp()` to the read outputs** (§4) so the sweep does not run the disk out.

3. **Decide the Ti/Tv question** (§2.1) *before* burning compute. Regenerating the truth sets
   afterwards invalidates every result.

4. **Then run the sweep**, ideally coverage-first — the coverage axis is where pipeline
   differences should appear, since all three aligners are at ~99% placement accuracy at 30×
   and the baseline barely separates them.

5. **Only then** start any modelling. With one seed per condition there is no variance
   estimate; with five there is a weak one.

---

## 6. The three questions I most want answered

1. **Is the deviation in §2.4 acceptable?** I scored with one vcfeval run per pipeline instead
   of the brief's pre-split method, because pre-splitting measurably distorts results
   (*E. coli* BWA+FreeBayes indels: 983 TP / 9 FP / 17 FN pre-split versus 994 / 0 / 6 in a
   single run, with 70.6% of the spurious FNs having another truth variant within 50 bp). It
   penalises FreeBayes hardest, so it changes conclusions rather than shifting all numbers.
   Both are in `results.tsv`. **If the marker expects the brief's method, the primary column
   must change.**

2. **Should the truth sets be regenerated with a realistic Ti/Tv (§2.1)?** It costs ~3 minutes
   to redo the baseline and makes absolute numbers transferable to real data. It must happen
   before the sweep or not at all.

3. **How many seeds are needed to call a difference real?** At baseline the nine *E. coli* SNV
   F1 values span 0.9906–0.9953 — a range of 0.0047. Nothing in this session establishes
   whether that exceeds seed-to-seed noise. If it does not, the headline finding
   ("aligner matters more than caller") is not yet supported, and 5 seeds may not be enough.

---

## 7. Known traps, so you do not rediscover them

| Trap | Symptom |
|---|---|
| Running Snakemake without `--use-conda` | every rule fails with **exit 127** |
| Snakemake cannot find `conda` | `Error running conda info` — looks like a broken install |
| `gatk` called by absolute path | `env: python: No such file or directory` — it is a Python launcher |
| `/usr/bin/time` wrapping a shell function | `time: gatk: No such file or directory` |
| `grep` on simuG's mutated FASTA | silently truncates the single 4.6 Mb line; use `awk` |
| `bcftools ... \| head` under `set -o pipefail` | SIGPIPE kills the whole script silently |
| ART truth SAM vs aligner BAM | different coordinate systems, same contig name — no error, wrong answers |
| RTG output directory already exists | refuses to overwrite; `rm -rf` the target first |
| `grep -c` with `\|\| echo 0` | appends a second line on zero matches |
| macOS `bash` is 3.2 | no associative arrays; `/usr/bin/time` is BSD (`-l`, bytes not KB) |
