# Read Alignment & Variant-Calling Benchmark

Benchmarks **3 read aligners × 3 variant callers = 9 pipelines** against genomes where the
correct variant list is known by construction.

| | |
|---|---|
| **Aligners** | BWA-MEM, Bowtie2, minimap2 (`-ax sr`) |
| **Callers** | GATK4 HaplotypeCaller, FreeBayes, BCFtools mpileup/call |
| **Genomes** | phiX174 (`NC_001422.1`, 5,386 bp) and *E. coli* K-12 MG1655 (`NC_000913.3`, 4,641,652 bp) — both **haploid** |
| **Scoring** | GA4GH standard via `rtg vcfeval` |

**How it works:** inject a known set of SNPs and indels into a reference → simulate reads
**from the mutated genome** → align them back **to the original reference** → call variants →
score against the injected truth. Because we wrote the mutation list ourselves, every call can
be graded as correct, spurious, or missed.

New to the field? **[NOTES.md](NOTES.md)** explains every decision from first principles.
Continuing the work? Start with **[HANDOFF.md](HANDOFF.md)**.

---

## Quick start

### 1. Prerequisites

- macOS or Linux, ~10 GB free disk for the baseline (~50 GB peak for the full sweep, ~20 GB retained)
- **Java 8+** (for RTG Tools): `java -version`
- **XeLaTeX** (e.g. MacTeX or TeX Live), only to render the report as PDF
- A conda-based package manager. If you have none, install
  [Miniforge](https://github.com/conda-forge/miniforge):

```bash
curl -L -o miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && bash miniforge.sh -b -p "$HOME/miniforge3"
```

### 2. One-time setup

Installs the third-party tools that are not conda packages, creates the five conda
environments, and records every version:

```bash
bash setup.sh
```

### 3. Reproduce the study

The workflow's settings (conda, cores, scheduler, fair-timing resources) live in
`profiles/default/`, which Snakemake picks up automatically — so a bare `snakemake` is enough.

**Baseline only** — both genomes at 30×, 150 bp, seed 1 (~5 minutes):

```bash
snakemake
```

**The full study** — 11 conditions × 5 seeds × 2 genomes × 9 pipelines = 990 runs, then the
variance analysis, error-mechanism tests, model, figures and the report, all in one command
(~2.5 hours on a 10-core Apple M4; 14,865 jobs):

```bash
snakemake --config run=all
```

The report is [docs/FINAL_REPORT.pdf](docs/FINAL_REPORT.pdf) (also
[Markdown](docs/FINAL_REPORT.md) and DOCX). Its prose lives in
`docs/report/FINAL_REPORT.template.md`; every table and every quoted number is filled in from
`results/` by `scripts/build_report.py`, so editing the template and re-running the command is
the only way to change the report.

### 4. Check everything is correct

```bash
bash scripts/verify_all.sh
```

Audits every requirement in the definition of done (96 checks at full-sweep scope). It
re-derives facts rather than checking for file presence: reads lengths from the `.fai`, re-greps
every genotype from all 990 VCFs, re-runs REF spot-checks, recomputes the truth-vs-mutated-genome
length bookkeeping, checks the sweep is complete and up to date under Snakemake, and checks that
the report and figures are newer than the results they were built from.

---

## Outputs

| Path | Contents |
|---|---|
| `docs/FINAL_REPORT.pdf` | **the final report** (also `.md`, `.docx`) — every table and number generated from the data |
| `results/results.tsv` | master results table — one row per run × pipeline × variant type × call set × scoring method |
| `results/analysis/` | variance analysis (per-condition summaries, ANOVA/Friedman effects), secondary experiments, error-mechanism tests (`fp_near_indel`, `fn_repeats`, `phix_errors`) |
| `results/model/` | decision tree, validation, permutation importances, recommendations |
| `results/figures/` | report figures F1–F7 |
| `results/read_metrics.tsv` | measured coverage and error rate for every read set |
| `results/runtime.tsv` | wall time, peak RSS and CPU per timed job (`exclusive = yes` rows are clean) |
| `results/roc_*.svg` | ROC curves, 9 pipelines each, SNV and indel separately |
| `results/align_metrics.tsv` | mapping rate, MAPQ, depth, runtime, peak RSS per aligner |
| `results/placement_accuracy.tsv` | pure aligner metric from the ART truth SAM |
| `results/qc/multiqc_report.html` | aggregated read QC |
| `results/workflow_rulegraph.svg` | pipeline structure (the readable graph) |
| `logs/ploidy_verification.txt` | R2 evidence — every genotype confirmed haploid |
| `logs/verification_report.txt` | full Definition-of-Done audit |

### Reading `results.tsv`

Filter on **`callset == raw`** for the primary analysis; `filt` rows are the hard-filtered call
sets (`QUAL >= 20 && DP >= 5`). **`timing_exclusive == yes`** marks runs whose timings were taken
with the machine to itself (seed 1); other timings were made under contention and are not used for
runtime comparisons.

The **`scoring_method`** column matters. Use **`single_run`** — these are the primary numbers,
taken from one `vcfeval` run per pipeline using RTG's own SNV/indel breakdown. The
`pre_split` rows split the VCF by type *before* scoring, which degrades vcfeval's
haplotype-aware comparison and understates performance (especially FreeBayes). They are kept
only so the difference is auditable. See [NOTES.md §7.6](NOTES.md).

---

## The parameter sweep

`config/conditions.tsv` defines a one-factor-at-a-time design:

| Axis | Levels |
|---|---|
| Coverage | 5, 10, 20, **30\***, 50, 100 |
| Read length | 75, 100, **150\*** |
| qs shift | **0\***, −2, −5, −10 |

\* = baseline. 11 unique conditions per genome × 5 seeds × 2 genomes × 9 pipelines =
**990 pipeline runs**.

```bash
snakemake -n --config run=all      # dry run
snakemake --config run=all         # execute
snakemake --config run=tags tags=TAG1,TAG2   # any subset, e.g. a pilot
```

The full sweep has been executed; see the report. `HANDOFF.md` records the state of the project
and what Phase 2 should do first.

---

## Repository layout

```
├── Snakefile                  workflow definition
├── setup.sh                   one-time tool + environment install
├── config/
│   ├── config.yaml            workflow settings (threads, filters, read group)
│   └── conditions.tsv         the 11-condition sweep design
├── envs/                      conda environments
│   ├── <name>.yaml            portable: top-level packages, pinned
│   └── <name>.lock.yaml       exact: full URL lock of the exact builds used
├── scripts/                   pipeline steps + verification tools
├── data/refs/                 references + all five index types
├── data/truth/                mutated genomes, truth VCFs, confident BEDs
├── docs/                      final report (prose template in docs/report/), slides, plans
├── work/                      intermediates (gitignored, deletable)
├── results/                   results table, analysis, model, figures, ROC curves, QC
└── logs/                      per-run stderr, versions, verification records
```

`work/`, `.snakemake/` and `tools/` are gitignored — all are regenerable, and together they
run to several GB.

---

## Correctness rules

The benchmark enforces ten rules whose violation produces *plausible-looking but wrong*
numbers. Each is explained in [NOTES.md](NOTES.md); the ones with the sharpest teeth:

- **R1** — reads are simulated FROM the mutated genome and aligned TO the original reference.
  `simulate_reads.sh` and `align_reads.sh` both refuse to run if this is inverted.
- **R2** — ploidy 1 in every caller. These organisms are haploid; all three callers default to
  diploid and emit heterozygous genotypes **with no error message**. Verified by grepping the
  `GT` field, not by inspecting scores — a caller emitting `1/1` scores a *perfect* F1 against
  haploid truth, so the score cannot detect the mistake.
- **R3** — truth and every call set are normalised with byte-identical arguments, including
  `--atomize`.
- **R4** — variant comparison is done only by `rtg vcfeval`. Never by string matching.
- **R5** — contig names are normalised once at download time. A mismatch does not error; it
  silently reports precision 0 / recall 0.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| Every rule fails with **exit 127** | Ran without `--use-conda` (the default profile sets it — check you are in the repo root) |
| Timings look noisy | Only `timing_exclusive == yes` rows were measured without contention |
| `collect_results` fails with "empty metric cells" | A metric join failed — this is `--strict` doing its job; check `results/*.tsv` |
| `Error running conda info` | Snakemake cannot find `conda` — add it to `PATH` |
| `env: python: No such file or directory` | `gatk` is a Python launcher and needs its env's `bin/` on `PATH` |
| `rtg ... already exists. Please remove it` | RTG refuses to overwrite output; delete the target first |
| `AmbiguousRuleException` | Snakemake wildcards match `/`, so a "more specific" path rule is ambiguous, not nested |

---

## References

- GA4GH benchmarking standard — [Krusche et al. 2019](https://doi.org/10.1038/s41587-019-0054-x)
- [RTG Tools / vcfeval](https://github.com/RealTimeGenomics/rtg-tools)
- [simuG](https://github.com/yjx1217/simuG) — Yue & Liti 2019
- ART read simulator — [Huang et al. 2012](https://doi.org/10.1093/bioinformatics/btr708)
- Bacterial benchmarking precedent — [Bush et al. 2020](https://academic.oup.com/gigascience/article/9/2/giaa007/5728470)
