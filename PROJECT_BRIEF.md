# PROJECT BRIEF — Read Alignment & Variant-Calling Benchmark

> **How to use this file:** save it as `PROJECT_BRIEF.md` in an empty project directory, then start Claude Code in that directory and say:
> *"Read PROJECT_BRIEF.md and begin. Work through the phases in order. Stop at every GATE and report to me before continuing."*

---

## 0. Your role and my situation

You are acting as a bioinformatics research engineer building a reproducible benchmarking pipeline. I am a student who is **new to DNA sequencing**. I will have to defend every design decision to a professor, and possibly in an interview.

This has two consequences that override normal coding-assistant behaviour:

1. **Explain as you go.** Maintain a running `NOTES.md`. Every time you make a non-obvious choice, add a short entry: what you did, why, and what would have gone wrong otherwise. Write it for someone who does not yet know the field. This file is a deliverable, not a courtesy.
2. **Never fabricate a result.** If a tool fails to install, a command errors, or output looks wrong, **stop and tell me**. Do not substitute a mock, do not generate placeholder numbers, do not write a script that "would work if the tool were installed." A silently wrong benchmark is worse than no benchmark. Report the blocker and wait.

---

## 1. What the project is

We benchmark **3 read aligners × 3 variant callers = 9 pipelines** against genomes where the correct variant list is known by construction.

**The logic:**
1. Take a reference genome.
2. Inject a known set of SNPs and indels into it → this creates a *mutated genome* plus a *truth VCF* recording exactly what was injected, in **reference coordinates**.
3. Simulate sequencing reads **from the mutated genome**.
4. Align those reads **to the original reference**, using each of 3 aligners.
5. Call variants from each alignment, using each of 3 callers → 9 VCFs.
6. Score all 9 against the truth VCF using the GA4GH benchmarking standard (`rtg vcfeval`) → precision, recall, F1.

Later phases (out of scope here) sweep coverage/read length/error rate and fit a model that predicts which pipeline wins under which conditions.

**Aligners:** BWA-MEM, Bowtie2, minimap2 (`-ax sr` preset)
**Callers:** GATK4 HaplotypeCaller, FreeBayes, BCFtools mpileup/call
**Genomes:** phiX174 (`NC_001422.1`, 5,386 bp) and *E. coli* K-12 MG1655 (`NC_000913.3` / `GCF_000005845.2`, 4,641,652 bp). **Both are haploid.**

---

## 2. NON-NEGOTIABLE CORRECTNESS RULES

These are silent-failure modes. Violating any of them produces plausible-looking numbers that are wrong. Verify each explicitly; do not assume.

**R1 — Direction of simulation.** Reads are simulated **FROM the mutated genome** and aligned **TO the original reference**. The truth VCF is in original-reference coordinates. If you ever align to the mutated genome, every pipeline finds zero variants.

**R2 — Ploidy must be set to 1 in every caller.** These organisms are haploid; all three callers default to diploid and will emit heterozygous genotypes with **no error message**.
- GATK: `--sample-ploidy 1`
- FreeBayes: `-p 1`
- BCFtools: `bcftools call --ploidy 1`

After the first run of each caller, **grep the `GT` field and confirm it reads `1`, not `0/1` or `1/1`.** A known failure mode is the flag being dropped in a multi-line shell command. Log the verification result.

**R3 — Normalise everything.** simuG does **not** normalise its output VCF. Run `bcftools norm -f <ref.fa> -m -any` on the truth VCF *and* on every caller's VCF, using identical settings. Skipping this inflates false positives and false negatives, especially for indels.

**R4 — Never hand-roll variant comparison.** The same variant has multiple valid VCF representations. Use `rtg vcfeval` only. Do not write position-and-allele string matching.

**R5 — Contig names must be identical** across reference FASTA, truth VCF, BED, and every call set. Fix this once at download time.

**R6 — Read groups are mandatory.** GATK refuses to run without them. Use `@RG\tID:s1\tSM:sim\tPL:ILLUMINA\tLB:lib1` — via `-R` for BWA-MEM and minimap2, and `--rg-id`/`--rg` for Bowtie2.

**R7 — No BQSR.** GATK Best Practices includes base quality score recalibration, which requires a known-variant database. None exists for *E. coli* or phiX. Skip it and record the reasoning in `NOTES.md`. Do not bootstrap it — that would give GATK a preprocessing step the other two callers don't get, confounding the comparison.

**R8 — Fairness across tools.** Identical thread count for every aligner (runtime is a reported metric). Identical filtering logic across callers. No read trimming — all pipelines must see byte-identical input.

**R9 — Record every seed and every version.** Seeds go in filenames and in the results table. Versions go in `logs/versions_*.txt` on day one.

**R10 — phiX first, always.** phiX is 5,386 bp and runs in seconds. Every phase must work end-to-end on phiX before touching *E. coli*.

---

## 3. Repository layout to create

```
.
├── PROJECT_BRIEF.md          # this file
├── NOTES.md                  # your running explanation log — a deliverable
├── README.md                 # how to reproduce, written last
├── Snakefile                 # workflow (Phase 7b)
├── envs/                     # one conda YAML per environment
├── config/
│   └── conditions.tsv        # sweep design table
├── scripts/                  # helper Python/shell scripts
├── data/
│   ├── refs/                 # downloaded references + indexes
│   └── truth/                # mutated genomes, truth VCFs, confident BEDs
├── work/                     # intermediates — gitignored, deletable
├── results/
│   ├── vcfeval/              # raw scoring output
│   └── results.tsv           # master results table
└── logs/                     # per-run stderr, versions, verification records
```

Initialise git. `.gitignore` must exclude `work/`, `*.fq`, `*.bam`, `*.sdf/`, and large FASTAs.

---

## 4. Phases to execute

Work in order. **Stop at each GATE**, report what you found, and wait for me.

### PHASE 0 — Environment

Install a conda-based package manager if absent (Miniforge preferred). Configure channels:
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

Create **separate** environments — GATK, FreeBayes and hap.py have conflicting dependencies:
```bash
conda create -n align   -c bioconda bwa bowtie2 minimap2 samtools
conda create -n callers -c bioconda gatk4 freebayes bcftools
conda create -n sim     -c bioconda art perl
conda create -n qc      -c bioconda fastqc multiqc
conda create -n ml      -c conda-forge python=3.11 scikit-learn pandas matplotlib seaborn jupyterlab
```

Install **RTG Tools** separately — it is a Java application, not conda-first. Get it from https://github.com/RealTimeGenomics/rtg-tools or https://www.realtimegenomics.com/products/rtg-tools. Requires Java 8+.

Clone **simuG** (single Perl script, no packaging): `git clone https://github.com/yjx1217/simuG.git`

Export every environment to `envs/*.yaml` with pinned versions. Write `logs/versions_*.txt`.

> **GATE 0.** Report a table: each of `bwa`, `bowtie2`, `minimap2`, `samtools`, `gatk`, `freebayes`, `bcftools`, `art_illumina`, `fastqc`, `rtg`, `simuG.pl` — installed yes/no, version, and how you verified it. **If anything fails to install, stop and tell me rather than working around it.** Network restrictions or a Java dependency may need my intervention.

---

### PHASE 1 — Reference genomes

Download phiX174 (`NC_001422.1`) and *E. coli* K-12 MG1655 (`GCF_000005845.2`, RefSeq — **not** GenBank `GCA_`). NCBI Datasets or NCBI Nucleotide.

Normalise contig names to short consistent identifiers and use them everywhere thereafter.

Index both:
```bash
samtools faidx ref.fa
bwa index ref.fa
bowtie2-build ref.fa <prefix>
gatk CreateSequenceDictionary -R ref.fa
rtg format -o ref.sdf ref.fa
```

> **GATE 1.** Confirm from the `.fai`: phiX = 5,386 bp, *E. coli* = 4,641,652 bp. Report the exact accession versions used. If lengths differ, stop — you have the wrong assembly version.

---

### PHASE 2 — Truth sets

Run simuG on each reference:
```bash
perl simuG.pl -refseq <ref.fa> -snp_count <N> -indel_count <M> -seed 20260814 -prefix data/truth/<name>
```
- phiX: 50 SNPs, 10 indels
- *E. coli*: 5,000 SNPs, 1,000 indels (≈1 variant per 930 bp — a realistic strain-divergence scale)

Check the simuG manual for exact output filenames in the cloned version; do not assume.

Then:
1. Concatenate the SNP and INDEL VCFs, sort, and **normalise** (R3).
2. **Handle the genotype column.** simuG's VCF may carry only the eight fixed columns with no sample/GT field, but `rtg vcfeval` compares genotypes. Inspect the actual output, then either (a) add a `FORMAT`/sample column with `GT=1`, or (b) plan to use `vcfeval --squash-ploidy` for allele-level matching. **Test both on phiX, pick one, document the choice in `NOTES.md`, and use it consistently.**
3. Build the confident-regions BED — for simulated data this is the whole genome:
   ```bash
   awk 'BEGIN{OFS="\t"}{print $1, 0, $2}' ref.fa.fai > data/truth/confident.bed
   ```
   In `NOTES.md`, explain *why this file exists* even though it's trivial here: with real benchmark data, truth is known only in part of the genome, and calls outside it must be scored as unknown rather than as false positives.
4. **Verify.** Count records against what was requested. Manually check 3 variants: pull the reference base with `samtools faidx ref.fa <contig>:<pos>-<pos>` and confirm it matches the VCF `REF` field. Report the check.

> **GATE 2.** Report variant counts, the indel length distribution, the genotype-column decision and why, and the results of the manual spot-check.

---

### PHASE 3 — Read simulation

Baseline condition: **30× coverage, 150 bp paired-end, HS25 profile, qs shift 0.**

```bash
art_illumina -ss HS25 -sam \
  -i <mutated_genome.fa> -p -l 150 -f 30 -m 350 -s 50 -rs <seed> \
  -o work/<genome>_cov30_len150_err0_seed<N>_
```

**Keep the `-sam` output.** It records the true origin of every read and enables a read-placement-accuracy metric in Phase 5 that isolates aligner performance from caller performance.

Note in `NOTES.md`: the HS25 profile supports read lengths only up to 150 bp. Reaching 250 bp requires MSv3 (MiSeq), which changes the entire error profile and would confound read length with platform chemistry. The read-length sweep will therefore be restricted to 75/100/150.

Verify read count ≈ requested coverage: `(reads × read_length) / genome_size`.

> **GATE 3.** Report read counts and computed actual coverage for both genomes.

---

### PHASE 4 — QC

Run FastQC on the phiX and *E. coli* baseline FASTQs, aggregate with MultiQC into `results/qc/`.

Report: does per-base quality show the characteristic Illumina 3′ decline? What is mean Q? Do **not** trim (R8).

Write a small script that extracts **mean Q per dataset** into a table — this converts the arbitrary `-qs` shift parameter into a physically interpretable error-rate feature for the later modelling phase.

> **GATE 4.** Report mean Q and confirm the quality profile looks like real Illumina data.

---

### PHASE 5 — Alignment

For each aligner, align the baseline reads to the **original reference** (R1), with read groups (R6), same thread count (R8):

```bash
RG='@RG\tID:s1\tSM:sim\tPL:ILLUMINA\tLB:lib1'
bwa mem -t 4 -R "$RG" ref.fa r1.fq r2.fq | samtools sort -o work/bwa.bam -
bowtie2 -p 4 --rg-id s1 --rg SM:sim --rg PL:ILLUMINA -x <idx> -1 r1.fq -2 r2.fq | samtools sort -o work/bowtie2.bam -
minimap2 -ax sr -t 4 -R "$RG" ref.fa r1.fq r2.fq | samtools sort -o work/minimap2.bam -
samtools index work/*.bam
gatk MarkDuplicates -I in.bam -O out.md.bam -M out.md.metrics
```

Note in `NOTES.md`: simulated reads have no PCR duplicates, so MarkDuplicates will find ~0%. It is run for pipeline realism and its inertness here should be stated in the report.

Collect per-BAM metrics: `samtools flagstat`, `samtools stats`, mean depth, MAPQ distribution, wall-clock time and peak RSS (`/usr/bin/time -v`).

**Write `scripts/placement_accuracy.py`:** parse the ART `-sam` truth file and each aligner's BAM, match reads by name, and compute the fraction placed within ±10 bp of true origin. This is a pure aligner metric, independent of any caller.

> **GATE 5.** Report a table: aligner × (mapping rate, properly-paired rate, mean MAPQ, mean depth, placement accuracy, runtime, peak memory).

---

### PHASE 6 — Variant calling

Run all 3 callers on all 3 BAMs → **9 VCFs**. Ploidy 1 everywhere (R2). No BQSR (R7).

```bash
gatk HaplotypeCaller -R ref.fa -I in.bam -O out.vcf.gz --sample-ploidy 1
freebayes -f ref.fa -p 1 in.bam > out.vcf
bcftools mpileup -f ref.fa -Ou in.bam | bcftools call -mv --ploidy 1 -Oz -o out.vcf.gz
```

**Immediately after the first run of each caller: verify the `GT` field is haploid.** Log the check explicitly to `logs/ploidy_verification.txt`. Do not proceed past a caller emitting diploid genotypes.

Produce **both** a raw and a hard-filtered callset per caller, using identical filtering logic across callers (R8). Retain `QUAL` for ROC analysis.

> **GATE 6.** Report the ploidy verification for all three callers, plus raw variant counts per pipeline. **Do not continue if any caller shows diploid genotypes.**

---

### PHASE 7a — Normalisation and GA4GH scoring

Normalise every call set exactly as the truth set was (R3). Then:

```bash
rtg vcfeval \
  -b data/truth/<truth>.vcf.gz \
  -c <normalised_calls>.vcf.gz \
  -t data/refs/<ref>.sdf \
  -e data/truth/confident.bed \
  --vcf-score-field=QUAL \
  -o results/vcfeval/<condition>_<aligner>_<caller>
```

**Score SNVs and indels separately** (`bcftools view -v snps` / `-v indels`, then vcfeval on each). Pooling them hides the most interesting differences.

Generate ROC curves: `rtg rocplot --svg results/roc.svg results/vcfeval/*/weighted_roc.tsv.gz`

**Write `scripts/collect_results.py`:** walk all vcfeval `summary.txt` files and emit `results/results.tsv` with columns:
```
genome, coverage, read_length, qs_shift, seed, aligner, caller, variant_type,
TP, FP, FN, precision, recall, f1, mapping_rate, mean_mapq, mean_depth,
placement_accuracy, align_seconds, call_seconds, peak_rss_mb
```

Add to `NOTES.md` the metric definitions in plain language, with the worked example: a caller reporting only its single most confident variant achieves precision 1.0 and recall ≈0.0002 — which is why F1 is the headline and why the full ROC curve is also shown.

> **GATE 7a.** Report the 3×3 matrix of F1 scores (separately for SNVs and indels) at baseline, for phiX and *E. coli*. Flag anything anomalous — an F1 of 0.0 or 1.0 usually indicates a bug, not a result.

---

### PHASE 7b — Snakemake workflow

Convert the working shell pipeline into a Snakemake workflow driven by `config/conditions.tsv`, with per-rule conda environments.

Populate `config/conditions.tsv` with the **full** one-factor-at-a-time design (11 conditions per genome), even though only the baseline will be executed now:

| Axis | Levels |
|---|---|
| Coverage | 5, 10, 20, 30\*, 50, 100 |
| Read length | 75, 100, 150\* |
| qs shift | 0\*, −2, −5, −10 |

(\* = baseline, counted once → 6 + 2 + 3 = 11 unique conditions per genome)

Verify with `snakemake -n` that the DAG resolves for the full 11 conditions × 5 seeds × 2 genomes × 9 pipelines. **Execute only the baseline.** Export the DAG image:
```bash
snakemake --dag | dot -Tsvg > results/workflow_dag.svg
```

> **GATE 7b.** Report the dry-run job count and confirm the baseline reproduces the Phase 7a numbers when run through Snakemake.

---

## 5. EXPLICITLY OUT OF SCOPE

Do **not** start these. They are Phases 8–10 and belong to the next work session.

- ❌ Executing the full parameter sweep (~50 hours of compute; needs my sign-off first)
- ❌ Any machine learning, model fitting, or decision-tree work
- ❌ hap.py installation or cross-validation
- ❌ Downloading real read data (ERX008638)
- ❌ GIAB / human data
- ❌ Writing the report or slides
- ❌ Adding aligners or callers beyond the specified 3×3

If you finish early, improve `NOTES.md` and the verification scripts instead.

---

## 6. DEFINITION OF DONE — the ~40% target

This session is complete when **all** of the following are true:

**Infrastructure**
- [ ] Git repo initialised with the layout in §3 and a working `.gitignore`
- [ ] All five conda environments created and exported to `envs/*.yaml` with pinned versions
- [ ] RTG Tools and simuG installed and verified
- [ ] `logs/versions_*.txt` written

**Data**
- [ ] Both reference genomes downloaded, contig-renamed, and fully indexed (fai, bwa, bowtie2, dict, sdf)
- [ ] Truth sets generated, normalised, genotype-column issue resolved and documented
- [ ] Confident-region BEDs created
- [ ] Truth sets manually spot-verified (3 variants each)

**Pipeline validated end-to-end**
- [ ] Baseline reads simulated for both genomes, with ART `-sam` truth retained
- [ ] FastQC/MultiQC report generated; mean-Q extraction script working
- [ ] All 3 aligners producing indexed, duplicate-marked BAMs with read groups
- [ ] Alignment metrics collected, including placement accuracy from the ART truth SAM
- [ ] All 9 aligner×caller combinations producing VCFs
- [ ] **Ploidy verified haploid for all 3 callers, logged**
- [ ] All call sets normalised
- [ ] `rtg vcfeval` scoring working, SNVs and indels separately
- [ ] `results/results.tsv` populated with the baseline (9 pipelines × 2 genomes × 2 variant types = 36 rows)
- [ ] ROC plot generated

**Ready for the sweep**
- [ ] Snakemake workflow reproduces the baseline result
- [ ] `config/conditions.tsv` contains the full 11-condition design
- [ ] `snakemake -n` dry-run resolves cleanly for the full sweep
- [ ] `results/workflow_dag.svg` exported

**Documentation**
- [ ] `NOTES.md` explains every non-obvious decision in beginner-accessible language
- [ ] `README.md` gives one-command reproduction instructions
- [ ] A `HANDOFF.md` listing: what works, what was assumed, anything unverified, estimated compute for the full sweep, and what Phase 8 should do first

---

## 7. Final report to me

When done, give me:
1. The baseline 3×3 F1 matrices (SNV and indel, both genomes)
2. Anything surprising or that looks like a bug rather than a result
3. Every assumption you made that I should verify
4. Measured wall-clock for one full *E. coli* pipeline run, extrapolated to the 990-run sweep
5. The three questions you most want me to answer before Phase 8 begins

---

## 8. Key references

- GA4GH benchmarking standard — Krusche et al. 2019, https://doi.org/10.1038/s41587-019-0054-x
- RTG Tools / vcfeval — https://github.com/RealTimeGenomics/rtg-tools
- simuG — https://github.com/yjx1217/simuG (Yue & Liti 2019)
- ART — https://doi.org/10.1093/bioinformatics/btr708
- GATK HaplotypeCaller — https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
- bcftools — https://samtools.github.io/bcftools/bcftools.html
- Bacterial benchmarking precedent — Bush et al. 2020, https://academic.oup.com/gigascience/article/9/2/giaa007/5728470
- VCF v4.3 spec — https://samtools.github.io/hts-specs/VCFv4.3.pdf
- SAM/BAM spec — https://samtools.github.io/hts-specs/SAMv1.pdf
