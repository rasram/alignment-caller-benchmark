# NOTES — running explanation log

This file explains **every non-obvious decision** made while building the benchmark,
written for someone who is new to DNA sequencing. Each entry says *what* was done,
*why*, and *what would have gone wrong otherwise*.

Entries are appended in phase order. Newest phase at the bottom.

---

## Background: what this project actually measures, in plain language

A **sequencing machine** does not read a genome end to end. It shreds many copies of the
DNA into short fragments and reads a few hundred letters from each fragment. Those short
strings are called **reads**. A typical read here is 150 letters long, and we generate
enough of them that every position in the genome is covered by ~30 different reads
("30× coverage").

Turning reads back into biology takes two steps, and this project benchmarks both:

1. **Alignment.** Each read is a 150-letter fragment from *somewhere* in a ~4.6-million-letter
   genome. An **aligner** finds where. This is hard because the read has sequencing errors
   in it, and because the genome contains repeated sequence that a short read may match in
   several places. We compare three aligners: BWA-MEM, Bowtie2, minimap2.
2. **Variant calling.** Once reads are placed, a **caller** looks at each genome position,
   sees the pile of reads stacked there, and decides whether the sample genuinely differs
   from the reference at that position — or whether the apparent difference is just
   sequencing error. We compare three callers: GATK HaplotypeCaller, FreeBayes, BCFtools.

3 aligners × 3 callers = **9 pipelines**, all scored on identical input.

### Why simulated data, and why that is the whole point

To grade a pipeline you must already know the right answer. For a real biological sample
nobody does. So we manufacture a sample where the answer is known **by construction**:

- Take the real reference genome.
- Deliberately inject a known list of mutations into it → a **mutated genome**, plus a
  **truth VCF** that records exactly what was injected.
- Simulate reads **from the mutated genome**.
- Align those reads back **to the original, unmutated reference**.

The pipeline's job is to rediscover, from the reads alone, the mutation list we injected.
Because we wrote that list ourselves, every call the pipeline makes can be graded as
correct (true positive), spurious (false positive), or missed (false negative).

The direction matters enormously and is the single easiest thing to get backwards
(rule R1 in the brief). Reads come **from the mutated genome**; they are aligned **to the
original reference**. If you aligned reads to the mutated genome instead, the reads would
match it perfectly, every pipeline would report zero variants, and all nine pipelines
would score identically and meaninglessly.

### What a "variant" is here

- **SNP / SNV** — a single-letter substitution: reference has `A`, sample has `G`.
- **Indel** — an insertion or deletion of one or more letters. Indels are much harder for
  both aligners and callers than SNPs, which is why the brief insists on scoring the two
  categories separately: pooling them lets good SNP performance mask bad indel performance.

### Why phiX and *E. coli*

- **phiX174** — a virus, 5,386 letters. Runs end to end in seconds, so it is the debugging
  organism. Rule R10: every phase must work on phiX before *E. coli* is touched.
- ***E. coli* K-12 MG1655** — a bacterium, 4,641,652 letters (~860× larger). Big enough to
  have real repeated sequence, so aligners can actually be distinguished from one another.

Both are **haploid** — one copy of the genome, not two. This drives correctness rule R2 and
is the most dangerous silent failure in the whole project (see the Phase 6 entry).

---

## Phase 0 — Environment

### 0.1 Why conda, and why five separate environments

Bioinformatics tools are old, numerous, and built against mutually incompatible versions of
the same libraries. `conda` installs each tool with its own pinned dependency tree.

They are split across five environments because **some of these tools genuinely cannot
coexist**. GATK4, FreeBayes and hap.py have conflicting dependency constraints; asking conda
to satisfy all of them at once either fails outright or silently downgrades something to an
old build. Splitting them means each environment gets the versions its tools actually want.

The five:

| Env | Contents | Role |
|---|---|---|
| `align` | bwa, bowtie2, minimap2, samtools | the three aligners |
| `callers` | gatk4, freebayes, bcftools | the three callers |
| `sim` | art, perl | read simulation + simuG (a Perl script) |
| `qc` | fastqc, multiqc | read quality control |
| `ml` | python 3.11, pandas, scikit-learn, matplotlib, seaborn, jupyterlab, snakemake, graphviz | analysis + workflow engine |

**What would have gone wrong otherwise:** one big environment resolves to whatever version
set happens to satisfy everything, which typically means an outdated GATK. The benchmark
would then be measuring a version of the tool nobody uses.

### 0.2 Miniforge rather than Anaconda

Installed **Miniforge3** (arm64 native). Miniforge defaults to the `conda-forge` channel and
carries no Anaconda commercial licence terms. Anaconda's `defaults` channel now requires
accepting Terms of Service for institutional use, which would make the pipeline
non-reproducible for anyone who has not clicked through that agreement.

After configuring channels I explicitly **removed** `defaults`, `repo.anaconda.com/pkgs/main`
and `pkgs/r` from `.condarc`, leaving exactly:

```yaml
channels:
  - bioconda
  - conda-forge
channel_priority: strict
```

`channel_priority: strict` means a package is taken from the highest-priority channel that
has it, and conda will not mix-and-match builds across channels to satisfy a solve. Without
`strict`, conda can assemble a working-looking environment from packages compiled against
different library versions — which produces crashes or, worse, silently wrong numeric output.
Bioconda is listed **above** conda-forge because that is the ordering bioconda itself
documents and tests against.

### 0.3 Apple Silicon: checked before assuming

This machine is an **arm64 Mac (macOS 26.6)**. Much bioinformatics software is distributed
only as x86-64, and the usual workaround is to force conda into `osx-64` mode and run every
tool under Rosetta 2 emulation — which is slower, and **runtime is a reported metric in this
benchmark**, so emulation would have contaminated the results.

So rather than assume, I queried the channel for every required package on both
architectures before creating anything. All ten (`bwa`, `bowtie2`, `minimap2`, `samtools`,
`gatk4`, `freebayes`, `bcftools`, `art`, `fastqc`, `multiqc`) have native `osx-arm64` builds.
**No emulation is used anywhere**, so the reported runtimes are native and comparable across
the three aligners.

This is worth stating explicitly in the report: had any single aligner required Rosetta while
the others ran native, the runtime comparison in Phase 5 would have been meaningless.

### 0.4 RTG Tools installed outside conda

`rtg vcfeval` is the GA4GH-standard variant comparison tool (rule R4). It is a **Java
application**, distributed as a zip from GitHub rather than as a conda-first package, so it
is installed into `tools/` and invoked by absolute path.

The release carries no macOS bundle, only Linux and Windows JRE bundles plus a **`nojre`**
zip that uses the system Java. This machine has **Java 23**; RTG requires Java 8+. Verified
by running `rtg version` — it reports cleanly under Java 23.

I also wrote a `rtg.cfg` next to the binary setting `RTG_TALKBACK=false` and
`RTG_USAGE=false`. By default RTG phones home with usage statistics and crash reports. Off
for reproducibility and because the pipeline should not depend on network availability at
scoring time.

### 0.5 simuG cloned, not installed

`simuG` is a single Perl script with no packaging, so it is cloned into `tools/simuG` and run
as `perl tools/simuG/simuG.pl`. The `sim` environment carries `perl` so the interpreter
version is pinned along with everything else rather than depending on whatever
`/usr/bin/perl` the OS ships.

### 0.6 A pre-existing Anaconda on this machine, deliberately not used

This Mac already has **Anaconda 25.5.1 at `/opt/anaconda3`**, but its shell-init block is
commented out in `~/.zshrc`, so it is not on `PATH`. Every script in `scripts/` therefore
hardcodes `CONDA_BASE="$HOME/miniforge3"` rather than calling a bare `conda`.

This matters for reproducibility: if the two installs were ever both on `PATH`, `conda run -n
callers` could resolve to a *different* `callers` environment depending on shell startup
order, and the benchmark could silently be run with different tool versions than the ones
recorded in `logs/versions_*.txt`.

### 0.7 The network here is intermittently unreliable — retries are load-bearing

The first environment build failed for 4 of 5 environments with
`Download error (6) Could not resolve hostname`. This was **not** a missing-package problem:
`align` had already installed cleanly, and the same URLs resolved fine minutes earlier and
minutes later. The system resolver on this machine points at `127.0.2.2` (a local DNS proxy —
VPN or filtering software), which appears to drop requests intermittently.

Rather than paper over it, conda was reconfigured to be patient:

```
remote_connect_timeout_secs 60
remote_read_timeout_secs    300
remote_max_retries          10
remote_backoff_factor       5
```

and each environment build is wrapped in an outer retry loop.

That fixed the download timeouts but exposed a second, different failure — worth writing
down because the error message points nowhere near the actual cause:

```
Failed to load subdir: Download error (7) Could not connect to server
    [https://prefix.dev/conda-forge/osx-arm64/repodata.json.zst]
critical libmamba Could not solve for environment specs
```

"Could not solve for environment specs" normally means *you asked for an impossible
combination of packages*. Here it meant nothing of the sort. Miniforge's base config at
`$HOME/miniforge3/.condarc` ships a mirror list:

```yaml
mirrored_channels:
  conda-forge:
    - https://conda.anaconda.org/conda-forge
    - https://prefix.dev/conda-forge
```

`prefix.dev` is unreachable from this host, and **mamba 2.x treats one dead mirror as a dead
subdir** rather than falling back to the mirror that works. So it never loaded the conda-forge
package list at all, and then reported that as an unsatisfiable dependency problem.

Fix: delete the `mirrored_channels` block so all fetches go to `conda.anaconda.org` directly.

**Why this is in NOTES rather than just fixed silently:** the misleading error would have sent
a reasonable person off editing package specs — downgrading Python, dropping `snakemake`,
relaxing versions — trying to satisfy a solver that was never actually solving. Any of those
"fixes" would have produced a working-but-wrong environment. The lesson generalises: when
conda reports an unsolvable environment, **check that it successfully downloaded the channel
metadata first**.

**Worth flagging for the sweep:** the full Phase 8 sweep involves no downloads, so this
affects setup only. But anything that fetches at runtime (RTG's usage reporting, which is why
`RTG_TALKBACK=false` was set in §0.4) is a liability on this host.

### 0.8 Verification was functional, not just `--version`

GATE 0 asks *how* each tool was verified. Printing a version string only proves a file exists
and is executable. So every tool was additionally exercised on a throwaway 2-contig, 50 bp
FASTA: `samtools faidx`, `bwa index`, `bowtie2-build`, `minimap2 -d`, `rtg format`,
`art_illumina` (simulated 8 reads), `gatk CreateSequenceDictionary`, and `simuG.pl`
(generated 3 SNPs). All succeeded and produced non-empty, correct-looking output.

This caught two things a version check would have missed, both recorded below.

### 0.9 GATK runs on Java 25, which is *not* its supported version

The `callers` environment resolved **OpenJDK 25** as GATK4's Java dependency. GATK 4.6 is
officially tested against **Java 17**. Nobody warns you about this; conda simply installed
the newest JDK satisfying the recipe.

It works so far — `gatk CreateSequenceDictionary` completed and produced a correct `.dict`
with proper `M5` checksums, exercising real HTSJDK and Picard code paths.

**But this is an open risk, not a cleared one.** `HaplotypeCaller` is far heavier than
`CreateSequenceDictionary`: it uses the Intel GKL native libraries for pair-HMM, and JVM
major-version mismatches have historically broken exactly that. If Phase 6 shows GATK
crashing or running suspiciously slowly, the first thing to try is pinning `openjdk=17` in
the `callers` environment. Flagged in `HANDOFF.md`.

Note also the *two different Javas* in play: system Java 23 runs RTG Tools, env-local Java 25
runs GATK. That is fine — they are isolated — but it means "what Java am I on?" has no single
answer here, so `logs/versions_*.txt` records both.

### 0.10 bcftools is 1.21 in `callers`, samtools is 1.24 in `align`

`bcftools` 1.24 exists in the channel, but the `callers` solve settled on **1.21** because of
htslib constraints shared with `gatk4` and `freebayes`.

This is acceptable but must be handled deliberately. `bcftools norm` is the normalisation step
mandated by rule R3, and R3's requirement is that the truth VCF and *every* caller's VCF are
normalised with **identical settings** — which implies an identical binary. So all
normalisation in this project is done with the single `bcftools` from the `callers`
environment (1.21), never the one that might arrive via `samtools`' env. If normalisation were
split across two bcftools versions, indel left-alignment could differ subtly between truth and
calls and would show up as phantom false positives that no amount of staring at the VCFs would
explain.

### 0.11 `tools/` is gitignored

`tools/` holds third-party software (RTG's ~9 MB jar, the simuG clone). Committing other
projects' binaries into this repo bloats history and confuses provenance. The README gives
the two commands that recreate `tools/` from scratch, and the exact versions are pinned in
`logs/versions_*.txt`.

---

## Phase 1 — Reference genomes

### 1.1 What a "reference genome" is, and why we align to it

The reference is one agreed-upon consensus sequence for a species, written down once so that
everybody can describe positions in the same coordinate system. "Position 1,234,567" only
means something relative to a named reference. That is why rule R5 (identical contig names
everywhere) matters so much: a coordinate without an agreed contig name is meaningless.

### 1.2 Accessions used — RefSeq, not GenBank

| Genome | Accession | Assembly | Length (verified from `.fai`) |
|---|---|---|---|
| phiX174 | `NC_001422.1` | — | **5,386 bp** ✓ |
| *E. coli* K-12 MG1655 | `NC_000913.3` | `GCF_000005845.2` | **4,641,652 bp** ✓ |

The original NCBI headers are preserved verbatim in `data/refs/*.original_header.txt`, because
the renaming step below destroys them and provenance must remain auditable:

```
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
>NC_001422.1 Escherichia phage phiX174, complete genome
```

`NC_` accessions are **RefSeq** — NCBI's curated copy. The brief requires RefSeq (`GCF_`) over
GenBank (`GCA_`). They are usually the same sequence, but GenBank is the depositor's original
submission and RefSeq is the maintained version; they can and do diverge across updates. Mixing
them across a project is a classic source of off-by-a-few coordinate errors.

Downloaded via NCBI E-utilities `efetch`, wrapped in a retry loop that also asserts the body
starts with `>`. Without that check, a rate-limited empty response would have been written out
as a plausible-looking but truncated FASTA.

### 1.3 Contig renaming — done once, at download time (R5)

NCBI FASTA headers carry a description after the accession:

```
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
```

Two hazards. First, tools disagree about that description — some truncate the name at the
first whitespace, some keep the whole line — so the "same" contig can end up named
`NC_000913.3` in one file and `NC_000913.3 Escherichia coli...` in another. Second,
versioned accessions are easy to mistype inconsistently across dozens of commands.

So the name is normalised **once, here**, to a short token: `phiX` and `ecoli`. Everything
downstream — mutated genome, truth VCF, confident BED, all nine call sets, the RTG SDF —
inherits it.

**What would have gone wrong otherwise:** `rtg vcfeval` matches truth to calls by contig name.
If truth says `NC_000913.3` and the calls say `ecoli`, vcfeval does not error. It finds zero
overlapping records and reports **precision 0, recall 0** — which looks exactly like a
catastrophically bad pipeline rather than a naming bug. This is the single most likely way to
waste a day on this project.

### 1.4 The five indexes, and why each exists

An index is a precomputed data structure that lets a tool find a sequence without scanning the
whole genome. Each tool wants its own format, and they are **not** interchangeable:

| Command | Produces | Consumed by |
|---|---|---|
| `samtools faidx` | `.fa.fai` | random access to reference bases; also the source for the confident BED |
| `bwa index` | `.amb .ann .bwt .pac .sa` | BWA-MEM (FM-index / Burrows-Wheeler) |
| `bowtie2-build` | 6 × `.bt2` | Bowtie2 (its own FM-index) |
| `gatk CreateSequenceDictionary` | `.dict` | GATK — refuses to run without it |
| `rtg format` | `.sdf/` directory | `rtg vcfeval` |

All five were built **from the same renamed FASTA**, so all agree on the contig name. Verified
explicitly: `.fai`, `.dict` and the SDF all report `phiX`/`ecoli` with the correct lengths.

### 1.5 Sanity checks beyond the length assertion

Matching lengths proves the right assembly version, but not that the file contains real
sequence. Two cheap extra checks:

- **GC content** — phiX **44.76%** (published ~44.8%), *E. coli* K-12 **50.79%**
  (published ~50.8%). Both match, so these are the genomes they claim to be.
- **Ambiguous bases** — **zero** non-ACGT characters in either genome. Worth knowing: runs of
  `N` in a reference create regions where no caller can call anything, which would otherwise
  show up later as unexplained false negatives. There are none here, so any FN in Phase 7a is
  a real pipeline limitation, not a masked-reference artefact.

---

## Phase 2 — Truth sets

### 2.1 What simuG gives you, and why it is not yet usable

`simuG` takes a reference and injects mutations, writing four files per run:

| File | Contents |
|---|---|
| `<p>.simseq.genome.fa` | the **mutated genome** — reads get simulated from this |
| `<p>.refseq2simseq.SNP.vcf` | the injected SNPs, in **reference** coordinates |
| `<p>.refseq2simseq.INDEL.vcf` | the injected indels, in **reference** coordinates |
| `<p>.refseq2simseq.map.txt` | reference ↔ mutated coordinate correspondence |

The brief says not to assume these names — verified against the cloned version (`0289e58`)
before writing any script against them.

"In reference coordinates" is the crucial property and the thing that makes rule R1 work. The
VCF describes the mutations *as offsets into the original reference*, which is the same
coordinate system the aligner will place reads into. If simuG had reported mutated-genome
coordinates, the truth set would be shifted by the cumulative indel length at every position
and nothing would line up.

Four things were wrong with the raw output for our purposes:

1. **Two separate files** (SNP and INDEL) — vcfeval scores one baseline.
2. **No `##contig` header lines** — verified: zero present. Downstream tools then have no
   declared contig order and either refuse to run or sort records unpredictably.
3. **No `FORMAT` column and no sample** — only the 8 fixed columns. See §2.3.
4. **Not normalised** — rule R3. See §2.2.

### 2.2 Normalisation, and hard evidence that it was needed

`bcftools norm -f <ref.fa> -m -any` was run on the truth set (and will be run identically on
every call set). Two operations:

- **`-f ref` left-align and trim.** The same indel can be written at several positions when it
  sits in a repeat. Deleting one `T` from `TTTT` can be described as deleting the T at any of
  four positions — all equally true, all different VCF records. Left-alignment forces one
  canonical choice: push it as far left as possible.
- **`-m -any` split multi-allelics.** One record listing two ALT alleles becomes two records.

**This was not theoretical.** Normalisation changed real records:

| Genome | Records | Realigned by `norm` |
|---|---|---|
| phiX | 60 | **2** |
| *E. coli* | 6,000 | **256** (25.6% of the 1,000 indels) |

Concretely, on phiX:

```
simuG wrote:      582  T    -> TT          2002  TTT -> T
normalised to:    581  C    -> CT          2000  GTT -> G
```

Both moved left into a homopolymer run. Had the truth kept simuG's coordinates while the
callers emitted normalised ones, those variants would have been counted as **one false
negative plus one false positive each** — a double penalty for a variant the pipeline found
perfectly well. At *E. coli* scale that is 256 indels, which would have made every pipeline's
indel recall look ~25% worse than reality.

### 2.3 The genotype column — decision and evidence

**Decision: option (a).** The truth VCF carries `FORMAT=GT` and a single haploid sample
`sim` with **`GT=1`**, and scoring uses vcfeval's **default genotype-aware matching — no
`--squash-ploidy`**.

The experiment is reproducible via `scripts/test_genotype_encoding.sh`; results in
`logs/genotype_decision.txt`. Synthetic call sets were built *from the truth itself*, so
variant content is identical everywhere and only the GT encoding differs — any departure from
F1 = 1.0 is therefore caused purely by genotype representation.

| Truth | Call GT | `--squash-ploidy` | Precision | Recall | F1 |
|---|---|---|---|---|---|
| GT=1 | `1` | no | 1.0000 | 1.0000 | **1.0000** |
| GT=1 | `1/1` | no | 1.0000 | 1.0000 | **1.0000** |
| GT=1 | `0/1` | no | 0.0000 | 0.0000 | **0.0000** |
| GT=1 | `1/1` | yes | 1.0000 | 1.0000 | 1.0000 |
| GT=1 | `0/1` | yes | 1.0000 | 1.0000 | 1.0000 |
| *no sample column* | `1` | either | — | — | **vcfeval refuses to run** |

Three findings, in order of importance:

**1. Option (b) does not exist.** vcfeval rejects a baseline with no sample column outright:
`Error: Record did not contain enough samples`. `--squash-ploidy` does not rescue it — that
flag relaxes *genotype* comparison, it does not conjure a genotype that was never there. So
the FORMAT/sample column is mandatory, and the choice collapses to "how should we score",
not "should we add GT".

**2. `--squash-ploidy` would hide the exact failure R2 warns about.** A call set with
heterozygous `0/1` genotypes scores **F1 = 0.0000** under default matching and **F1 = 1.0000**
with squashing. Rule R2 exists because all three callers default to diploid and will silently
emit het genotypes on a haploid organism. Turning on `--squash-ploidy` would make that
misconfiguration invisible and report a perfect score. Default matching makes it deafening.
That is why squashing is **off**.

**3. vcfeval will not catch every ploidy error by itself — hence R2's grep.** Note row two:
`GT=1/1` matches haploid truth *perfectly* without squashing, because vcfeval treats
homozygous-ALT as equivalent to haploid-ALT. So a caller emitting `1/1` everywhere would score
1.0 and look fine while being configured wrongly. **The score cannot be used as a ploidy
check.** This is exactly why R2 demands grepping the `GT` field directly in Phase 6 rather
than inferring correctness from F1.

`GT=1` — a single allele index with no slash — is the correct haploid encoding. `1/1` means
"diploid, both copies ALT" and `0/1` means "diploid, one copy each", and neither is a true
statement about an organism with one genome copy.

### 2.4 The confident-regions BED, and why a trivial file is not pointless

```
phiX    0   5386
ecoli   0   4641652
```

One line per contig, spanning the whole genome. Built from the `.fai` so contig names are
inherited, never retyped (R5). BED is **0-based half-open** while VCF is **1-based inclusive**
— the `0` is not an off-by-one error, it is a different coordinate convention, and mixing the
two is the classic bug in this area.

**Why it exists even though it covers everything.** vcfeval sorts calls into three bins, not
two: true positive, false positive, and *ignored*. The confident-region BED decides which
calls are eligible to be judged at all. With real benchmark data (GIAB and similar) truth is
only established in part of the genome — the rest is repetitive or structurally messy — and a
call outside those regions must be scored as **unknown**, not as a false positive, because
nobody knows whether it is right.

Here, truth is known *everywhere* by construction: we wrote the mutations ourselves, so any
call not in the truth set genuinely is wrong. The whole genome is confident.

Keeping the file anyway means (i) the command line is identical to what a real benchmark uses,
so nothing has to change when someone later swaps in GIAB data, and (ii) the scoring is
explicit about its scope rather than relying on a default. It also matters for the Phase 8
sweep: if a later condition needs regions masked out, the mechanism is already wired in.

### 2.5 Verification — three independent checks

`scripts/verify_truth.sh` runs three checks that fail in different ways, logged to
`logs/truth_verification_*.txt`:

1. **Counts** — do we have what we asked simuG for? phiX 50/10, *E. coli* 5000/1000. Both exact.
2. **REF bases** — for 3 SNPs and 3 indels per genome, pull the reference base with
   `samtools faidx` and compare to the VCF `REF` field. All 12 matched. This catches
   coordinate-system errors, which are otherwise invisible.
3. **Length bookkeeping** — the strongest check, and one the brief did not ask for.

The third deserves explanation. Sum the length change of every indel in the truth VCF
(`len(ALT) - len(REF)`) and compare it to the actual size difference between the mutated
genome and the reference:

| Genome | Reference | Mutated | Observed Δ | Net indel sum in VCF | |
|---|---|---|---|---|---|
| phiX | 5,386 | 5,400 | +14 | +14 | ✓ |
| *E. coli* | 4,641,652 | 4,641,839 | +187 | +187 | ✓ |

This ties the truth VCF to *the actual FASTA the Phase 3 reads are generated from*. Checks 1
and 2 both pass even if simuG's VCF and its mutated genome disagree — and if they disagreed,
every number in the entire benchmark would be wrong, with no other symptom than mysteriously
poor scores across all nine pipelines. Now confirmed to the base.

### 2.6 Properties of the truth sets, and one caveat to flag

| | phiX | *E. coli* |
|---|---|---|
| SNPs | 50 | 5,000 |
| Indels | 10 | 1,000 |
| Density | 1 per 89 bp | 1 per 773 bp |
| Insertions : deletions | 4 : 6 | 508 : 492 |
| Indel size range | −2 to +13 | −45 to +50 |
| Ti/Tv | 0.667 | 0.480 |

**The Ti/Tv caveat.** Transitions (A↔G, C↔T) and transversions (everything else) do not occur
equally in real genomes — real bacterial genomes run roughly Ti/Tv ≈ 1–2 because transitions
are chemically easier. Our truth sets sit at **~0.5**, which is simuG's default and is exactly
what uniformly-random base substitution produces: each base has one transition partner and two
transversion partners, so random picking gives 1:2 = 0.5.

So **these truth sets are less realistic than real strain divergence in their mutation
spectrum.** It does not bias the comparison — all nine pipelines are scored against the same
truth, and no aligner or caller here is tuned for a particular Ti/Tv — so the *ranking* is
unaffected. But absolute recall numbers are not directly transferable to a real resequencing
project. Setting `-titv_ratio 2.0` would fix it; the brief specifies neither, so simuG's
default was kept and this is flagged as an assumption to confirm.

Indel sizes reach ±50 bp, which is genuinely hard for 150 bp reads — a 50 bp deletion leaves
only ~50 bp of anchor on each side. Expect indel recall to fall off sharply with size, and
expect that to be one of the more interesting aligner differences.

---

## Phase 3 — Read simulation

### 3.1 What ART does, and the direction that must not be reversed

ART takes a genome and produces the FASTQ files a real Illumina machine would have produced
from it: fragments of a chosen length, read from both ends, with realistic position-dependent
sequencing errors and quality scores drawn from an empirical profile of real instrument data.

**Rule R1 lives here.** The input is `data/truth/<g>.simseq.genome.fa` — the **mutated**
genome — never `data/refs/<g>.fa`. `scripts/simulate_reads.sh` asserts this: it refuses to run
if the input file is byte-identical to the reference. That guard is cheap insurance against
the one mistake that would silently make every pipeline score zero.

### 3.2 Baseline parameters

```
art_illumina -ss HS25 -sam -na -i <mutated>.fa -p -l 150 -f 30 -m 350 -s 50 \
             -qs 0 -qs2 0 -rs <seed> -o work/<tag>_
```

| Flag | Meaning |
|---|---|
| `-ss HS25` | HiSeq 2500 empirical error profile |
| `-p -l 150` | paired-end, 150 bp per read |
| `-f 30` | 30× fold coverage |
| `-m 350 -s 50` | DNA fragment length: mean 350 bp, sd 50 |
| `-qs 0 -qs2 0` | quality-score shift for R1/R2 — 0 leaves the profile untouched |
| `-rs <seed>` | random seed (R9) |
| `-sam` | **also write the truth alignment** — see §3.4 |
| `-na` | skip the `.aln` files; the SAM carries the same information |

**Fragment length vs read length.** A 350 bp fragment sequenced 150 bp from each end leaves a
~50 bp unsequenced gap in the middle. This is normal and is what makes paired-end data useful:
the aligner knows the two reads should land ~350 bp apart in the correct orientation, which
helps place reads that are individually ambiguous.

Verified against what ART actually produced (TLEN in the truth SAM):

| Genome | pairs | mean fragment | sd | min | max |
|---|---|---|---|---|---|
| phiX | 540 | 347.7 | 49.4 | 191 | 491 |
| *E. coli* | 464,175 | 349.5 | 50.0 | 150 | 581 |

Both match the requested 350 ± 50.

### 3.3 Coverage verification

"30× coverage" means each base is covered by ~30 reads on average:
`(number of reads × read length) / genome size`.

Measured against the **mutated** genome length, because that is the template the reads were
drawn from. Using the reference length would be wrong by the net indel balance — negligible
here (187 bp in 4.6 Mb) but wrong on principle, and it would grow if a future condition
injected more indels.

| Genome | Mutated length | Pairs | Reads | Actual coverage |
|---|---|---|---|---|
| phiX | 5,400 | 540 | 1,080 | **30.00×** |
| *E. coli* | 4,641,839 | 464,175 | 928,350 | **30.00×** |

Coverage is an *average*. Real per-base depth follows roughly a Poisson distribution around
30, so some positions get 15× and some 45× purely by chance. That variance is exactly why low
coverage hurts variant calling: at 5× a position can easily receive 1–2 reads, and no caller
can distinguish a real variant from a sequencing error with that little evidence.

### 3.4 Why the `-sam` truth file is retained — and a coordinate trap in it

ART's `-sam` output records, for every read, **where in the input genome that read actually
came from**. This is ground truth no real experiment ever has, and it enables a metric in
Phase 5 that isolates the aligner completely: *what fraction of reads did the aligner put back
where they belong?* That is a pure aligner property, measurable without running any variant
caller, so it cleanly separates "BWA placed reads better" from "GATK called variants better".

**The trap.** The truth SAM is in **mutated-genome coordinates**:

```
@SQ  SN:phiX   LN:5400        <- ART truth SAM (mutated genome)
@SQ  SN:phiX   LN:5386        <- aligner BAMs in Phase 5 (original reference)
```

The contig is named `phiX` in both, and the lengths differ by only the net indel balance. So
comparing the two directly **produces no error, no warning — just wrong answers.** Positions
drift apart as you move along the genome, accumulating the indel offset seen so far. For
*E. coli* the drift reaches 187 bp by the end of the genome, far beyond the ±10 bp tolerance
the placement metric uses, so reads near the end would nearly all be scored as misplaced and
the aligners would look far worse than they are — worse toward one end of the genome, which is
a bizarre and hard-to-diagnose signature.

`data/truth/<g>.refseq2simseq.map.txt` gives the reference↔mutated coordinate correspondence
for every injected variant, so the conversion is computable. **Phase 5's
`placement_accuracy.py` must convert mutated → reference coordinates before comparing.**
Recorded here so it cannot be forgotten.

### 3.5 Read length is capped at 150 bp by the platform profile

The HS25 (HiSeq 2500) profile supports read lengths only up to **150 bp**. Going to 250 bp
would require switching to `MSv3` (MiSeq), which is a different instrument with a different
chemistry and a different error profile entirely.

Doing so would **confound read length with platform chemistry**: any change in results between
150 bp and 250 bp could be caused by either, and the design could not separate them. The
read-length sweep is therefore restricted to **75 / 100 / 150**, all on HS25, so read length
varies alone.

### 3.6 Seeds — and which seed controls what

Two distinct random processes, deliberately seeded separately (R9):

- **simuG seed = 20260814**, fixed across the entire project. The mutated genome and truth set
  are the *experimental subject*; they must stay identical across every condition, or a
  coverage comparison would also be comparing two different genomes.
- **ART seed = 1** for the baseline; the sweep will use **1–5**. Each ART seed is a replicate
  *sequencing run* of the same sample, capturing run-to-run variation in which fragments
  happened to be sequenced and where errors happened to land.

The seed appears in every filename (`..._seed1_`) and will appear as a column in
`results/results.tsv`.

### 3.7 No trimming (R8)

Real workflows often trim low-quality read ends before alignment. **We do not**, so all nine
pipelines receive byte-identical input. If reads were trimmed, differences between aligners
could come from how the trimmer interacted with each one rather than from the aligners
themselves. Mean quality is ~Q36.6 (see Phase 4), so there is little to trim anyway.

### 3.8 A tooling trap worth recording: `grep` and very long lines

While verifying coverage interactively I got a nonsensical *E. coli* genome length of 447,535 bp
instead of 4,641,839 — and therefore a coverage of 311× instead of 30×.

Cause: the mutated-genome FASTA that simuG writes is **two lines** — a header plus one
4.6-million-character sequence line. The idiom `grep -v '^>' file | tr -d '\n' | wc -c`
truncated it. The interactive shell on this machine resolves `grep` to a `ugrep` shim, which
does not handle a 4.6 MB line; the scripts, which do not inherit interactive shell functions,
were using `/usr/bin/grep` and had been correct all along.

Both possible outcomes here are bad: a wrong number that looks plausible, or a correct number
that depends on which `grep` is first on `PATH`. So all sequence-length computation now uses
`awk '!/^>/{n+=length($0)} END{print n+0}'`, which sums line lengths natively, has no long-line
limit, and cannot vary with the environment.

**The general lesson:** genomics files routinely contain single lines megabytes long. Line-based
UNIX text tools are not all safe on them, and when they fail they usually truncate silently
rather than erroring.

---

## Phase 4 — Quality control

### 4.1 What FastQC is for here, and what it is *not* for

FastQC reads a FASTQ and reports diagnostics: quality per cycle, GC content, adapter
contamination, duplication. MultiQC merges the four per-file reports into one page at
`results/qc/multiqc_report.html`.

**This phase changes nothing.** It is a sanity check that the simulated reads look like real
sequencing data, so that when a pipeline performs badly later we know it is the pipeline's
fault and not because the input was malformed. Per rule R8 nothing is trimmed — see §4.5.

### 4.2 A Phred score, in plain language

Each base in a FASTQ carries a quality character encoding a **Phred score** `Q`, which is the
sequencer's estimate of the probability that base is wrong:

```
P(error) = 10 ^ (-Q / 10)
```

| Q | P(error) | in words |
|---|---|---|
| 10 | 1 in 10 | terrible |
| 20 | 1 in 100 | poor |
| 30 | 1 in 1,000 | good — the industry "%≥Q30" threshold |
| 40 | 1 in 10,000 | excellent |

This is the entire reason variant calling is hard. At 30× coverage a position is covered by
~30 reads; at Q30 roughly 1 base in 1,000 is wrong, so across a 4.6 Mb genome at 30× there
are ~139 million sequenced bases and therefore **~140,000 erroneous bases**. The truth set
contains only 6,000 real variants. A caller must therefore separate 6,000 real signals from
~140,000 pieces of noise — which is why callers weigh evidence probabilistically rather than
just looking for mismatches.

### 4.3 Results — and the quality profile does look like real Illumina data

| Dataset | mean Q | mean P(error) | Q_effective | % ≥ Q30 | cycle 1 | peak | last 10 cycles | drop from peak |
|---|---|---|---|---|---|---|---|---|
| phiX R1 | 36.60 | 1.57e-03 | 28.03 | 91.8% | 33.1 | 39.3 @ c20 | 35.58 | 3.74 |
| phiX R2 | 36.18 | 2.26e-03 | 26.46 | 89.6% | 33.1 | 38.9 @ c32 | 35.03 | 3.82 |
| *E. coli* R1 | 36.58 | 1.66e-03 | 27.79 | 91.7% | 33.0 | 39.1 @ c16 | 35.52 | 3.63 |
| *E. coli* R2 | 36.18 | 2.23e-03 | 26.51 | 89.7% | 32.8 | 38.6 @ c19 | 35.03 | 3.56 |

**Mean Q ≈ 36.6 (R1) / 36.2 (R2).** All FastQC quality modules PASS.

Three signatures of genuine Illumina data are present:

1. **Low at cycle 1, rising to a peak around cycle 16–32, then declining.** Cycle 1 is ~Q33,
   peaks near Q39, ends near Q35.5. Real Illumina reads behave exactly this way: the first
   cycles are noisy while cluster identification stabilises, and quality then decays as the
   run proceeds because of phasing/pre-phasing (molecules in a cluster gradually fall out of
   sync) and reagent depletion.
2. **The 3′ decline is present: ~3.6–3.8 Phred points from peak to the last 10 cycles.**
3. **R2 is consistently worse than R1** — lower mean Q (36.18 vs 36.58) and a lower 3′ end
   (35.03 vs 35.52). This is real Illumina behaviour: the second read is sequenced after the
   template has spent longer on the flow cell.

Honest caveat: at ~3.7 points, this decline is **milder than many real HiSeq runs**, where the
3′ end can fall to Q30 or below. ART's HS25 profile is empirical but represents one
well-behaved instrument run. So these reads are, if anything, slightly *easier* than typical
real data — worth stating when comparing absolute numbers to a real experiment.

### 4.4 The measurement subtlety: mean Q is not the mean error rate

`scripts/extract_mean_q.py` reports both, and the gap is large:

```
mean Q       = 36.58      -> if taken at face value implies P = 2.2e-04
mean P       = 1.66e-03   -> the actual average error probability
Q_effective  = 27.79      -> mean P expressed back on the Phred scale
```

**A ~9 Phred point gap, i.e. the true error rate is ~7.5× higher than "mean Q 36.6" suggests.**

The cause is that Q is logarithmic, so averaging Q values is averaging exponents. A handful of
very bad bases dominate the true error rate but barely move the arithmetic mean:

> Two bases at Q40 and Q10. Mean Q = 25 (implying P = 0.0032). But the actual mean error
> probability is (0.0001 + 0.1)/2 = 0.05 — Q13. Twelve Phred points apart, a ~16× difference.

**This is why the script exists.** The brief asks to convert ART's `-qs` shift into "a
physically interpretable error-rate feature" for the later modelling phase. `-qs = -5` is an
arbitrary knob; `mean_p` is a physical rate that can be measured on any real dataset and
compared. **The modelling phase should use `mean_p`, not `mean_q`** — a model fitted on mean Q
is fitted on a quantity that systematically understates the noise it is trying to explain.

Cross-validated against FastQC's own per-cycle table (cycle 1: 32.99 vs my 33.0; cycle 111:
35.58 vs my 35.6) so the script is not quietly computing something else.

Output: `results/qc/mean_q.tsv`, one row per FASTQ, with `mean_q`, `mean_p`, `q_effective`,
`frac_q30`, peak cycle, and the 3′ decline.

### 4.5 No trimming, deliberately (R8)

The standard next step in a real workflow would be trimming low-quality 3′ ends with
Trimmomatic or fastp. **We do not trim.**

If we did, every aligner would receive reads shaped by the trimmer's decisions, and an
observed difference between BWA-MEM and Bowtie2 could be caused by how each responds to
variable-length reads rather than by the aligners themselves. Rule R8 requires byte-identical
input to all nine pipelines, and untrimmed reads are the only way to guarantee it.

There is little to trim in any case: 91.7% of bases are ≥ Q30 and the worst cycles average
Q35.

### 4.6 The two FastQC warnings are both expected

| Warning | Where | Explanation |
|---|---|---|
| Per sequence GC content | *E. coli*, phiX R2 | FastQC compares the observed GC distribution against a theoretical normal curve fitted to the data. A single-organism sample has a narrow, sharply-peaked GC distribution (here 50%, matching *E. coli*'s 50.79%) that deviates from that model. The warning fires on essentially every clean single-genome dataset. |
| Overrepresented sequences | phiX only | phiX is 5,386 bp covered at 30× by 150 bp reads. Reads necessarily overlap heavily and identical reads recur by chance. On a 5 kb genome this is arithmetic, not contamination. |

Neither is a data problem, and neither warrants action.

### 4.7 A number to carry into Phase 5

FastQC reports **97.64% "Total Deduplicated Percentage"** for *E. coli* — i.e. ~2.4% of reads
share a sequence with another read.

These are **not PCR duplicates**. They are coincidental collisions: at 30× coverage of a
4.6 Mb genome, two independently sampled fragments occasionally start at the same position.
No PCR amplification was simulated, so no true duplicates exist.

This sets the expectation for Phase 5: `gatk MarkDuplicates` should mark ~0%. It is run
anyway, for pipeline realism, and its inertness here must be stated in the report rather than
presented as a meaningful result.

---

## Phase 5 — Alignment

### 5.1 What an aligner does

Each read is a 150-letter string that came from somewhere in a 4.6-million-letter genome. The
aligner finds where. It is hard for two reasons: the read contains sequencing errors and real
variants (so it will not match exactly anywhere), and the genome contains repeated sequence
(so a short read may match several places equally well).

The output is a **BAM** file: one record per read giving its position, orientation, a CIGAR
string describing how it lines up (matches, insertions, deletions, clipping), and a **MAPQ**
score — the aligner's own confidence, on a Phred scale, that this is the right location.

**R1 again:** reads are aligned to `data/refs/<g>.fa`, the *original* reference, never to the
mutated genome. `align_reads.sh` asserts this before doing anything.

### 5.2 Read groups — three tools, three different spellings (R6)

A read group tags reads with which sample and sequencing run they came from. GATK refuses to
run without one. All three aligners were given the identical group
`ID:s1 SM:sim PL:ILLUMINA LB:lib1`, but they will not accept it the same way:

| Aligner | Syntax |
|---|---|
| BWA-MEM | `-R '@RG\tID:s1\tSM:sim\tPL:ILLUMINA\tLB:lib1'` — one tab-delimited string |
| minimap2 | same `-R` string as BWA |
| **Bowtie2** | `--rg-id s1 --rg SM:sim --rg PL:ILLUMINA --rg LB:lib1` — ID separately, then one `--rg` per field |

Passing BWA's tab-delimited string to Bowtie2 produces a malformed header that GATK rejects
later, in a completely different phase, with an error that does not mention read groups.
Verified after alignment that all three BAMs carry byte-identical `@RG` lines.

### 5.3 A launcher trap: `gatk` is a Python script, not a binary

`MarkDuplicates` failed with:

```
env: python: No such file or directory
```

`gatk` is a **Python wrapper** that locates and launches the GATK jar. Calling it by absolute
path (`$CONDA/envs/callers/bin/gatk`) is not enough — the wrapper itself runs `env python`, so
its environment's `bin/` must be on `PATH`. Everywhere else in this project, invoking tools by
absolute path is the right call (it avoids `conda run` overhead and its broken-pipe behaviour,
see §2 tooling notes), but **gatk is the exception** and needs `PATH="$CONDA/envs/callers/bin:$PATH"`.
This applies to Phase 6's HaplotypeCaller too.

### 5.4 Timing methodology — the aligner alone

The brief's example pipes each aligner straight into `samtools sort`. That is the normal
production idiom, but it is a poor measurement: sorting costs roughly the same for all three
aligners (identical read counts) and would dilute the very difference being measured.

So each aligner is timed **alone**, writing SAM to disk; sorting, duplicate marking and
indexing happen afterwards, untimed. Applied identically to all three, so the comparison stays
fair (R8), and all three get 4 threads.

`/usr/bin/time` on macOS is BSD, not GNU: the flag is **`-l`**, not `-v`, and **peak RSS is
reported in bytes**, whereas GNU `time -v` reports kilobytes. Getting that wrong misreports
memory by 1024×.

### 5.5 Fairness fix: count primary reads, not "in total" (R8)

`samtools flagstat`'s "in total" line includes **supplementary** alignments — extra records an
aligner emits when it splits a chimeric read across two locations. BWA-MEM emits them (36 on
*E. coli*); Bowtie2 and minimap2 `-ax sr` do not.

Using "in total" gave BWA a denominator of 928,386 against 928,350 for the other two, for the
*same 928,350 input reads* — so the three "mapping rates" were not fractions of the same
quantity. The metric now uses `primary` / `primary mapped`, which is exactly one record per
input read for every aligner. The supplementary count is reported as its own column rather
than being folded into a rate.

This is small (0.004%) but it is precisely the class of error that makes a benchmark
indefensible: not wrong enough to notice, entirely wrong in principle.

### 5.6 Placement accuracy — and the coordinate conversion, quantified

`scripts/placement_accuracy.py` compares each aligner's BAM against ART's truth SAM and reports
the fraction of reads placed within ±10 bp of where they really came from. No caller is
involved, so this is a pure aligner metric.

As predicted in §3.4, the two files are in different coordinate systems (truth SAM = mutated
genome, BAM = original reference). The script builds an explicit mutated→reference map from
simuG's injected indels, treating three cases separately: unchanged stretches (a constant
shift), deleted reference bases (no mutated equivalent), and **inserted bases, which have no
reference position at all** and are clamped to the anchor base rather than shifted (otherwise
bases inside a 50 bp insertion land up to 50 bp away — larger than the tolerance, producing
fake errors concentrated at insertion sites).

The map **validates itself**: simuG records each indel's own mutated coordinate as `sim_start`,
and the independently-derived map must reproduce it. All 10 phiX and all 1,000 *E. coli* indels
agree. If they did not, the script hard-errors rather than reporting numbers.

**What the conversion was worth**, BWA on *E. coli*:

| | Placement accuracy | Median \|offset\| |
|---|---|---|
| With conversion (correct) | **99.028%** | 0 bp |
| Without conversion (naive) | **8.903%** | 38 bp |

Drift along the genome, and note it is not monotonic — it wanders as insertions and deletions
alternate, ending at the net +187:

```
mutated 500,000 -> reference   499,968   (drift  +32 bp)
mutated 2,000,000 -> reference 2,000,053   (drift  -53 bp)
mutated 4,641,000 -> reference 4,640,813   (drift +187 bp)
```

A naive comparison would have reported all three aligners at ~9% placement accuracy — a
catastrophic-looking result, with no error message, that is entirely an artefact of the
measurement.

### 5.7 Results

| Genome | Aligner | Mapping | Properly paired | Mean MAPQ | MAPQ0 | Mean depth | Placement | Runtime | Peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| phiX | BWA-MEM | 100% | 100% | 60.00 | 0% | 29.89 | 99.815% | 0.01 s | 3.1 MB |
| phiX | Bowtie2 | 100% | 100% | 40.90 | 0% | 29.93 | 100.000% | 0.09 s | 57.5 MB |
| phiX | minimap2 | 100% | 100% | 59.98 | 0% | 29.89 | 99.815% | 0.00 s | 3.9 MB |
| *E. coli* | BWA-MEM | 100% | 100% | 59.06 | 1.274% | 29.98 | 99.028% | 5.16 s | 382 MB |
| *E. coli* | Bowtie2 | 99.981% | 99.833% | 41.17 | 0.020% | 29.98 | 99.014% | 17.55 s | 67.1 MB |
| *E. coli* | minimap2 | 100% | 100% | 59.08 | 1.275% | 29.98 | 99.019% | 1.85 s | 490.5 MB |

Observations worth defending in a viva:

- **Placement accuracy is essentially identical across all three (99.01–99.03% on *E. coli*).**
  The remaining ~1% is not aligner weakness; it is reads drawn from repeated sequence where the
  correct location is not recoverable from a 150 bp read. All three hit the same information
  limit. **Do not expect the aligner to be the discriminating factor at baseline** — with 30×
  coverage, 150 bp reads and a small bacterial genome, this is an easy alignment problem.
  Differences should emerge at lower coverage and shorter reads, which is what the sweep is for.
- **Mean MAPQ differs a lot (BWA/minimap2 ≈ 59, Bowtie2 ≈ 41) but this is a scale difference,
  not a quality difference.** MAPQ is defined per-tool: BWA caps at 60, Bowtie2 at 42. Comparing
  the raw numbers across tools is meaningless. It matters anyway, because **callers filter and
  weight on MAPQ using tool-agnostic thresholds** — so an identical threshold is a stricter
  filter on Bowtie2 output than on BWA output. This is a real confound for Phase 6 and one of
  the more interesting things this benchmark can quantify.
- **MAPQ0 (multi-mapping) reads: BWA and minimap2 flag ~1.27%, Bowtie2 flags 0.02%.** They are
  looking at the same repeats; they differ in how they report ambiguity. Bowtie2 by default
  reports one alignment for a multi-mapping read with a low-but-nonzero MAPQ, rather than
  marking it 0. Since most callers discard MAPQ0 reads outright, BWA and minimap2 effectively
  hand the caller ~1.25% less usable coverage in repeats.
- **Runtime: minimap2 (1.85 s) < BWA (5.16 s) < Bowtie2 (17.55 s)** — minimap2 ~9.5× faster
  than Bowtie2. But memory is inverted: Bowtie2 uses 67 MB where minimap2 uses 490 MB, a ~7×
  difference. That is a genuine engineering trade-off, not a defect.
- phiX runtimes (0.00–0.09 s) are **too short to be meaningful** and should not be reported as
  a speed comparison. They are below timer resolution and dominated by process startup.

### 5.8 MarkDuplicates is inert here, as predicted

Duplicate rate: **0.00% on phiX, ~0.027% on *E. coli*** for all three aligners.

PCR duplicates arise when library amplification copies the same original fragment many times,
producing reads that look like independent evidence but are not. ART simulates no PCR, so
there are none. The ~0.027% on *E. coli* is coincidental collision — at 30× coverage two
independently drawn fragments occasionally share a start position — which matches FastQC's
independent 97.64%-deduplicated figure from §4.7.

It is run for pipeline realism (a real workflow has this step, and omitting it would make the
pipeline unrepresentative), but **its inertness must be stated in the report** rather than
presented as a finding. It also means MarkDuplicates cannot be a source of difference between
the nine pipelines here.

---

## Phase 6 — Variant calling

### 6.1 What a variant caller does

The aligner produced a pile of reads stacked over every genome position. The caller looks at
each position and asks: *do the reads here disagree with the reference in a way that is better
explained by a real difference than by sequencing error?*

At 30× coverage with ~1-in-1,000 base error, a position covered by 30 reads will show a
spurious mismatch reasonably often. A single mismatching read is almost certainly error; 28 of
30 reads agreeing on a different base is almost certainly real. The callers differ in how they
do this reasoning:

- **GATK HaplotypeCaller** — locally *reassembles* reads into candidate haplotypes and
  realigns against them. Slow, but robust around indels because it does not trust the aligner's
  original per-read alignment in messy regions.
- **FreeBayes** — Bayesian, haplotype-based over short windows.
- **BCFtools mpileup/call** — the classic pileup model, position by position. Fastest.

### 6.2 R2 — the ploidy verification, and why it is not optional

**All 18 runs (2 genomes × 3 aligners × 3 callers) PASS: every genotype is haploid.**
Logged to `logs/ploidy_verification.txt`.

The flags, spelled differently by each tool:

| Caller | Flag |
|---|---|
| GATK | `--sample-ploidy 1` |
| FreeBayes | `-p 1` |
| BCFtools | `bcftools call --ploidy 1` |

`bcftools`' flag needed checking: `--ploidy` normally takes a *predefined assembly name*
(e.g. `GRCh37`) or a ploidy file, so `--ploidy 1` looks like it might be silently misparsed.
Verified empirically before the real run — it is accepted and yields `GT=1` with `AN=1`.

The verification greps the `GT` field directly and counts any genotype containing `/` or `|`
(the diploid separators). The script **aborts** on a single diploid genotype rather than
continuing.

**Why this cannot be replaced by looking at the F1 score:** Phase 2 established that a caller
emitting `1/1` everywhere scores a *perfect* F1 against haploid truth, because vcfeval treats
homozygous-ALT as equivalent to haploid-ALT. A ploidy misconfiguration of that kind is
completely invisible in the results table. Only the direct GT check catches it.

### 6.3 R7 — no BQSR, and why that is the fair choice

Base Quality Score Recalibration is part of GATK Best Practices. It learns systematic biases in
the sequencer's quality scores by assuming that mismatches at *known* variant sites are real
and everything else is error.

That requires a database of known variants. **None exists for phiX or *E. coli*.** It could be
bootstrapped — call variants, treat confident calls as "known", recalibrate, re-call — but that
would give GATK an extra data-driven preprocessing step that FreeBayes and BCFtools do not get,
and any GATK advantage afterwards could not be attributed to the caller. Skipped deliberately
(R7). This should be stated in the report: **the GATK arm here is deliberately not the full
Best Practices pipeline**, and that is a fairness decision, not an oversight.

### 6.4 A launcher trap, part two

Phase 5 established that `gatk` needs its environment's `bin/` on `PATH`. Wrapping it in a
shell function was not enough here, because `/usr/bin/time` **execs a real binary and cannot
run a shell function**:

```
time: gatk: No such file or directory
```

which reads like a missing installation rather than a quoting problem. Fixed by invoking
`env PATH=... /path/to/gatk`.

### 6.5 Raw counts and runtimes

Truth: phiX 60 variants, *E. coli* 6,000.

| Genome | Aligner | GATK | FreeBayes (GT=1) | BCFtools |
|---|---|---|---|---|
| phiX | bwa | 60 | 57 (+4 GT=0) | 60 |
| phiX | bowtie2 | 60 | 57 (+7 GT=0) | 60 |
| phiX | minimap2 | 60 | 57 (+4 GT=0) | 60 |
| *E. coli* | bwa | 5,946 | 5,919 (+2,339 GT=0) | 5,948 |
| *E. coli* | bowtie2 | 5,896 | 5,939 (+2,858 GT=0) | 5,946 |
| *E. coli* | minimap2 | 5,938 | 5,920 (+2,340 GT=0) | 5,944 |

Runtime, *E. coli* (seconds, 4 threads where supported):

| Caller | bwa | bowtie2 | minimap2 |
|---|---|---|---|
| GATK | 31.25 | 31.00 | 31.72 |
| FreeBayes | 10.61 | 10.84 | 10.66 |
| BCFtools | 9.92 | 9.72 | 9.82 |

GATK is ~3× slower than the other two, which is the expected cost of local reassembly. Caller
runtime is essentially independent of which aligner produced the BAM.

**FreeBayes' `GT=0` records.** FreeBayes emits candidate sites it evaluated and *rejected* —
they carry an ALT allele but a genotype of `0` (reference) and `QUAL=0`. They are not variant
calls. vcfeval treats `GT=0` as non-variant, and the shared QUAL≥20 filter removes them anyway,
so they do not inflate FreeBayes' false positives. But they do mean **raw record counts are not
comparable across callers** — FreeBayes' 8,258 records represent 5,919 actual calls.

### 6.6 The important discovery: FreeBayes writes complex variants, and it breaks type-splitting

FreeBayes appeared to *miss* 3 of 60 phiX variants. It does not. It **represents them
differently**:

| Truth (atomic) | GATK | FreeBayes |
|---|---|---|
| `215 A>T` and `217 A>T` | two SNP records | **one** record `215 AAA>TAT` |
| `1234 T>C` and `1235 G>A` | two SNP records | **one** record `1234 TG>CA` |
| `1699 C>G` and `1702 TC>T` | SNP + indel | **one** record `1699 CCGTCCTT>GCGTCTT` |

FreeBayes is haplotype-based and merges nearby variants into a single MNP or complex record.
Both representations describe exactly the same sequence. This is precisely the situation rule
R4 exists for.

**Why this matters far more than it first appears.** The brief specifies scoring SNVs and
indels *separately*, via `bcftools view -v snps` / `-v indels`. But bcftools classifies a
record by its overall shape, and a complex record is neither a SNP nor an indel:

```
FreeBayes phiX raw:  -v snps 49   -v indels 9   -v mnps 2   -v other 1
truth:                  50 snps      10 indels
```

Splitting by type *before* vcfeval would silently discard the 2 MNP and 1 complex record —
along with the 6 real variants inside them. FreeBayes' SNP recall would read 49/50 and its
indel recall 9/10, and the deficit would look like a genuine sensitivity difference. It is
purely a representation artefact.

**Fix: add `--atomize` to the normalisation step (R3), applied identically to truth and to all
nine call sets.** `bcftools norm -f ref -m -any --atomize` decomposes MNVs and complex records
into consecutive atomic SNVs and indels. Verified on phiX:

```
FreeBayes after --atomize:  -v snps 54   -v indels 10   -v mnps 0   -v other 0
```

and the six previously-hidden variants reappear at exactly the truth positions and alleles
(215 A>T, 217 A>T, 1234 T>C, 1235 G>A, 1699 C>G, 1702 TC>T). Indels now match truth exactly at
10. (The 54 SNPs include the 4 rejected `GT=0` candidate sites, which are not calls.)

This changes the Phase 7a normalisation command from the brief's

```
bcftools norm -f <ref.fa> -m -any
```

to

```
bcftools norm -f <ref.fa> -m -any --atomize
```

applied to **the truth set and every call set with identical settings**, as R3 demands. The
truth set contains no MNPs so atomising it changes nothing — but it must still be atomised,
because R3's requirement is identical *treatment*, not identical *outcome*.

**Generalisable lesson:** rule R4 says never hand-roll variant comparison because the same
variant has many valid representations. This shows the same hazard reaching *upstream* of the
comparison: even when using vcfeval correctly, a `bcftools view -v snps` in the preparation
step can reintroduce exactly the representation-sensitivity that vcfeval was chosen to avoid.

### 6.7 Hard filtering — identical logic, with a stated caveat (R8)

Both raw and filtered call sets are produced. The filter is identical for every caller:

```
QUAL >= 20 && INFO/DP >= 5
```

Only `QUAL` and `DP` are used, because they are the only fields all three callers emit with
comparable meaning. GATK Best Practices would filter on `QD`, `FS`, `MQRankSum` and similar —
but FreeBayes and BCFtools do not produce those annotations, so using them would apply a
better-tuned filter to GATK than to its competitors and confound the comparison.

**The caveat, stated plainly: `QUAL` is not calibrated identically across these three tools.**
A QUAL of 20 does not mean the same thing to GATK as to FreeBayes. So a single threshold is
*procedurally* identical but not *statistically* equivalent — it is the fairest available
choice, not a perfect one.

This is exactly why the raw call sets retain QUAL and why **ROC curves are the honest
comparison**: the ROC sweeps the threshold across its whole range, removing the arbitrariness
of any single cut-off. The hard-filtered numbers should be read as one operating point on that
curve, not as the result.

Filtering effect on *E. coli* is small for GATK and BCFtools (~0.1–0.4% removed) but large for
FreeBayes (8,258 → 5,914), because it removes the rejected `GT=0` candidate records.

---

## Phase 7a — Normalisation and GA4GH scoring

### 7.1 The metrics, in plain language

vcfeval sorts every variant into one of three bins:

- **TP (true positive)** — in the truth set *and* called. Correct.
- **FP (false positive)** — called but not in truth. A variant invented from noise.
- **FN (false negative)** — in truth but not called. A variant missed.

From these:

```
precision = TP / (TP + FP)   "of the variants I reported, what fraction were real?"
recall    = TP / (TP + FN)   "of the variants that exist, what fraction did I find?"
F1        = 2 * precision * recall / (precision + recall)
```

**Why F1 is the headline, with the worked example the brief asks for.** Precision and recall
trade off, and either alone is trivially gameable. Consider a caller on *E. coli* that reports
only its single most confident variant and nothing else:

```
TP = 1, FP = 0, FN = 5,999
precision = 1 / 1        = 1.0000     <- perfect!
recall    = 1 / 6,000    = 0.000167
F1        = 2(1)(0.000167)/(1.000167) = 0.000333
```

Precision 1.0 looks flawless while the caller found 0.017% of the variants. F1 is the harmonic
mean, which is dominated by the *smaller* of the two, so it collapses to ~0.0003 and correctly
calls this useless. That is why F1 is reported as the headline.

But F1 is still **one operating point**, determined by whatever QUAL threshold was applied —
which is why the full ROC curve is also produced (§7.7).

### 7.2 Normalisation — identical treatment, verified (R3)

Every call set is normalised with the command applied to the truth set, character for character:

```
bcftools norm -f <ref.fa> -m -any --atomize
```

`--atomize` was added after the Phase 6 discovery (§6.6). The truth set was **re-normalised**
with it too — it contains no MNVs so nothing changed, but R3 requires identical *treatment*,
not identical *outcome*.

Comparison is done **only** by `rtg vcfeval` (R4). No position-and-allele string matching
anywhere in this project.

### 7.3 phiX: all nine pipelines score F1 = 1.0000 — and why that is real

Every pipeline: **50/50 SNVs, 10/10 indels, zero FP, zero FN.**

The brief warns that an F1 of exactly 1.0 usually indicates a bug. That warning is correct and
was taken seriously, so `scripts/scoring_negative_control.sh` deliberately breaks the call set
three ways and confirms the machinery punishes each:

| Call set | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| unmodified | 60 | 0 | 0 | 1.0000 | 1.0000 | **1.0000** |
| all positions shifted +5 bp | 0 | 60 | 60 | 0.0000 | 0.0000 | **0.0000** |
| every SNV ALT changed | 10 | 50 | 50 | 0.1667 | 0.1667 | **0.1667** |
| half the calls removed | 30 | 0 | 30 | 1.0000 | 0.5000 | **0.6667** |

The third row is a good internal check: only the 10 indels were left untouched, and exactly 10
TPs survive. The fourth behaves exactly as arithmetic demands. **The scoring genuinely
discriminates**, so the perfect phiX scores are a real result, not a broken comparison.

They are also unsurprising. phiX is 5,386 bp with no repetitive sequence; at 30× with 150 bp
reads every variant is unambiguously recoverable, and Phase 5 already showed ~100% placement
accuracy. **phiX cannot discriminate between pipelines and should not be presented as if it
does.** Its role is exactly what R10 says: prove the pipeline works end to end, fast.

(A negative-control note worth keeping: the first version of the ALT-mutation control rotated
A→C→G→T→A, which sometimes produced `REF == ALT`. vcfeval correctly *refuses* such a record
rather than scoring it. The control now picks a base different from both REF and the original
ALT.)

### 7.4 *E. coli* baseline results (PRIMARY — single-run scoring)

**SNV F1**

| | GATK | FreeBayes | BCFtools |
|---|---|---|---|
| **BWA-MEM** | **0.9953** | **0.9953** | 0.9948 |
| **Bowtie2** | 0.9909 | 0.9906 | 0.9914 |
| **minimap2** | 0.9946 | **0.9953** | 0.9944 |

**Indel F1**

| | GATK | FreeBayes | BCFtools |
|---|---|---|---|
| **BWA-MEM** | **0.9975** | 0.9970 | 0.9970 |
| **Bowtie2** | 0.9935 | 0.9890 | 0.9815 |
| **minimap2** | 0.9970 | **0.9975** | 0.9970 |

What the numbers actually say:

- **The aligner matters more than the caller at this baseline.** Every Bowtie2 row is worse
  than the corresponding BWA or minimap2 row, for both variant types. The spread across
  aligners (SNV 0.9906–0.9953) is wider than across callers within an aligner. Note this is
  *despite* Phase 5 showing all three aligners at ~99.0% placement accuracy — so the
  difference is not mainly about where reads land.
- **The likely mechanism is MAPQ, not placement.** Bowtie2's MAPQ scale tops out at 42 versus
  60 for the others (§5.7). Callers apply tool-agnostic MAPQ thresholds, so an identical
  internal cut-off is a *stricter* filter on Bowtie2 output. Bowtie2+GATK has the highest SNV
  FN count of any pipeline (90 vs 47 for BWA+GATK) with **zero** false positives — the
  signature of a caller discarding evidence, not of an aligner misplacing reads. This is a
  hypothesis consistent with the data, not something this experiment has yet proven; the
  coverage sweep should test it.
- **Precision is near-perfect almost everywhere; recall is what separates pipelines.** Six of
  nine SNV pipelines have precision exactly 1.0000. The variation is essentially all in recall
  — these callers are conservative, and at 30× on a small genome false positives are rare.
- **Indels are harder than SNVs but not dramatically so at 30×.** The worst indel pipeline
  (Bowtie2+BCFtools, 0.9815) is meaningfully worse than the worst SNV pipeline (0.9906), and
  the indel spread (0.9815–0.9975) is wider than the SNV spread. Expect this gap to widen as
  coverage drops.
- **Nothing is 0.0 or 1.0 on *E. coli*.** Every pipeline misses 5–18 indels and 44–90 SNVs out
  of 1,000 and 5,000. That is a plausible result, not a bug signature.

### 7.5 Anomalies flagged

1. **phiX F1 = 1.0000 everywhere.** Investigated with the negative control above. Genuine, but
   phiX has no discriminating power and must not be reported as a comparison.
2. **Bowtie2+FreeBayes has the *highest* SNV recall of any pipeline (0.9912) yet a
   middling F1**, because it also has by far the most false positives (50, versus 0 for six
   other pipelines). It is the one genuinely aggressive pipeline in the set — a real
   precision/recall trade-off, visible properly only on the ROC.
3. **Bowtie2+BCFtools indels: 19 FP and 18 FN**, the worst on both axes simultaneously. Worth
   watching in the sweep; it may be a Bowtie2 MAPQ/indel-alignment interaction.

### 7.6 IMPORTANT — scoring method: pre-splitting by type distorts the result

The brief specifies scoring SNVs and indels separately by splitting the VCF first
(`bcftools view -v snps` / `-v indels`, then vcfeval on each). Both that method and RTG's
native alternative were run, and **they disagree systematically**.

`rtg vcfeval` matches variants **haplotype-aware**: it reconstructs the local haplotype implied
by a set of variants and compares *sequence*, which is exactly what makes it immune to
representation differences (R4). Pre-splitting breaks that. Remove the indels from a call set
and a SNV sitting next to an indel can no longer be reconciled with truth, so vcfeval charges
it as **both** a false negative and a false positive.

Measured on *E. coli*, BWA + FreeBayes, indels:

| Method | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| Pre-split (brief's method) | 983 | 9 | 17 | 0.9909 | 0.9830 | 0.9869 |
| Single run, RTG's own per-type split | **994** | **0** | **6** | 1.0000 | 0.9940 | **0.9970** |

11 true positives destroyed and 9 false positives invented, purely by splitting the file.
**70.6% of the pre-split method's extra false negatives have another truth variant within
50 bp** — precisely the neighbours whose haplotype context was removed.

The effect is not uniform, so it changes conclusions rather than just shifting all numbers:
it penalises FreeBayes hardest (haplotype-based calling depends most on context), making
FreeBayes' indel performance look distinctly worse than GATK's when in fact they are
comparable (0.9970 vs 0.9975).

**Resolution.** A single vcfeval run per pipeline on the full call set, taking the per-type
breakdown from RTG's own `snp_roc.tsv.gz` / `non_snp_roc.tsv.gz`. Validated: per-type TP/FP/FN
sum **exactly** to the combined `summary.txt` for every pipeline checked.

`results/results.tsv` carries a **`scoring_method`** column with both — `single_run` (primary)
and `pre_split` (the brief's method, retained for comparison) — so the difference is auditable
rather than hidden.

**The generalisable lesson**, and the deepest one in this project: rule R4 says do not
hand-roll variant comparison, because one variant has many valid representations. Phase 6
showed that hazard reaching upstream into *preparation* (type-splitting hides variants inside
MNVs). This shows it again at a subtler level — even with correct atomisation, splitting the
call set before scoring silently degrades a haplotype-aware comparator into a
context-free one. **`bcftools view -v snps` is not a neutral operation before vcfeval.**

### 7.7 ROC curves

`results/roc_{ecoli,phiX}_{SNV,INDEL}.svg`, 9 curves each, generated from the single-run ROC
files with `rtg rocplot`.

The ROC is the honest comparison and the F1 table is a summary of it. A single F1 depends on
one QUAL threshold, and **QUAL is not calibrated identically across GATK, FreeBayes and
BCFtools** (§6.7) — so comparing single F1 values partly compares the callers' QUAL scales
rather than their ability to find variants. Sweeping the threshold removes that. The
Bowtie2+FreeBayes case in §7.5 is exactly the situation only a ROC resolves: highest recall,
most false positives, unremarkable F1.

(`rtg rocplot` refuses to overwrite an existing SVG, so `make_roc.sh` removes the target first;
otherwise every re-run fails.)

### 7.8 `results/results.tsv`

Written by `scripts/collect_results.py`, which joins four sources: vcfeval output,
`results/align_metrics.tsv`, `results/placement_accuracy.tsv`, and `logs/call_timing.tsv`.

**72 rows** = 2 genomes × 9 pipelines × 2 variant types × 2 scoring methods. The
`scoring_method == single_run` subset is exactly the **36 baseline rows** the definition of
done requires.

The condition axes (genome, coverage, read_length, qs_shift, seed) are parsed out of the run
tag, which is why tags encode them as `<genome>_cov30_len150_err0_seed1` — the sweep in Phase 8
adds rows without any schema change.

---

## Phase 7b — Snakemake workflow

### 7b.1 What Snakemake adds, and why it is not just a shell script with extra steps

The Phase 3–7a scripts run a phase at a time, in order, by hand. Snakemake instead describes
each step as a **rule** with declared inputs and outputs, then works out the dependency graph
itself. Three things follow that the scripts cannot give:

1. **Resumability.** A failure in run 700 of the sweep does not cost the first 699.
   Snakemake re-runs only what is missing or out of date.
2. **Parallelism.** With `--cores 8` it schedules independent jobs concurrently without any
   manual job management.
3. **Provenance.** The graph *is* the documentation: every output states exactly which inputs
   and which command produced it.

### 7b.2 Fine-grained rules, not script wrappers

The scripts each do a whole phase in one invocation (all three aligners; all nine caller runs).
Wrapping them would have produced a five-node DAG with no real parallelism and no
per-pipeline resumability — one failed caller would force re-running all nine.

So the Snakefile re-expresses the pipeline as **one rule per tool invocation**. The cost is
duplicated command lines: the Snakefile and the scripts must stay in step. **GATE 7b guards
exactly this** — the baseline is rebuilt from scratch through Snakemake and the numbers must
match Phase 7a. That check is meaningful precisely *because* the two are independent
implementations (see §7b.5).

### 7b.3 The sweep design table

`config/conditions.tsv`, generated by `scripts/make_conditions.py`:

| Axis | Levels | Conditions |
|---|---|---|
| Coverage | 5, 10, 20, **30\***, 50, 100 | 6 |
| Read length | 75, 100, **150\*** | 3 |
| qs shift | **0\***, −2, −5, −10 | 4 |

13 minus the baseline counted three times = **11 unique conditions per genome**.
× 5 seeds × 2 genomes = **110 rows**; × 9 pipelines = **990 pipeline runs**.

**This is one-factor-at-a-time (OFAT), and its limitation should be stated in the report.**
A full factorial would be 6 × 3 × 4 = 72 conditions per genome (6,480 pipeline runs), far
beyond budget. OFAT answers "what does each factor do on its own" for 11 conditions — but it
**cannot detect interactions**. If a pipeline only degrades when coverage is low *and* reads
are short, this design will never see it, because it never moves two axes at once. A sensible
Phase 8+ refinement is to add a small factorial patch around whatever region the OFAT sweep
flags as interesting.

### 7b.4 Making per-rule conda environments real, not decorative

This was the substantive work of the phase, and it exposed two problems.

**Problem 1: the scripts ignored PATH.** They call tools by absolute path into named conda
environments — deliberate, because `conda run` costs ~1s per call and breaks pipes (§2). But
Snakemake's `conda:` directive works by *activating* an environment and putting its tools on
PATH. A script that ignores PATH would silently keep using the developer's local
environments, and the workflow would *appear* portable while not being portable.

Fixed with `scripts/lib/tools.sh`: `resolve_tool` prefers whatever is already on PATH and
falls back to the named environment. The same scripts are now correct under
`snakemake --use-conda`, under plain `snakemake`, and when run directly.

**Problem 2: the exported environment YAMLs could not be re-solved.** `envs/*.yaml` was a full
`conda env export` — every transitive dependency pinned to an exact version. That is a fine
*record* but a poor *specification*: it over-constrains packages we never chose, and
re-solving it on another machine frequently fails.

`scripts/export_envs.sh` now writes two files per environment, with different jobs:

| File | Contents | Purpose |
|---|---|---|
| `envs/<n>.yaml` | only the packages we asked for, pinned (e.g. `bwa=0.7.19`) | what Snakemake installs; portable |
| `envs/<n>.lock.yaml` | full `--explicit` URL lock of this exact build set | proves what *these* numbers came from |

A separate minimal `envs/python.yaml` serves the collector rule, which uses only the Python
standard library — pointing it at the full `ml` environment would install jupyterlab,
scikit-learn and matplotlib to run a script that needs none of them.

**Two traps worth recording:**

- Snakemake must be able to find `conda` itself. Launching it from an environment where
  `conda` is not on PATH fails with `Error running conda info`, which reads like a broken
  conda install rather than a PATH problem.
- Running the workflow **without** `--use-conda` fails every rule with **exit 127**
  (command not found), because the rules use bare tool names. That is the correct behaviour —
  it is the workflow refusing to silently use whatever happens to be installed — but the exit
  code alone does not say so.

### 7b.5 GATE 7b result — the baseline reproduces exactly

The Phase 7a `results.tsv` was moved aside and the baseline rebuilt from scratch:

```
snakemake --use-conda --conda-frontend conda --forceall --cores 8
252 of 252 steps (100%) done      0 errors
```

This re-simulated the reads, re-ran all three aligners, re-marked duplicates, re-ran all nine
caller combinations per genome, re-normalised, and re-scored — inside **freshly created conda
environments**, using command lines written independently of the shell scripts.

**Result: `results/results.tsv` is byte-for-byte identical to Phase 7a. All 72 rows match on
TP, FP, FN, precision, recall and F1.**

The workflow's own R2 gate also fired: 18 `logs/ploidy/*.txt` files, **18 PASS, zero diploid
genotypes**. The `verify_ploidy` rule is a hard dependency of `normalise_calls`, so scoring
cannot proceed past a caller that emits diploid genotypes.

A second dry run reports "Nothing to be done", confirming the workflow is idempotent rather
than rebuilding on every invocation.

**Honest scope note:** the vcfeval-derived columns (TP/FP/FN/precision/recall/F1) were fully
regenerated by Snakemake. The joined metric columns (`mapping_rate`, `mean_mapq`,
`mean_depth`, `placement_accuracy`, `align_seconds`, `call_seconds`, `peak_rss_mb`) come from
`results/align_metrics.tsv`, `results/placement_accuracy.tsv` and `logs/call_timing.tsv`,
which are still produced by the Phase 5/6 shell scripts and are **not yet Snakemake rules**.
They were therefore carried over, not recomputed. Adding those as rules is the first item in
`HANDOFF.md`.

### 7b.6 Dry-run job counts

| Target | Jobs |
|---|---|
| Baseline, forced (2 genomes × 9 pipelines, seed 1) | **252** |
| **Full sweep, forced (11 conditions × 5 seeds × 2 genomes × 9 pipelines)** | **13,320** |
| Full sweep from current state (baseline already built) | 13,070 |

The full-sweep DAG **resolves cleanly** — no ambiguity, no missing inputs, no cycles. Internal
consistency of the graph against the design:

```
simulate_reads   110  = 11 conditions x 5 seeds x 2 genomes
call_gatk        330  = 110 x 3 aligners          (same for freebayes, bcftools)
verify_ploidy    990  = 110 x 9 pipelines         <- the "990-run sweep"
vcfeval_all      990  = one primary scoring run per pipeline
```

`verify_ploidy = 990` is the direct confirmation that the graph encodes the intended design.

**Only the baseline has been executed.** The sweep needs sign-off (out of scope, §5 of the brief).

### 7b.7 DAG exports

- `results/workflow_dag.svg` — the full baseline job graph. Complete but dense.
- `results/workflow_rulegraph.svg` — the rule-level graph. **This is the one to put
  in a report**: it shows the pipeline's structure rather than every individual job.

(Both were regenerated at the end of Phase 8, after the metrics, analysis and report rules were
added: the rule graph now has 34 rules, ending `collect_results → analyse_sweep / fit_model /
diagnose_errors → make_figures → build_report → render_report`.)

(`rtg rocplot` and `rtg vcfeval` both refuse to write over existing output, so the relevant
rules `rm -rf` their target first; otherwise every re-run fails on the second invocation.)

### 7b.8 One rule-design trap

`index_vcf` (`work/{prefix}.vcf.gz` → `.tbi`) and a separate `index_norm_vcf`
(`work/norm/{prefix}.vcf.gz` → `.tbi`) looked like a general rule plus a more specific one.
They are not: **Snakemake wildcards match `/` by default**, so `{prefix}` in the general rule
already covers `norm/...` and the two rules are *ambiguous*, not nested. Snakemake refuses to
build the DAG with `AmbiguousRuleException`. Resolved by deleting the redundant rule.

---

## Phase 8 — Completing the benchmark: truth-set realism, workflow hardening, the sweep

Phase 8 executes the plan in `docs/PHASE1_COMPLETION_PLAN.md`: fix what blocked the sweep,
run all 990 pipeline runs, and analyse them. Several things went wrong on the way, and each is
recorded here because each is a way a benchmark can quietly produce wrong numbers.

### 8.1 Truth sets regenerated with a realistic Ti/Tv — and a free controlled experiment

The original truth sets used simuG's default `-titv_ratio 0.5`, which is what *uniformly random*
substitution produces (each base has one transition partner and two transversion partners, so
random picking gives 1:2). Real bacterial genomes run Ti/Tv ≈ 1–2. Regenerated with
`-titv_ratio 2.0`, **same seed** (20260814). Measured: phiX 2.57 (only 50 SNPs, so noisy),
*E. coli* **2.04**.

**The mid-semester deck quotes the Ti/Tv 0.5 numbers**, so before regenerating, those results
and truth sets were archived to `results/archive/titv0.5_baseline/` with a README. Every figure
in that deck stays traceable.

An unexpected bonus: with the same seed, simuG placed **every SNP at the same position and
kept every indel identical** — verified by hashing positions and alleles. Only the SNP *alleles*
changed. So the regeneration is a perfectly controlled perturbation: same sites, same indels,
different mutation spectrum. Comparing the two baselines therefore isolates the effect of Ti/Tv
on the pipelines:

| | max \|ΔF1\| over all callers and both variant types |
|---|---|
| BWA-MEM | 0.0001 |
| minimap2 | 0.0002 |
| Bowtie2 | 0.0024 |

**No ranking changed.** The earlier claim — "Ti/Tv shouldn't bias the comparison, because every
pipeline sees the same truth" — was an argument; it is now a measurement.

### 8.2 The blocking defect: metrics were not part of the workflow

Six metric producers (`align_metrics.sh`, `placement_accuracy.py`, caller timing,
`extract_mean_q.py`, …) ran only by hand. `collect_results.py` joins them by `(tag, aligner)`,
and their tables held only the two baseline tags — so a sweep would have produced 3,960 rows
whose **feature columns were empty** (mapping rate, MAPQ, depth, placement accuracy, runtime,
memory, error rate) — precisely the features the model needs.

The fix is not "wrap each script in a rule". Every one of those scripts **appended** to a shared
TSV. Under Snakemake's parallel execution, two jobs appending to one file at once interleave
their writes and corrupt rows. The standard remedy is **scatter–gather**: each unit of work
writes its own one-row file (`work/metrics/<tag>.<aligner>.align.tsv`, …), and a single
aggregation rule concatenates them. Appending is safe serially and unsafe in parallel, and a
workflow engine *is* parallel.

The collector now runs with `--strict`, which **fails the workflow** if any feature cell is
empty, `nan` or `NA`. (The first version only checked for empty strings; a `nan` slipped through
on the first run — see 8.3 — which is why the check now treats all three as missing.)

### 8.3 Timing: sampling versus accounting

Runtime and peak memory are reported metrics, so they must be measured correctly.

**Attempt 1 — Snakemake's `benchmark:` directive.** It recorded wall time but every resource
column came back `NA`. The cause, found by reading Snakemake's source: it *samples* process memory
with `psutil` on a timer, and on macOS psutil generally cannot read another process's memory, so
every sample failed silently. Even where sampling works, a poll can miss a short memory peak.

**Attempt 2 — kernel accounting.** The operating system already records each process's
true high-water mark (`getrusage`), and `/usr/bin/time` reports it. `scripts/lib/measure.sh`
wraps exactly one tool invocation and records wall time, peak RSS and CPU time. It handles both
`time` dialects — BSD `-l` (RSS in bytes) and GNU `-f` (kilobytes) — because mixing those units
misreports memory by 1024×. Verified on a known 200 MB allocation (reported 243 MB: Python's
~40 MB baseline + 200 MB) and against the Phase 5 manual numbers (BWA 373 vs 382 MB; Bowtie2 67
vs 67 MB; minimap2 483 vs 491 MB).

Snakemake's own whole-job time is kept alongside as `job_seconds`. The difference is a uniform
~0.24 s of conda activation — small, and now visible rather than silently inside every number.

**Fairness (R8): who else is running?** A tool timed while seven other jobs share the CPU is
measured under arbitrary contention. So timed jobs claim a `machine` resource of 8 units against
a budget of 8 — they can only start when nothing else runs, and nothing else can start while they
do. That reproduces the conditions the Phase 5/6 baseline was timed under.

**But exclusivity for every run was too expensive.** While a *single-threaded* FreeBayes or
BCFtools job held the machine, seven cores sat idle with hundreds of short jobs queued behind it;
extrapolated, the sweep would have taken 6+ hours. The resolution separates two questions that
need different designs:

- **Accuracy** varies with the seed — that variance is the whole point — so it needs all 5 seeds.
- **Runtime** varies with the *condition* (depth, read length), not the seed. One clean
  measurement per condition is enough.

So only **seed 1** runs timed jobs exclusively (`config.yaml: timing_seeds`); seeds 2–5 run fully
parallel. Their timings are still recorded but flagged `timing_exclusive = no` and excluded from
runtime analysis — real measurements, honestly labelled, never presented as clean. Snakemake
supports this directly: a rule's `resources` can be a function of the wildcards.

A second, larger speed-up came from the scheduler. Snakemake's default solves an integer linear
program every scheduling round; with ~16,000 jobs that became the bottleneck. `scheduler: greedy`
in the profile cut the pilot's remaining 388 jobs to **78 seconds**.

### 8.4 Edge cases the baseline never exercised

The sweep reaches 5× coverage and a −10 quality shift, where call sets can be tiny or empty. Two
things were tested *before* running rather than assumed:

**What does vcfeval do with an empty call set?** It exits 0 and reports TP = 0, FN = every truth
variant, precision `NaN` (0/0). That is correct behaviour — and it exposed a latent bug in
`score_variants.sh`, which had a fallback that, *if vcfeval failed*, wrote a placeholder summary
with **TP = FP = FN = 0**. That fallback was written on the assumption that vcfeval errors on empty
input. It does not, so the fallback could only ever fire on a *genuine* failure — and would then
have **fabricated a result** (and FN = 0 would have been wrong even in its own terms). It now fails
loudly.

**F1 when precision is undefined.** vcfeval computes F1 from precision and recall, so with no calls
it reports F1 = `NaN`. The collector computes the standard count form, `F1 = 2TP / (2TP + FP + FN)`,
which equals the harmonic mean wherever that is defined and correctly gives 0 when TP = 0.
Precision is left as `NaN` — it genuinely is undefined.

**A silent-failure bug in the FreeBayes and BCFtools rules.** Both are pipelines
(`freebayes … | bgzip`). Run through plain `sh -c`, a pipeline's exit status is the *last*
command's, so a FreeBayes crash would let `bgzip` "succeed" on empty input and leave behind a
valid-looking empty VCF. They now run under `bash -o pipefail`, so any stage's failure fails the
job.

### 8.5 The filtered call sets were never built

Running the sweep-aware audit on the pilot found **zero hard-filtered call sets**. The
`hard_filter` rule existed in the Snakefile, but *no target ever requested its output*, so
Snakemake — which only builds what something asks for — never ran it. The old audit had passed
because it counted `*.filt.vcf.gz` files left behind by the Phase 6 shell scripts; clearing
`work/` exposed the gap. The brief requires both raw and filtered call sets, so the collector now
requests filtered sets and they are scored (primary single-run method). `results.tsv` carries both
under the `callset` column; the primary analysis uses `raw`, because the filtered set is one
operating point on the raw set's ROC curve.

**General lesson:** in a pull-based build system, a rule that nothing depends on is dead code
that looks alive. A check that counts files on disk cannot tell "the workflow built this" from
"something else left this here".

### 8.6 Disk: measured, then managed

Measured per 30× *E. coli* run: FASTQs 284 MB, ART truth SAM 315 MB, three BAMs 249 MB — 0.83 GB
retained (the HANDOFF's earlier 1.32 GB estimate guessed 250 MB per BAM; they are 83 MB). Across
the sweep's coverage-weighted total that is ~50 GB. FASTQs and the truth SAM are now `temp()`:
Snakemake deletes them once their last consumer has run (three aligners and `read_metrics` for
the reads; three placement jobs for the truth SAM). ART is deterministic for a fixed seed, so
anything deleted is exactly regenerable. Retained steady state: ~15 GB.

### 8.7 Repository hygiene at sweep scale

The sweep writes ~16,000 per-job logs, ~1,300 benchmark files and ~5,000 vcfeval directories.
After the `.snakemake/` incident (Phase 7b), these are gitignored up front: `logs/run/`,
`logs/ploidy/`, `benchmarks/`, `results/vcfeval/`. What is committed is the distilled output —
`results/*.tsv`, `results/analysis/`, `results/model/`, `results/figures/` — plus summary records
(`logs/ploidy_verification.txt`, versions, verification reports). The Phase 5–7 per-run logs for
the Ti/Tv 0.5 run were moved to `logs/legacy_titv0.5/`; they remain in git history.

### 8.8 The pilot gate

Before 990 runs, a pilot ran the harshest corners: *E. coli* at 5× (all 5 seeds), the −10 quality
shift on both genomes, phiX at 5×, plus both baselines. The sweep-aware audit
(`VERIFY_SCOPE=present`) passed every check except the filtered-set gap in 8.5, which was fixed
before launching. Coverage landed within ±2% of the request for every run, every timed job had
both wall time and peak RSS, and all 90 call sets were haploid.

### 8.9 Three problems that only appeared at sweep scale

**A SIGPIPE race in `align_metrics`.** Two of ~90 jobs failed with exit status 141 (= 128 + 13,
SIGPIPE). The script parsed `samtools flagstat` output with `echo "$FS" | awk '/primary/{print;
exit}'`. awk exits the moment it finds its line, closing the pipe; if `echo` is still writing it is
killed by SIGPIPE, `pipefail` propagates it and `set -e` aborts. A race — never seen in the pilot,
twice in the sweep. It could only cause a *failure*, never a wrong value (a successful run has read
the correct line). Fixed with here-strings (`awk '…' <<< "$FS"` — no pipe, no SIGPIPE). The script
was replaced by write-then-rename rather than edited in place, because bash reads scripts
incrementally and a running job could otherwise execute a half-old, half-new file. This is the
same failure class as the `bcftools | head` bug in Phase 2: **any pipe whose reader can exit early
is unsafe under `pipefail`.**

**A thread reservation is not thread usage.** Mid-sweep the machine was 36% idle with only two
GATK jobs running, each at ~100% CPU. Each *reserved* 4 cores (`--native-pair-hmm-threads 4`), but
HaplotypeCaller is single-threaded apart from short PairHMM bursts — so two jobs held the whole
8-core budget while using about 2. Lowering GATK's threads would have fixed the accounting but
changed its command, making new timings inconsistent with every seed-1 timing already taken. The
fix was the budget instead: `--cores 12`. This cannot compromise clean timing, because exclusivity
is enforced by the separate `machine` resource (8 of 8), not by cores.

**`kill -INT` did nothing.** Restarting the sweep needed the running Snakemake stopped, and SIGINT
was silently ignored. POSIX rule: a command started with `&` from a *non-interactive* shell runs
with SIGINT set to "ignore", and that disposition survives `exec`; Python, seeing SIGINT already
ignored at start-up, never installs its KeyboardInterrupt handler. SIGTERM worked, and Snakemake
cleaned up its in-flight jobs. Afterwards the process tree was checked for orphaned job processes —
an orphan still writing an output while the restarted run scheduled the same job would have been a
genuine race.

A smaller one: adding the report rules to `rule all` made the restart fail immediately, because the
report template did not exist yet and Snakemake validates the whole DAG before running anything.
The sweep was relaunched with explicit targets (`results/results.tsv`,
`logs/ploidy_verification.txt`).

**External processes and clean timing.** The `machine` resource only excludes *Snakemake's own*
jobs from running beside a timed job; it cannot see anything else on the computer. Work run outside
the workflow during the sweep (model fitting, figure generation, a conda install) was therefore
checked against the log: none of it overlapped a seed-1 timed job, and nothing CPU-heavy was run
outside Snakemake for the remainder of the sweep.

### 8.10 What the sweep found, in plain language

The full numbers are in `docs/FINAL_REPORT.pdf`; this is the shape of the answer.

- **Everything is accurate, and the differences are still real.** On *E. coli* at 30× every
  pipeline has SNV F1 between 0.989 and 0.995. That sounds like "they're all the same", but the
  seed-to-seed standard deviation is about 0.0005, so a gap of 0.005 is about ten standard
  deviations: a reproducible difference, not noise. This is why five seeds were worth running:
  with one seed (Phase 7) there was no way to tell.
- **The aligner matters more than the caller.** In a blocked two-way ANOVA the aligner explains a
  median 78% of the between-pipeline variation in SNV F1 and the caller 2%. Most of that is one
  aligner, Bowtie2 (8.11).
- **Depth is the only data property that changes the answer.** Below 10× the ranking reorders:
  BCFtools, the least conservative caller, wins SNVs because it calls on thinner evidence.
  Raising the error rate seven-fold barely moves F1. Independent random errors rarely agree on the
  same wrong base at the same site, so callers filter them easily. The errors that hurt are
  *systematic* ones from alignment, which repeat across reads.
- **A model adds little above 10×.** One pipeline family (BWA-MEM or minimap2 with FreeBayes or
  GATK) is within 0.001 of the best everywhere, so "always use BWA-MEM + GATK" is nearly as good as
  any prediction. That is itself a finding, and it is reported as one rather than dressed up.

### 8.11 Testing explanations instead of telling stories

A surprising pattern invites a plausible story. The F1-vs-coverage figure showed Bowtie2 with
FreeBayes or BCFtools getting *worse* with more data. The story wrote itself: Bowtie2 aligns
end-to-end, cannot soft-clip, so reads ending just past an indel get mismatches instead of a gap,
and with more depth those mismatches become confident false SNVs. A story is not evidence, so
`scripts/diagnose_errors.py` turns each one into a test with a number that could have come out
wrong:

| Story | Test | Result |
|---|---|---|
| Bowtie2's false SNVs are indel artefacts | distance from each FP SNV to the nearest *true* indel, vs the fraction of genome that close by chance | 99% within 150 bp at 100× (median 12 bp) vs 6.3% background |
| …and they grow with depth | same, at every coverage | 24 → 137 FPs from 5× to 100× (Bowtie2 + FreeBayes) |
| BCFtools' slow F1 decline above 30× is the same thing, milder | same test, BWA-MEM + BCFtools | 5 → 21 FPs, median 1 bp from a true indel |
| the ~50 variants nobody finds are in repeats | share of reads with MAPQ ≥ 20 at missed sites vs all true sites | 96–100% of misses are low-MAPQ sites vs 1.4–2.3% of all sites |
| phiX misses are a coverage artefact | list every phiX error by position | every miss is position 51, in the first read length of a linearised circular genome |

Two details are worth noticing. At 5× the false calls are *not* near indels (6% vs 6.3%
background): that is a different mechanism, low-depth calls on sequencing errors. Without the
test the indel explanation would have been wrongly applied there too. And Bowtie2 + GATK misses
about twice as many variants as other GATK pipelines, because Bowtie2 gives low MAPQ at more
sites and GATK, unlike FreeBayes (threshold 1) and BCFtools (0), drops reads below MAPQ 20 by
default. That last link is inferred from the tools' documented defaults, not tested by changing
the threshold, and HANDOFF says so.

The strongest test, rerunning Bowtie2 with `--local` and watching the artefacts disappear, was
not run. A test that shows a prediction holding is weaker than an intervention that removes the
cause; that experiment is first on the Phase 2 list.

### 8.12 A bug in the model's scoring: ties broken alphabetically

The first recommendation table claimed the decision tree's pick for 5× indels was Bowtie2 +
BCFtools. That is the *worst* pipeline at 5×. The cause was not the model but the code scoring it.

A regression tree predicts one value per leaf. Any pipelines that land in the same leaf get the
same prediction: they are tied, and the model is indifferent between them. The scoring code then
chose `idxmax()` of the predictions, and pandas' `idxmax` returns the *first* maximum. The rows
were sorted alphabetically, so every tie went to whichever pipeline name sorts first:
`bowtie2+bcftools`. The model was being blamed, or credited, for choices it never made.

The fix is to score what the model actually says. When k pipelines tie, the model's choice is a
uniform random pick among them, so its regret is the *expected* regret: best F1 minus the mean F1
of the tied set. The recommendation table now reports the tied set itself ("BWA-MEM/minimap2 +
GATK/FreeBayes, 4 tied"). This changed the conclusions. The tree's held-out-seed regret went from
worse than the trivial rule (0.00155 vs 0.00093) to slightly better (0.00086). The true best
pipeline turned out to be inside the tree's set in all 22 condition × type cases.

**The general lesson: any argmax over a model's output needs an explicit tie policy.** Tree
ensembles average away most ties, but a single tree, a rule list or a rounded score produces
them constantly, and `idxmax`/`argmax` will silently resolve them by row order.

### 8.13 What the tree's "error rate" splits really mean

The fitted tree splits on measured error rate twice, and neither split means what it says.

- **Under 5× indels**, "error ≤ 0.19%" separates *seeds of one condition*. 5× was only simulated
  at the baseline error, so its measured error varies only from read set to read set (about 0.5%
  relative). The split fits seed-level noise. The figure labels it as noise; it was not removed,
  because re-tuning the tree after seeing an embarrassing split is exactly the kind of
  after-the-fact adjustment that makes models look better than they are.
- **"Error ≤ 0.15%" elsewhere is "the 75 bp reads".** ART's error rate rises along a read, so
  shorter reads have a lower *mean* error (0.145% vs 0.195%). In a one-factor-at-a-time design,
  read length and error rate are therefore partly confounded, and the tree can use either to
  separate the 75 bp condition. The forest's separate importances for the two should not be
  compared with each other.

The figure code works out what each error-rate split separates from the training rows: one
condition's seeds, or exactly one condition. It does not hard-code the explanation.

### 8.14 Building the report so it cannot drift from the data

Every table in `docs/FINAL_REPORT.md` and every number quoted in its prose is generated by
`scripts/build_report.py` from `results/`. The prose template uses `{{T_...}}` placeholders for
tables and `{{V_...}}` for numbers, and the builder refuses to write a report with an unknown or
unfilled placeholder. A hand-typed number in a report is a number that will silently go stale the
first time anything is re-run.

The renderer surfaced four problems worth knowing:

- **pandas attribute traps, twice.** `x.cov == 5` compared the DataFrame's `.cov()` *method* to 5
  and silently matched nothing, so the alignment table came out empty. Later `g.pipe == p` did
  the same with `.pipe()`. Any column whose name is also a DataFrame method must be accessed as
  `df["name"]`.
- **Equal-width table columns.** Pandoc sizes a wrapping pipe table's columns by the *relative
  number of dashes* in the separator row, so the conventional `|---|---|` makes every column the
  same width. The 22-row recommendation table spanned 2.5 pages. Sizing dashes by content, floored
  at each column's longest unbreakable word and capped so one long column cannot starve the rest,
  brought the report from 23 to 21 pages without dropping anything.
- **Silent missing glyphs.** XeLaTeX with Helvetica has no `→`. It does not fail; it drops the
  character and writes "Missing character" to a log nobody reads. The fix maps `→` to a math arrow
  and uses Menlo for code, which has `≥`. The render was checked by counting those warnings
  (zero) and by looking at every page.
- **Tables reported to the reader must be honest about scope.** For example, Table 5's
  false-positive counts are seed 1 only; the caption says so.

### 8.15 Corrections to earlier claims

- The mid-semester deck's numbers come from the Ti/Tv 0.5 truth sets at seed 1. They remain
  traceable in `results/archive/titv0.5_baseline/` but are superseded by the report.
- "All nine pipelines are perfect on phiX" was true for seed 1 at baseline. Across five seeds,
  16 of 18 baseline cells are perfect, and every exception is explained (8.11).
- Phase 7's question "is the SNV spread larger than noise?" is answered: yes, by about ten
  standard deviations.
