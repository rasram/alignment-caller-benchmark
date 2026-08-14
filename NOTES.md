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
