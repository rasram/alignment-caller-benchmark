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
