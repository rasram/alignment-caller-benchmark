#!/usr/bin/env bash
# PHASE 5 — align simulated reads with all three aligners.
#
# RULE R1: reads are aligned to the ORIGINAL reference (data/refs/<g>.fa), never
#          to the mutated genome they were simulated from. Asserted below.
# RULE R6: every BAM carries a read group. GATK refuses to run without one, and
#          the three aligners spell read groups differently (see below).
# RULE R8: identical thread count (4) for all three aligners, since runtime is a
#          reported metric. No trimming; all three read the same FASTQ bytes.
#
# TIMING METHODOLOGY
# The aligner is timed ALONE, writing SAM to disk, rather than piped into
# `samtools sort`. Piping is the normal production idiom, but sorting costs
# roughly the same for all three aligners (same read count) and would dilute the
# very difference we are trying to measure. Sorting and indexing happen after,
# untimed. This is applied identically to all three, so the comparison is fair.
#
# Usage: bash scripts/align_reads.sh <genome> [tag] [threads]
#        bash scripts/align_reads.sh phiX phiX_cov30_len150_err0_seed1 4

set -euo pipefail

GEN="${1:?usage: align_reads.sh <genome> [tag] [threads]}"
TAG="${2:-${GEN}_cov30_len150_err0_seed1}"
THREADS="${3:-4}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tools.sh"
BWA="$(resolve_tool bwa align)"
BOWTIE2="$(resolve_tool bowtie2 align)"
MINIMAP2="$(resolve_tool minimap2 align)"
SAMTOOLS="$(resolve_tool samtools align)"
GATK="$(resolve_tool gatk callers)"
GATK_DIR="$(dirname "$GATK")"
REF="$REPO/data/refs/${GEN}.fa"
BT2IDX="$REPO/data/refs/${GEN}"
W="$REPO/work"
LOGS="$REPO/logs"
mkdir -p "$W" "$LOGS"

R1="$W/${TAG}_1.fq"
R2="$W/${TAG}_2.fq"
for f in "$R1" "$R2" "$REF"; do
  [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 1; }
done

# --- R1 guard ---------------------------------------------------------------
# Aligning to the mutated genome would make every pipeline find zero variants.
if cmp -s "$REF" "$REPO/data/truth/${GEN}.simseq.genome.fa"; then
  echo "FATAL: reference is identical to the mutated genome — R1 violated." >&2
  exit 1
fi

# --- Read group (R6) --------------------------------------------------------
# ID=s1 read-group id, SM=sim sample name, PL=ILLUMINA platform, LB=lib1 library.
# GATK requires SM; MarkDuplicates uses LB to decide what counts as a duplicate.
RG_ID="s1"; RG_SM="sim"; RG_PL="ILLUMINA"; RG_LB="lib1"
RG_STR="@RG\tID:${RG_ID}\tSM:${RG_SM}\tPL:${RG_PL}\tLB:${RG_LB}"

TIMEFILE_DIR="$LOGS/timing"; mkdir -p "$TIMEFILE_DIR"

# BSD /usr/bin/time -l reports max RSS in BYTES (GNU time -v reports kilobytes).
# Getting this backwards inflates or deflates memory by 1024x.
parse_time() { # timefile -> "seconds<TAB>peak_rss_mb"
  local tf="$1"
  awk '
    / real / { secs=$1 }
    /maximum resident set size/ { rss=$1 }
    END { printf "%.2f\t%.1f", secs, rss/1048576 }
  ' "$tf"
}

run_aligner() {
  local name="$1"; shift
  local sam="$W/${TAG}.${name}.sam"
  local tf="$TIMEFILE_DIR/${TAG}.${name}.time"
  echo "  [$name] aligning with $THREADS threads..."
  /usr/bin/time -l "$@" > "$sam" 2> "$tf" || {
    echo "FATAL: $name failed; see $tf" >&2; tail -20 "$tf" >&2; exit 1; }
  local stats; stats=$(parse_time "$tf")
  echo -e "  [$name] wall=$(echo "$stats" | cut -f1)s peak_rss=$(echo "$stats" | cut -f2)MB"
  printf '%s\t%s\t%s\n' "$TAG" "$name" "$stats" >> "$LOGS/align_timing.tsv"
}

echo "[$TAG] aligning to ORIGINAL reference: $REF"

# --- 1. BWA-MEM -------------------------------------------------------------
run_aligner bwa "$BWA" mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2"

# --- 2. Bowtie2 -------------------------------------------------------------
# Bowtie2 will not accept a single @RG string; it takes the ID separately via
# --rg-id and each additional field as its own --rg. Passing the tab-delimited
# @RG string here produces a malformed header that GATK later rejects.
run_aligner bowtie2 "$BOWTIE2" -p "$THREADS" \
  --rg-id "$RG_ID" --rg "SM:${RG_SM}" --rg "PL:${RG_PL}" --rg "LB:${RG_LB}" \
  -x "$BT2IDX" -1 "$R1" -2 "$R2"

# --- 3. minimap2 ------------------------------------------------------------
# -ax sr = short-read preset. Without it minimap2 uses long-read defaults and
# places short reads badly.
run_aligner minimap2 "$MINIMAP2" -ax sr -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2"

# --- sort, mark duplicates, index (untimed; identical for all three) --------
for aln in bwa bowtie2 minimap2; do
  sam="$W/${TAG}.${aln}.sam"
  sorted="$W/${TAG}.${aln}.sorted.bam"
  md="$W/${TAG}.${aln}.md.bam"

  echo "  [$aln] sort -> markdup -> index"
  "$SAMTOOLS" sort -@ "$THREADS" -o "$sorted" "$sam" 2> "$LOGS/sort_${TAG}_${aln}.log"
  "$SAMTOOLS" index "$sorted"

  # Simulated reads contain no PCR duplicates, so this marks ~0%. Run anyway for
  # pipeline realism; its inertness is reported rather than hidden (see NOTES).
  #
  # NOTE: `gatk` is a Python launcher script, not a binary. Invoking it by
  # absolute path is NOT enough — it calls `env python`, so its environment's
  # bin/ must be on PATH or it dies with "env: python: No such file or
  # directory". Hence the PATH prefix here (and everywhere gatk is used).
  PATH="$GATK_DIR:$PATH" "$GATK" MarkDuplicates \
      -I "$sorted" -O "$md" -M "$LOGS/${TAG}.${aln}.md.metrics" \
      --VALIDATION_STRINGENCY LENIENT \
      > "$LOGS/markdup_${TAG}_${aln}.log" 2>&1
  "$SAMTOOLS" index "$md"

  rm -f "$sam"        # uncompressed SAM is large and fully regenerable
done

echo "[$TAG] done. BAMs: $W/${TAG}.{bwa,bowtie2,minimap2}.md.bam"
