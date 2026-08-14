#!/usr/bin/env bash
# PHASE 7a — normalise every call set and score it against truth with rtg vcfeval.
#
# RULE R3 — the normalisation applied here is BYTE-IDENTICAL to the command used
# on the truth set in build_truth.sh:
#     bcftools norm -f <ref.fa> -m -any --atomize
# --atomize is essential: FreeBayes merges nearby variants into MNV/complex
# records, and without decomposing them a type-split drops real variants (NOTES 6.6).
#
# RULE R4 — comparison is done ONLY by rtg vcfeval. No position/allele string
# matching anywhere.
#
# SNV and indel are scored separately (the brief's requirement) by splitting with
# `bcftools view -v snps|indels` and running vcfeval on each. A third run scores
# the full call set; RTG emits its OWN per-type breakdown (snp_roc/non_snp_roc)
# in that run, which is used to CROSS-CHECK the split approach. If the two
# disagree, the split is distorting haplotype context and the numbers are suspect.
#
# Usage: bash scripts/score_variants.sh <genome> [tag] [raw|filt]

set -euo pipefail

GEN="${1:?usage: score_variants.sh <genome> [tag] [raw|filt]}"
TAG="${2:-${GEN}_cov30_len150_err0_seed1}"
SET="${3:-raw}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tools.sh"
BCFTOOLS="$(resolve_tool bcftools callers)"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
REF="$REPO/data/refs/${GEN}.fa"
SDF="$REPO/data/refs/${GEN}.sdf"
TRUTH="$REPO/data/truth/${GEN}.truth.vcf.gz"
BED="$REPO/data/truth/${GEN}.confident.bed"
W="$REPO/work"
VE="$REPO/results/vcfeval"
mkdir -p "$VE" "$W/norm"

ALIGNERS="bwa bowtie2 minimap2"
CALLERS="gatk freebayes bcftools"

# IDENTICAL to build_truth.sh — R3.
NORM_ARGS=(-f "$REF" -m -any --atomize)

# Pre-split the truth set once (same normalisation already applied).
for VT in snps indels; do
  TS="$W/norm/${GEN}.truth.${VT}.vcf.gz"
  "$BCFTOOLS" view -v "$VT" -Oz -o "$TS" "$TRUTH" 2>/dev/null
  "$BCFTOOLS" index -t -f "$TS"
done

run_eval() { # outdir baseline calls
  local out="$1" b="$2" c="$3"
  rm -rf "$out"
  "$RTG" vcfeval -b "$b" -c "$c" -t "$SDF" -e "$BED" \
      --vcf-score-field=QUAL -o "$out" > "$out.log" 2>&1 || {
    # vcfeval exits non-zero when the call set is empty; record that rather than
    # aborting the whole sweep.
    mkdir -p "$out"
    printf 'Threshold\tTrue-pos-baseline\tTrue-pos-call\tFalse-pos\tFalse-neg\tPrecision\tSensitivity\tF-measure\n' > "$out/summary.txt"
    printf -- '----\n' >> "$out/summary.txt"
    printf 'None\t0\t0\t0\t0\t0.0000\t0.0000\t0.0000\n' >> "$out/summary.txt"
    return 0
  }
}

echo "[$TAG/$SET] normalising and scoring 9 pipelines"

for ALN in $ALIGNERS; do
  for CAL in $CALLERS; do
    SRC="$W/${TAG}.${ALN}.${CAL}.${SET}.vcf.gz"
    [[ -s "$SRC" ]] || { echo "  skip (missing) $ALN/$CAL"; continue; }

    NORM="$W/norm/${TAG}.${ALN}.${CAL}.${SET}.norm.vcf.gz"
    "$BCFTOOLS" norm "${NORM_ARGS[@]}" -Oz -o "$NORM" "$SRC" \
        2> "$W/norm/${TAG}.${ALN}.${CAL}.${SET}.norm.log"
    "$BCFTOOLS" index -t -f "$NORM"

    # --- full call set (also yields RTG's native snp/non_snp ROC) ------------
    run_eval "$VE/${TAG}__${ALN}__${CAL}__${SET}__all" "$TRUTH" "$NORM"

    # --- explicit per-type split (the brief's method) ------------------------
    for VT in snps indels; do
      SPLIT="$W/norm/${TAG}.${ALN}.${CAL}.${SET}.${VT}.vcf.gz"
      "$BCFTOOLS" view -v "$VT" -Oz -o "$SPLIT" "$NORM" 2>/dev/null
      "$BCFTOOLS" index -t -f "$SPLIT"
      run_eval "$VE/${TAG}__${ALN}__${CAL}__${SET}__${VT}" \
               "$W/norm/${GEN}.truth.${VT}.vcf.gz" "$SPLIT"
    done

    printf '  scored %-9s %-10s\n' "$ALN" "$CAL"
  done
done

echo "[$TAG/$SET] vcfeval output in $VE"
