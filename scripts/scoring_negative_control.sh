#!/usr/bin/env bash
# PHASE 7a — negative control for the scoring machinery.
#
# WHY THIS EXISTS
# All nine phiX pipelines score F1 = 1.0000. The brief warns that an F1 of
# exactly 1.0 usually indicates a bug rather than a result — a scoring setup that
# is accidentally comparing a file against itself, or matching everything, would
# produce the same number.
#
# A perfect score is only believable if the same machinery can be shown to FAIL
# when it should. This deliberately corrupts a call set in three different ways
# and confirms vcfeval punishes each one:
#
#   1. positions shifted by +5 bp   -> right variants, wrong places
#   2. ALT alleles mutated          -> right places, wrong alleles
#   3. half the calls deleted       -> recall must drop to ~0.5, precision stay 1.0
#
# If any of these still scored 1.0, the F1 = 1.0 results would be meaningless.
#
# Usage: bash scripts/scoring_negative_control.sh
# Writes: logs/scoring_negative_control.txt

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
C="$CONDA_BASE/envs/callers/bin"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
GEN=phiX
TAG=phiX_cov30_len150_err0_seed1
TRUTH="$REPO/data/truth/${GEN}.truth.vcf.gz"
SDF="$REPO/data/refs/${GEN}.sdf"
BED="$REPO/data/truth/${GEN}.confident.bed"
GOOD="$REPO/work/norm/${TAG}.bwa.gatk.raw.norm.vcf.gz"
LOG="$REPO/logs/scoring_negative_control.txt"
W="$REPO/work/negctl"; rm -rf "$W"; mkdir -p "$W"

score() { # label calls
  local lab="$1" calls="$2" slug
  slug="$(echo "$lab" | tr -c 'A-Za-z0-9._-' '_')"
  rm -rf "$W/out_$slug"
  if "$RTG" vcfeval -b "$TRUTH" -c "$calls" -t "$SDF" -e "$BED" \
        --vcf-score-field=QUAL -o "$W/out_$slug" > "$W/$slug.log" 2>&1; then
    tail -1 "$W/out_$slug/summary.txt" | awk -v l="$lab" \
      '{printf "  %-34s TP=%-5s FP=%-5s FN=%-5s precision=%-8s recall=%-8s F1=%s\n", l, $2, $4, $5, $6, $7, $8}'
  else
    printf "  %-34s vcfeval failed\n" "$lab"
  fi
}

{
echo "# Scoring negative control — $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "# Baseline: phiX truth (50 SNP + 10 indel). Calls: bwa/gatk, normalised."
echo

# 0. the real, uncorrupted call set
score "0. unmodified (expect F1=1.0)" "$GOOD"

# 1. shift every position by +5 bp
"$C/bcftools" view "$GOOD" \
  | awk 'BEGIN{OFS="\t"} /^#/{print;next} {$2=$2+5; print}' \
  | "$C/bgzip" > "$W/shift.vcf.gz"
"$C/bcftools" index -t -f "$W/shift.vcf.gz" 2>/dev/null
score "1. positions +5bp (expect ~0)" "$W/shift.vcf.gz"

# 2. change every ALT base to a DIFFERENT base (SNVs only; indels left alone).
#    A naive A->C->G->T->A rotation is not safe: it can land on the REF base
#    (e.g. REF=A ALT=T rotates to A), and vcfeval correctly rejects such a record
#    with "ALT allele is the same as REF" rather than scoring it. So step through
#    the cycle until the result differs from REF.
"$C/bcftools" view "$GOOD" \
  | awk 'BEGIN{OFS="\t"; split("A,C,G,T",B,",")}
      /^#/{print;next}
      { if (length($4)==1 && length($5)==1) {
          for (i=1; i<=4; i++) {
            if (B[i] != $5 && B[i] != $4) { $5 = B[i]; break }
          }
        }
        print }' \
  | "$C/bgzip" > "$W/altmut.vcf.gz"
"$C/bcftools" index -t -f "$W/altmut.vcf.gz" 2>/dev/null
score "2. ALT alleles rotated (expect low)" "$W/altmut.vcf.gz"

# 3. keep only every second call
"$C/bcftools" view "$GOOD" \
  | awk 'BEGIN{OFS="\t"} /^#/{print;next} {n++; if(n%2==1) print}' \
  | "$C/bgzip" > "$W/half.vcf.gz"
"$C/bcftools" index -t -f "$W/half.vcf.gz" 2>/dev/null
score "3. half the calls (expect rec~0.5)" "$W/half.vcf.gz"

echo
echo "INTERPRETATION"
echo "  If rows 1-3 score well below 1.0 while row 0 scores 1.0, the scoring"
echo "  pipeline genuinely discriminates and the perfect phiX results are real:"
echo "  phiX is 5,386 bp with no repetitive sequence, so at 30x coverage with"
echo "  150 bp reads every variant is unambiguously recoverable."
} | tee "$LOG"

echo
echo "Wrote $LOG"
