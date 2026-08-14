#!/usr/bin/env bash
# PHASE 2 — decide how the truth VCF should carry genotypes.
#
# The brief offers two options and says to test both on phiX:
#   (a) add a FORMAT/sample column with GT=1
#   (b) leave simuG's 8-column output and use `vcfeval --squash-ploidy`
#
# This script runs the experiment. It builds synthetic "call sets" from the truth
# itself (so the variant content is identical and ONLY the genotype encoding
# differs), then scores each with vcfeval. Any deviation from F1=1.0 is therefore
# caused purely by genotype representation, not by variant-calling accuracy.
#
# Call sets tested:
#   calls_hap  GT=1     what a correctly-configured haploid caller emits
#   calls_dip  GT=1/1   what a diploid-defaulted caller emits at a clean site
#   calls_het  GT=0/1   what a diploid-defaulted caller emits when error/noise
#                       makes the ALT look like only half the reads (rule R2)
#
# Usage: bash scripts/test_genotype_encoding.sh
# Writes: logs/genotype_decision.txt

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
BCFTOOLS="$CONDA_BASE/envs/callers/bin/bcftools"
BGZIP="$CONDA_BASE/envs/callers/bin/bgzip"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
GEN=phiX
TRUTH="$REPO/data/truth/${GEN}.truth.vcf.gz"
NOGT="$REPO/data/truth/${GEN}.truth.nogt.vcf.gz"
SDF="$REPO/data/refs/${GEN}.sdf"
BED="$REPO/data/truth/${GEN}.confident.bed"
LOG="$REPO/logs/genotype_decision.txt"
W="$REPO/work/gt_test"; rm -rf "$W"; mkdir -p "$W"

mk() { # name gt_value
  "$BCFTOOLS" view "$TRUTH" \
    | awk -v g="$2" 'BEGIN{OFS="\t"} /^#/{print;next} {$10=g; print}' \
    | "$BGZIP" > "$W/$1.vcf.gz"
  "$BCFTOOLS" index -t -f "$W/$1.vcf.gz"
}

# vcfeval summary.txt columns (verified against real output):
#   1 Threshold  2 True-pos-baseline  3 True-pos-call  4 False-pos
#   5 False-neg  6 Precision          7 Sensitivity    8 F-measure
# Note there are TWO true-positive columns: TP-baseline counts matched TRUTH
# records, TP-call counts the CALL records that matched them. They differ when
# one truth variant is represented by several call records (or vice versa).
score() { # label baseline calls [extra flags...]
  local lab="$1" b="$2" c="$3"; shift 3
  # Labels contain '/' and spaces, which cannot go in a path — slugify them.
  local slug; slug="$(echo "$lab" | tr -c 'A-Za-z0-9._-' '_')"
  rm -rf "$W/out_$slug"
  if "$RTG" vcfeval -b "$b" -c "$c" -t "$SDF" -e "$BED" -o "$W/out_$slug" "$@" \
       > "$W/$slug.log" 2>&1; then
    tail -1 "$W/out_$slug/summary.txt" | awk -v l="$lab" \
      '{printf "  %-22s TP=%-5s FP=%-5s FN=%-5s precision=%-8s recall=%-8s F1=%s\n", l, $2, $4, $5, $6, $7, $8}'
  else
    printf "  %-22s REFUSED TO RUN: %s\n" "$lab" \
      "$(grep -iE '^Error' "$W/$slug.log" | head -1 | cut -c1-64)"
  fi
}

{
echo "# Genotype-encoding decision experiment — $GEN — $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "# Variant content is IDENTICAL across all call sets; only GT encoding differs."
echo

mk calls_hap "1"
mk calls_dip "1/1"
mk calls_het "0/1"

echo "## Option (a): truth HAS a sample column with GT=1"
score "calls GT=1"      "$TRUTH" "$W/calls_hap.vcf.gz"
score "calls GT=1/1"    "$TRUTH" "$W/calls_dip.vcf.gz"
score "calls GT=0/1"    "$TRUTH" "$W/calls_het.vcf.gz"
echo
echo "## Option (a) + --squash-ploidy"
score "GT=1/1 squashed"  "$TRUTH" "$W/calls_dip.vcf.gz" --squash-ploidy
score "GT=0/1 squashed"  "$TRUTH" "$W/calls_het.vcf.gz" --squash-ploidy
echo
echo "## Option (b): truth has NO sample column (raw simuG shape)"
score "no-GT baseline"           "$NOGT" "$W/calls_hap.vcf.gz"
score "no-GT + squash-ploidy"    "$NOGT" "$W/calls_hap.vcf.gz" --squash-ploidy
echo
cat << 'EOF'
## DECISION

Option (b) is not available. vcfeval REFUSES a baseline with no sample column,
with or without --squash-ploidy:
    Error: Record did not contain enough samples
So the truth VCF must carry a FORMAT/sample column. Option (a) it is: GT=1.

Chosen scoring mode: option (a) with DEFAULT genotype-aware matching,
i.e. NO --squash-ploidy.

Reasoning — look at the GT=0/1 row above. Against the haploid truth it scores
F1 = 0.0000 (every variant counted as both FP and FN), but WITH --squash-ploidy
the same file scores F1 = 1.0000.

That is exactly the rule-R2 failure mode. If a caller's ploidy flag were dropped
and it emitted heterozygous genotypes, --squash-ploidy would hide the mistake
completely and report a perfect score. Default matching makes the failure loud
and impossible to miss.

Note also the GT=1/1 row: it matches the haploid truth perfectly even WITHOUT
squashing, because vcfeval treats homozygous-ALT as equivalent to haploid ALT.
So vcfeval will NOT catch every ploidy misconfiguration on its own. This is why
rule R2 requires grepping the GT field directly (Phase 6) rather than trusting
the score.
EOF
} | tee "$LOG"

echo
echo "Wrote $LOG"
