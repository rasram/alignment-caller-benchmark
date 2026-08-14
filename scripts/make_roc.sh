#!/usr/bin/env bash
# PHASE 7a — ROC curves for the 9 pipelines, SNV and indel separately.
#
# WHY THE ROC MATTERS MORE THAN THE F1 TABLE
# A single F1 number is one operating point: it depends on whatever QUAL
# threshold happened to be applied. Since QUAL is not calibrated identically
# across GATK, FreeBayes and BCFtools (NOTES 6.7), comparing single F1 values
# partly compares the callers' QUAL scales rather than their ability to find
# variants. The ROC sweeps the threshold across its entire range, so it shows the
# whole precision/recall trade-off and removes that arbitrariness.
#
# Curves come from the SINGLE-RUN vcfeval output (per NOTES 7.6), using RTG's own
# snp_roc / non_snp_roc files rather than pre-split call sets.
#
# Usage: bash scripts/make_roc.sh [tag]

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
VE="$REPO/results/vcfeval"
OUTDIR="$REPO/results"
mkdir -p "$OUTDIR"

TAGS="${*:-ecoli_cov30_len150_err0_seed1 phiX_cov30_len150_err0_seed1}"

for TAG in $TAGS; do
  for TYPE in snp non_snp; do
    label=$([ "$TYPE" = "snp" ] && echo SNV || echo INDEL)
    args=()
    for ALN in bwa bowtie2 minimap2; do
      for CAL in gatk freebayes bcftools; do
        f="$VE/${TAG}__${ALN}__${CAL}__raw__all/${TYPE}_roc.tsv.gz"
        [[ -s "$f" ]] && args+=(--curve "${f}=${ALN}+${CAL}")
      done
    done
    if [[ ${#args[@]} -eq 0 ]]; then
      echo "  no ROC data for $TAG $label"; continue
    fi
    out="$OUTDIR/roc_${TAG%%_*}_${label}.svg"
    # rtg rocplot refuses to overwrite an existing file, so the script would fail
    # on every re-run rather than regenerating.
    rm -f "$out"
    if "$RTG" rocplot --svg "$out" --title "${TAG%%_*} — ${label} — 30x 150bp PE" \
         "${args[@]}" > /dev/null 2>&1; then
      echo "  wrote $out  ($((${#args[@]}/2)) curves)"
    else
      echo "  FAILED to plot $TAG $label"
    fi
  done
done
