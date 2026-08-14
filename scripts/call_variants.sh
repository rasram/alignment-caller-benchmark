#!/usr/bin/env bash
# PHASE 6 — call variants with all three callers on all three BAMs -> 9 VCFs.
#
# RULE R2 — PLOIDY. phiX and E. coli are HAPLOID. All three callers default to
# diploid and will emit 0/1 and 1/1 genotypes with NO error message. Each caller
# spells the flag differently:
#     GATK       --sample-ploidy 1
#     FreeBayes  -p 1
#     BCFtools   bcftools call --ploidy 1
# The GT field is grepped and logged after every run (logs/ploidy_verification.txt).
# As Phase 2 proved, vcfeval will NOT reliably catch this: a caller emitting 1/1
# scores a perfect F1 against haploid truth. Only the direct GT check catches it.
#
# RULE R7 — NO BQSR. GATK Best Practices includes base quality score
# recalibration, which needs a database of known variants to tell real variation
# apart from systematic error. No such database exists for E. coli or phiX.
# Bootstrapping one would hand GATK a preprocessing step the other two callers do
# not get, confounding the comparison. Skipped deliberately.
#
# RULE R8 — identical treatment. Same BAMs, same reference, same filter
# expression for every caller. QUAL is retained in the raw callset for ROC.
#
# Usage: bash scripts/call_variants.sh <genome> [tag] [threads]

set -euo pipefail

GEN="${1:?usage: call_variants.sh <genome> [tag] [threads]}"
TAG="${2:-${GEN}_cov30_len150_err0_seed1}"
THREADS="${3:-4}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
C="$CONDA_BASE/envs/callers/bin"
REF="$REPO/data/refs/${GEN}.fa"
W="$REPO/work"; LOGS="$REPO/logs"
PLOIDY_LOG="$LOGS/ploidy_verification.txt"
mkdir -p "$W" "$LOGS"

ALIGNERS="bwa bowtie2 minimap2"
CALLERS="gatk freebayes bcftools"

# gatk is a Python launcher and needs its env bin on PATH (see NOTES 5.3).
# It must be invoked via `env`, NOT via a shell function: /usr/bin/time execs a
# real binary and cannot run a shell function, failing with
# "time: gatk: No such file or directory" — which looks like a missing install
# rather than a quoting problem.
GATK_RUN=(env "PATH=$C:$PATH" "$C/gatk")
secs_of() { awk '/ real /{printf "%.2f", $1}' "$1"; }

# ---------------------------------------------------------------------------
# R2 verification: read the GT field straight out of the VCF.
# ---------------------------------------------------------------------------
verify_ploidy() { # vcf aligner caller
  local vcf="$1" aln="$2" cal="$3"
  local gts n_hap n_dip total
  # Collect distinct GT values actually present.
  gts=$("$C/bcftools" query -f '[%GT]\n' "$vcf" 2>/dev/null | sort | uniq -c \
        | awk '{printf "%s(%s) ", $2, $1}')
  total=$("$C/bcftools" view -H "$vcf" 2>/dev/null | wc -l | tr -d ' ')
  n_dip=$("$C/bcftools" query -f '[%GT]\n' "$vcf" 2>/dev/null | grep -c '[/|]' || true)
  n_hap=$(( total - n_dip ))

  local verdict
  if [[ "$total" -eq 0 ]]; then
    verdict="NO CALLS — cannot verify"
  elif [[ "$n_dip" -eq 0 ]]; then
    verdict="PASS (haploid)"
  else
    verdict="FAIL (${n_dip}/${total} diploid genotypes)"
  fi

  printf '%-8s %-10s %-10s records=%-8s haploid=%-8s diploid=%-8s GTs: %-24s %s\n' \
    "$GEN" "$aln" "$cal" "$total" "$n_hap" "$n_dip" "$gts" "$verdict" \
    | tee -a "$PLOIDY_LOG"

  if [[ "$n_dip" -ne 0 ]]; then
    echo "FATAL: $cal emitted diploid genotypes on a haploid organism (R2)." >&2
    exit 1
  fi
}

echo "# ploidy verification — $(date -u '+%Y-%m-%d %H:%M:%S UTC') — $TAG" >> "$PLOIDY_LOG"

for ALN in $ALIGNERS; do
  BAM="$W/${TAG}.${ALN}.md.bam"
  [[ -s "$BAM" ]] || { echo "FATAL: missing $BAM" >&2; exit 1; }

  # ---- GATK HaplotypeCaller ------------------------------------------------
  OUT="$W/${TAG}.${ALN}.gatk.raw.vcf.gz"
  TF="$LOGS/timing/${TAG}.${ALN}.gatk.time"; mkdir -p "$LOGS/timing"
  echo "  [$ALN/gatk] calling..."
  /usr/bin/time -l "${GATK_RUN[@]}" HaplotypeCaller \
      -R "$REF" -I "$BAM" -O "$OUT" \
      --sample-ploidy 1 \
      --native-pair-hmm-threads "$THREADS" \
      > "$LOGS/gatk_${TAG}_${ALN}.log" 2> "$TF" || {
        echo "FATAL: GATK failed; see $TF" >&2
        tail -25 "$TF" >&2; exit 1; }
  "$C/bcftools" index -t -f "$OUT"
  verify_ploidy "$OUT" "$ALN" gatk
  printf '%s\t%s\t%s\t%s\n' "$TAG" "$ALN" gatk "$(secs_of "$TF")" >> "$LOGS/call_timing.tsv"

  # ---- FreeBayes -----------------------------------------------------------
  OUT="$W/${TAG}.${ALN}.freebayes.raw.vcf.gz"
  TF="$LOGS/timing/${TAG}.${ALN}.freebayes.time"
  echo "  [$ALN/freebayes] calling..."
  /usr/bin/time -l "$C/freebayes" -f "$REF" -p 1 "$BAM" \
      > "$W/.fb.$$.vcf" 2> "$TF" || {
        echo "FATAL: FreeBayes failed" >&2; tail -20 "$TF" >&2; exit 1; }
  "$C/bgzip" -c "$W/.fb.$$.vcf" > "$OUT"; rm -f "$W/.fb.$$.vcf"
  "$C/bcftools" index -t -f "$OUT"
  verify_ploidy "$OUT" "$ALN" freebayes
  printf '%s\t%s\t%s\t%s\n' "$TAG" "$ALN" freebayes "$(secs_of "$TF")" >> "$LOGS/call_timing.tsv"

  # ---- BCFtools mpileup/call ----------------------------------------------
  OUT="$W/${TAG}.${ALN}.bcftools.raw.vcf.gz"
  TF="$LOGS/timing/${TAG}.${ALN}.bcftools.time"
  echo "  [$ALN/bcftools] calling..."
  # -a AD,DP annotates depth so the shared hard filter can use DP (see below).
  /usr/bin/time -l sh -c \
      "'$C/bcftools' mpileup -f '$REF' -a AD,DP -Ou '$BAM' \
       | '$C/bcftools' call -mv --ploidy 1 -Oz -o '$OUT'" 2> "$TF" || {
        echo "FATAL: BCFtools failed" >&2; tail -20 "$TF" >&2; exit 1; }
  "$C/bcftools" index -t -f "$OUT"
  verify_ploidy "$OUT" "$ALN" bcftools
  printf '%s\t%s\t%s\t%s\n' "$TAG" "$ALN" bcftools "$(secs_of "$TF")" >> "$LOGS/call_timing.tsv"
done

# ---------------------------------------------------------------------------
# Hard-filtered callsets — IDENTICAL expression for every caller (R8).
#
# Only QUAL and DP are used, because they are the only fields all three callers
# emit with the same meaning. GATK Best Practices would filter on QD/FS/MQRankSum,
# but FreeBayes and BCFtools do not produce those, so using them would apply a
# different (and better-tuned) filter to GATK than to the others.
#
# CAVEAT, stated plainly: QUAL is NOT calibrated identically across these three
# tools. A QUAL of 20 does not mean the same thing to GATK as to FreeBayes. So a
# single threshold is *procedurally* identical but not *statistically* equivalent.
# This is exactly why the raw callsets retain QUAL and why ROC curves are produced
# in Phase 7a — the ROC sweeps the threshold and removes the arbitrariness.
# ---------------------------------------------------------------------------
FILTER_EXPR='QUAL>=20 && INFO/DP>=5'
echo "  applying shared hard filter: $FILTER_EXPR"
for ALN in $ALIGNERS; do
  for CAL in $CALLERS; do
    RAW="$W/${TAG}.${ALN}.${CAL}.raw.vcf.gz"
    FLT="$W/${TAG}.${ALN}.${CAL}.filt.vcf.gz"
    [[ -s "$RAW" ]] || continue
    "$C/bcftools" view -i "$FILTER_EXPR" -Oz -o "$FLT" "$RAW" 2>/dev/null
    "$C/bcftools" index -t -f "$FLT"
    printf '    %-9s %-10s raw=%-7s filtered=%s\n' "$ALN" "$CAL" \
      "$("$C/bcftools" view -H "$RAW" | wc -l | tr -d ' ')" \
      "$("$C/bcftools" view -H "$FLT" | wc -l | tr -d ' ')"
  done
done

echo "[$TAG] 9 raw + 9 filtered VCFs written to $W"
