#!/usr/bin/env bash
# PHASE 2 verification — prove the truth set is what we think it is.
#
# Three independent checks, because each catches a different class of error:
#
#   1. COUNTS      — do we have the number of variants we asked simuG for?
#                    Catches silent truncation / a failed concat.
#   2. REF BASES   — does every VCF REF field actually match the reference FASTA
#                    at that coordinate? Catches off-by-one and coordinate-system
#                    errors (VCF is 1-based inclusive; BED is 0-based half-open;
#                    mixing them is the classic bug here).
#   3. LENGTH BOOK — does (mutated genome length - reference length) equal the net
#                    indel balance? This is an end-to-end check that the truth VCF
#                    describes the SAME edits that were actually applied to the
#                    mutated genome the reads will be simulated from.
#
# Check 3 is the strongest of the three: it ties the truth VCF to the FASTA that
# Phase 3 reads are generated from. If simuG's VCF and its mutated genome ever
# disagreed, every downstream number would be wrong, and checks 1 and 2 would
# both still pass.
#
# Usage: bash scripts/verify_truth.sh <genome> [n_spotcheck]

set -euo pipefail

GEN="${1:?usage: verify_truth.sh <genome> [n_spotcheck]}"
NSPOT="${2:-3}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
REF="$REPO/data/refs/${GEN}.fa"
T="$REPO/data/truth"
TRUTH="$T/${GEN}.truth.vcf.gz"
MUT="$T/${GEN}.simseq.genome.fa"
LOG="$REPO/logs/truth_verification_${GEN}.txt"

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
BCFTOOLS="$CONDA_BASE/envs/callers/bin/bcftools"
SAMTOOLS="$CONDA_BASE/envs/align/bin/samtools"
bcf(){ "$BCFTOOLS" "$@"; }
sam(){ "$SAMTOOLS" "$@"; }

{
echo "# Truth-set verification — $GEN — $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo

# ---- 1. counts -------------------------------------------------------------
NSNP=$(bcf view -H -v snps   "$TRUTH" | wc -l | tr -d ' ')
NIND=$(bcf view -H -v indels "$TRUTH" | wc -l | tr -d ' ')
echo "## 1. Variant counts"
echo "SNPs   : $NSNP"
echo "INDELs : $NIND"
echo "TOTAL  : $((NSNP + NIND))"
GLEN=$(cut -f2 "$REF.fai")
echo "Genome : $GLEN bp  -> 1 variant per $(( GLEN / (NSNP + NIND) )) bp"
echo

# ---- 2. indel size distribution -------------------------------------------
echo "## 2. Indel length distribution (ALT length - REF length; +=insertion)"
bcf view -H -v indels "$TRUTH" \
  | awk '{d=length($5)-length($4); print d}' \
  | sort -n | uniq -c \
  | awk '{printf "  size %+4d : %5d\n", $2, $1}'
echo
echo "  insertions: $(bcf view -H -v indels "$TRUTH" | awk 'length($5)>length($4)' | wc -l | tr -d ' ')"
echo "  deletions : $(bcf view -H -v indels "$TRUTH" | awk 'length($5)<length($4)' | wc -l | tr -d ' ')"
echo

# ---- 3. Ti/Tv ratio --------------------------------------------------------
echo "## 3. SNP transition/transversion ratio"
bcf view -H -v snps "$TRUTH" | awk '
  function ti(r,a){ return (r=="A"&&a=="G")||(r=="G"&&a=="A")||(r=="C"&&a=="T")||(r=="T"&&a=="C") }
  { if (ti($4,$5)) t++; else v++ }
  END { printf "  transitions=%d transversions=%d Ti/Tv=%.3f\n", t, v, (v?t/v:0) }'
echo

# ---- 4. manual REF spot-check ---------------------------------------------
echo "## 4. Manual spot-check: VCF REF vs samtools faidx on the reference"
printf "  %-6s %-10s %-14s %-14s %s\n" "TYPE" "POS" "VCF_REF" "FASTA_REF" "MATCH"
SPOT_FAIL=0
for vt in snps indels; do
  # NOTE: `bcftools ... | head -N` makes bcftools die on SIGPIPE, which under
  # `set -o pipefail` aborts this whole script — silently skipping the indel
  # checks. Capture to a variable with `|| true` instead of piping into head.
  recs="$(bcf view -H -v "$vt" "$TRUTH" 2>/dev/null | head -"$NSPOT" || true)"
  while IFS=$'\t' read -r c p i r a rest; do
    [ -n "${p:-}" ] || continue
    end=$(( p + ${#r} - 1 ))
    fa=$(sam faidx "$REF" "${c}:${p}-${end}" 2>/dev/null | grep -v '^>' | tr -d '\n')
    if [ "$r" = "$fa" ]; then m="OK"; else m="MISMATCH"; SPOT_FAIL=1; fi
    printf "  %-6s %-10s %-14s %-14s %s\n" "$vt" "$p" "$r" "$fa" "$m"
  done <<< "$recs"
done
if [ "$SPOT_FAIL" -eq 0 ]; then
  echo "  RESULT: OK — every checked REF matches the reference FASTA."
else
  echo "  RESULT: MISMATCH — coordinate or reference error. STOP."
fi
echo

# ---- 5. length bookkeeping (the strong check) ------------------------------
echo "## 5. Length bookkeeping: truth VCF vs the actual mutated genome"
REFLEN=$(grep -v '^>' "$REF" | tr -d '\n' | wc -c | tr -d ' ')
MUTLEN=$(grep -v '^>' "$MUT" | tr -d '\n' | wc -c | tr -d ' ')
NETVCF=$(bcf view -H -v indels "$TRUTH" | awk '{s += length($5)-length($4)} END{print s+0}')
DELTA=$(( MUTLEN - REFLEN ))
echo "  reference length      : $REFLEN"
echo "  mutated genome length : $MUTLEN"
echo "  observed difference   : $DELTA"
echo "  net indel sum in VCF  : $NETVCF"
if [ "$DELTA" -eq "$NETVCF" ]; then
  echo "  RESULT: OK — the truth VCF accounts for every base of length change."
else
  echo "  RESULT: MISMATCH — truth VCF does NOT describe the mutated genome. STOP."
  FAILED=1
fi
echo

echo "## Contig name consistency (R5)"
printf "  ref.fai : %s\n" "$(cut -f1 "$REF.fai" | tr '\n' ' ')"
printf "  truth   : %s\n" "$(bcf view -H "$TRUTH" | cut -f1 | sort -u | tr '\n' ' ')"
printf "  bed     : %s\n" "$(cut -f1 "$T/${GEN}.confident.bed" | tr '\n' ' ')"
printf "  mutated : %s\n" "$(grep '^>' "$MUT" | tr -d '>' | tr '\n' ' ')"
} | tee "$LOG"

echo
echo "Wrote $LOG"
