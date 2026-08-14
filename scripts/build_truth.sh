#!/usr/bin/env bash
# PHASE 2 — turn simuG's raw output into a scoreable truth set.
#
# simuG writes two separate VCFs (SNP and INDEL) that are NOT ready for vcfeval:
#   * no ##contig header lines            -> bcftools/vcfeval reject or mis-sort
#   * no FORMAT column and no sample       -> vcfeval has no genotype to compare
#   * not normalised (rule R3)             -> indels mismatch on representation
#   * QUAL is '.'                          -> fine for truth, it is never scored
#
# This script produces, per genome:
#   data/truth/<g>.truth.vcf.gz        normalised, GT=1, indexed   <- the truth set
#   data/truth/<g>.truth.nogt.vcf.gz   normalised, no sample col   <- for the
#                                       --squash-ploidy comparison in NOTES 2.3
#   data/truth/<g>.confident.bed       confident regions (whole genome here)
#
# Usage: bash scripts/build_truth.sh <genome>      e.g. phiX | ecoli

set -euo pipefail

GEN="${1:?usage: build_truth.sh <genome>}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tools.sh"
REF="$REPO/data/refs/${GEN}.fa"
FAI="$REF.fai"
T="$REPO/data/truth"

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# R3/§0.10: ALL normalisation in this project uses this one bcftools binary.
BCFTOOLS="$(resolve_tool bcftools callers)"
bcf() { "$BCFTOOLS" "$@"; }

SNPVCF="$T/${GEN}.refseq2simseq.SNP.vcf"
INDVCF="$T/${GEN}.refseq2simseq.INDEL.vcf"
[[ -s "$SNPVCF" ]] || { echo "missing $SNPVCF" >&2; exit 1; }
[[ -s "$INDVCF" ]] || { echo "missing $INDVCF" >&2; exit 1; }

RAW_SNP=$(grep -vc '^#' "$SNPVCF" || true)
RAW_IND=$(grep -vc '^#' "$INDVCF" || true)
echo "[$GEN] simuG raw records: SNP=$RAW_SNP INDEL=$RAW_IND"

work="$T/.tmp_${GEN}"
rm -rf "$work"; mkdir -p "$work"

# ---------------------------------------------------------------------------
# 1. Build a proper VCF header.
#    ##contig lines come from the .fai, so contig naming is inherited from the
#    reference rather than retyped (rule R5).
# ---------------------------------------------------------------------------
{
  echo '##fileformat=VCFv4.2'
  echo "##source=simuG.pl (via build_truth.sh)"
  echo "##reference=file://${REF}"
  awk 'BEGIN{OFS=""}{print "##contig=<ID=", $1, ",length=", $2, ">"}' "$FAI"
  # Carry over simuG's INFO definitions so its annotations stay valid.
  grep '^##INFO' "$SNPVCF"
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
} > "$work/header_common.txt"

# ---------------------------------------------------------------------------
# 2. Concatenate SNP + INDEL bodies.
# ---------------------------------------------------------------------------
grep -v '^#' "$SNPVCF"  > "$work/body.txt"
grep -v '^#' "$INDVCF" >> "$work/body.txt"

# ---------------------------------------------------------------------------
# 3a. Variant WITH genotype: append FORMAT=GT and a haploid sample "sim" (GT=1).
#     GT=1 (single allele index, no slash) is the correct haploid encoding.
#     "1/1" would be diploid-homozygous and "0/1" diploid-het — both wrong for
#     these organisms and both a source of vcfeval genotype mismatches.
# ---------------------------------------------------------------------------
{
  cat "$work/header_common.txt"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsim\n'
  awk 'BEGIN{OFS="\t"}{print $0, "GT", "1"}' "$work/body.txt"
} > "$work/gt.vcf"

# ---------------------------------------------------------------------------
# 3b. Variant WITHOUT genotype: the 8 fixed columns only (as simuG emits).
#     Used only to test whether vcfeval --squash-ploidy accepts a sample-less
#     baseline. See NOTES 2.3.
# ---------------------------------------------------------------------------
{
  cat "$work/header_common.txt"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  cat "$work/body.txt"
} > "$work/nogt.vcf"

# ---------------------------------------------------------------------------
# 4. Sort -> normalise -> compress -> index.
#
#    `bcftools norm -f ref -m -any` does two things (rule R3):
#      -f ref     left-align and trim indels against the reference. The same indel
#                 in a repeat can be written at several positions; left-alignment
#                 forces one canonical choice.
#      -m -any    split any multi-allelic record into one record per ALT allele.
#      --atomize  decompose MNVs/complex records into consecutive atomic SNVs and
#                 indels. Added after Phase 6 found FreeBayes merges nearby
#                 variants into single records (215 AAA>TAT covering two SNPs;
#                 1699 CCGTCCTT>GCGTCTT covering a SNP and an indel). Without
#                 atomising, `bcftools view -v snps` drops those records and the
#                 variants inside them, understating FreeBayes recall as a pure
#                 representation artefact. See NOTES 6.6.
#    The truth set and every call set get IDENTICAL treatment, so representation
#    differences cannot masquerade as FP/FN. The truth contains no MNVs, so
#    atomising it changes nothing — but R3 requires identical TREATMENT, not
#    identical outcome, so it is atomised too.
# ---------------------------------------------------------------------------
norm_one() {
  local in="$1" out="$2" label="$3"
  bcf sort "$in" -Oz -o "$work/${label}.sorted.vcf.gz" 2>"$work/${label}.sort.log"
  bcf index -t -f "$work/${label}.sorted.vcf.gz"
  bcf norm -f "$REF" -m -any --atomize \
      -Oz -o "$out" "$work/${label}.sorted.vcf.gz" 2>"$work/${label}.norm.log"
  bcf index -t -f "$out"
  echo "  norm($label): $(grep -E 'total|realigned|split|changed' "$work/${label}.norm.log" | tr '\n' ' ')"
}

norm_one "$work/gt.vcf"   "$T/${GEN}.truth.vcf.gz"      "gt"
norm_one "$work/nogt.vcf" "$T/${GEN}.truth.nogt.vcf.gz" "nogt"

# ---------------------------------------------------------------------------
# 5. Confident-regions BED — the whole genome, because this data is simulated.
#    BED is 0-based half-open, hence "0, length".
# ---------------------------------------------------------------------------
awk 'BEGIN{OFS="\t"}{print $1, 0, $2}' "$FAI" > "$T/${GEN}.confident.bed"

# ---------------------------------------------------------------------------
# 6. Report final counts by type.
# ---------------------------------------------------------------------------
FIN_SNP=$(bcf view -H -v snps   "$T/${GEN}.truth.vcf.gz" | wc -l | tr -d ' ')
FIN_IND=$(bcf view -H -v indels "$T/${GEN}.truth.vcf.gz" | wc -l | tr -d ' ')
FIN_ALL=$(bcf view -H           "$T/${GEN}.truth.vcf.gz" | wc -l | tr -d ' ')
echo "[$GEN] truth set: SNP=$FIN_SNP INDEL=$FIN_IND TOTAL=$FIN_ALL"
echo "[$GEN] confident BED: $(cat "$T/${GEN}.confident.bed")"

rm -rf "$work"
