#!/usr/bin/env bash
# PHASE 4 — read quality control.
#
# FastQC inspects the raw reads and reports per-base quality, GC content,
# adapter content, duplication and so on. MultiQC then merges the per-file
# FastQC reports into one HTML page.
#
# NOTE (R8): this is diagnostic ONLY. Nothing here modifies the reads. Real
# workflows often trim low-quality 3' ends at this point; we deliberately do not,
# because every pipeline must see byte-identical input or differences between
# aligners could be caused by the trimmer rather than by the aligners.
#
# Usage: bash scripts/run_qc.sh [fastq ...]
#        with no arguments, runs on all baseline FASTQs in work/

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
FASTQC="$CONDA_BASE/envs/qc/bin/fastqc"
MULTIQC="$CONDA_BASE/envs/qc/bin/multiqc"
QC="$REPO/results/qc"
mkdir -p "$QC" "$REPO/logs"

if [[ $# -gt 0 ]]; then
  FQ=("$@")
else
  FQ=()
  while IFS= read -r f; do FQ+=("$f"); done < <(ls "$REPO"/work/*_cov30_len150_err0_seed1_[12].fq 2>/dev/null)
fi

[[ ${#FQ[@]} -gt 0 ]] || { echo "FATAL: no FASTQ files found" >&2; exit 1; }

echo "Running FastQC on ${#FQ[@]} file(s)..."
# --threads: one worker per file, capped at 4 to match the thread budget used
# for the aligners (R8 — keep resource use comparable across phases).
"$FASTQC" --outdir "$QC" --threads 4 --quiet "${FQ[@]}" \
  > "$REPO/logs/fastqc.log" 2>&1

echo "Aggregating with MultiQC..."
"$MULTIQC" "$QC" --outdir "$QC" --force --quiet \
  > "$REPO/logs/multiqc.log" 2>&1

echo "FastQC + MultiQC written to $QC"
ls "$QC"/*.html 2>/dev/null | sed 's|.*/|  |'
