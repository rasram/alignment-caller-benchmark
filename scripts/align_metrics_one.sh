#!/usr/bin/env bash
# Alignment metrics for ONE (tag, aligner) — Snakemake rule `align_metrics`.
#
# Same definitions as align_metrics.sh (see its header for why each metric
# exists), restructured for the workflow:
#   * one unit in, one row out, written to its own file. The old script looped
#     over aligners and APPENDED to a shared TSV, which is unsafe under parallel
#     execution — concurrent appends interleave and corrupt rows.
#   * no timing here. Runtime and peak RSS come from Snakemake's `benchmark:`
#     directive on the aligner rule itself (scripts/aggregate_runtime.py).
#
# Usage: align_metrics_one.sh <tag> <aligner> <md.bam> <markdup.metrics> > out.tsv
set -euo pipefail
TAG="$1"; ALN="$2"; BAM="$3"; MET="$4"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$REPO/scripts/lib/tools.sh"
S="$(resolve_tool samtools align)"
GEN="${TAG%%_*}"
REFLEN=$(cut -f2 "$REPO/data/refs/${GEN}.fa.fai")

# FAIRNESS (R8): PRIMARY counts, never flagstat "in total" — that line includes
# supplementary alignments, which BWA-MEM emits and the other two do not, giving
# BWA a larger denominator for the same input reads (NOTES 5.5).
#
# Here-strings, NOT `echo "$FS" | awk '...exit'`. awk exits as soon as it has its
# line and closes the pipe; if echo is still writing it dies of SIGPIPE, and under
# pipefail + set -e the whole script aborts with status 141. A race: it struck once
# in ~90 jobs during the sweep, never in the pilot. No pipe, no SIGPIPE.
FS="$("$S" flagstat "$BAM")"
TOTAL=$(awk '/ primary$/{print $1; exit}'          <<< "$FS")
MAPPED=$(awk '/primary mapped \(/{print $1; exit}' <<< "$FS")
PPAIR=$(awk '/properly paired/{print $1; exit}'    <<< "$FS")
SUPPL=$(awk '/supplementary/{print $1; exit}'      <<< "$FS")

# printf needs the trailing \n, or `read` returns non-zero under set -e.
read -r MEANQ MAPQ0F < <("$S" view -F 0x900 "$BAM" \
  | awk '{s+=$5; n++; if($5==0) z++} END{if(n) printf "%.3f %.6f\n", s/n, z/n; else print "NA NA"}')

MEANDEP=$("$S" depth -a "$BAM" | awk -v L="$REFLEN" '{s+=$3} END{printf "%.3f", (L? s/L : 0)}')
DUP=$(awk -F'\t' '$1=="LIBRARY"{h=1;next} h&&NF>7{print $9; exit}' "$MET" 2>/dev/null || true)

rate() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b>0) printf "%.6f", a/b; else print "NA" }'; }

printf 'tag\taligner\tprimary_reads\tmapped\tmapping_rate\tproperly_paired\tproperly_paired_rate\tmean_mapq\tmapq0_frac\tmean_depth\tsupplementary\tdup_rate\n'
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$TAG" "$ALN" "$TOTAL" "$MAPPED" "$(rate "$MAPPED" "$TOTAL")" "$PPAIR" "$(rate "$PPAIR" "$TOTAL")" \
  "$MEANQ" "$MAPQ0F" "$MEANDEP" "$SUPPL" "${DUP:-NA}"
