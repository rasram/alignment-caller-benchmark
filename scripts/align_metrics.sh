#!/usr/bin/env bash
# PHASE 5 — collect per-BAM alignment metrics into results/align_metrics.tsv.
#
# Metrics and why each is here:
#   mapping_rate       fraction of reads the aligner placed anywhere at all.
#   properly_paired    fraction of pairs placed in the expected orientation AND
#                      at roughly the expected distance apart. A read can be
#                      "mapped" while its pair relationship is nonsense, so this
#                      is a stricter and more informative measure than mapping rate.
#   mean_mapq          MAPQ is the aligner's own confidence that a read is in the
#                      right place, on a Phred scale. Callers WEIGHT reads by it,
#                      so a systematically lower mean MAPQ propagates into
#                      variant calls even when placement is identical.
#   mapq0_frac         reads the aligner explicitly flags as multi-mapping
#                      (MAPQ 0 = "could be several places"). Most callers ignore
#                      these entirely, so they are effectively lost coverage.
#   mean_depth         average number of reads covering each reference base.
#   dup_rate           MarkDuplicates output; expected ~0 for simulated reads.
#
# Usage: bash scripts/align_metrics.sh <tag> [tag ...]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
S="$CONDA_BASE/envs/align/bin/samtools"
W="$REPO/work"; LOGS="$REPO/logs"
OUT="$REPO/results/align_metrics.tsv"
mkdir -p "$REPO/results"

printf 'tag\taligner\tprimary_reads\tmapped\tmapping_rate\tproperly_paired\tproperly_paired_rate\tmean_mapq\tmapq0_frac\tmean_depth\tsupplementary\tdup_rate\talign_seconds\tpeak_rss_mb\n' > "$OUT"

for TAG in "$@"; do
  GEN="${TAG%%_*}"
  REFLEN=$(cut -f2 "$REPO/data/refs/${GEN}.fa.fai")

  for A in bwa bowtie2 minimap2; do
    BAM="$W/${TAG}.${A}.md.bam"
    [[ -s "$BAM" ]] || { echo "skip: no $BAM" >&2; continue; }

    # flagstat gives counts. Use PRIMARY counts, not "in total".
    #
    # FAIRNESS (R8): "in total" includes supplementary alignments, which are the
    # extra records an aligner emits when it splits a chimeric read across two
    # locations. BWA-MEM emits them (36 here); bowtie2 and minimap2 -ax sr do
    # not. Using "in total" therefore gives BWA a LARGER denominator than the
    # others for the same 928,350 input reads, so the three mapping rates would
    # not be computed over the same quantity. "primary" is exactly one record
    # per input read for every aligner, which is the comparable basis.
    FS="$($S flagstat "$BAM")"
    TOTAL=$(echo "$FS"  | awk '/ primary$/{print $1; exit}')
    MAPPED=$(echo "$FS" | awk '/primary mapped \(/{print $1; exit}')
    PPAIR=$(echo "$FS"  | awk '/properly paired/{print $1; exit}')
    SUPPL=$(echo "$FS"  | awk '/supplementary/{print $1; exit}')

    # MAPQ mean and MAPQ0 fraction over primary alignments only.
    # NOTE: awk's printf must end with \n. Without a trailing newline `read`
    # returns non-zero, and under `set -e` that silently kills the whole script.
    read -r MEANQ MAPQ0F < <($S view -F 0x900 "$BAM" \
      | awk '{s+=$5; n++; if($5==0) z++} END{if(n) printf "%.3f %.6f\n", s/n, z/n; else print "0 0"}')

    # Mean depth across the reference (samtools depth -a includes zero-depth bases).
    MEANDEP=$($S depth -a "$BAM" | awk -v L="$REFLEN" '{s+=$3} END{if(L) printf "%.3f", s/L; else print 0}')

    MET="$LOGS/${TAG}.${A}.md.metrics"
    DUP=$(awk -F'\t' '$1=="LIBRARY"{h=1;next} h&&NF>7{print $9; exit}' "$MET" 2>/dev/null || echo "NA")
    [[ -n "$DUP" ]] || DUP="NA"

    SECS=NA; RSS=NA
    if [[ -f "$LOGS/align_timing.tsv" ]]; then
      read -r SECS RSS < <(awk -v t="$TAG" -v a="$A" \
        '$1==t && $2==a {print $3, $4; found=1; exit} END{if(!found) print "NA NA"}' \
        "$LOGS/align_timing.tsv") || true
    fi

    MRATE=$(python3 -c "print(f'{$MAPPED/$TOTAL:.6f}')" 2>/dev/null || echo NA)
    PRATE=$(python3 -c "print(f'{$PPAIR/$TOTAL:.6f}')" 2>/dev/null || echo NA)

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$TAG" "$A" "$TOTAL" "$MAPPED" "$MRATE" "$PPAIR" "$PRATE" \
      "$MEANQ" "$MAPQ0F" "$MEANDEP" "$SUPPL" "$DUP" "$SECS" "$RSS" >> "$OUT"
    echo "  collected $TAG / $A"
  done
done

echo
column -t -s$'\t' "$OUT"
echo
echo "Wrote $OUT"
