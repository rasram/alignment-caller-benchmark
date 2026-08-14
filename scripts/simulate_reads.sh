#!/usr/bin/env bash
# PHASE 3 — simulate sequencing reads with ART.
#
# RULE R1, the one that invalidates everything if reversed:
#   reads are simulated FROM the MUTATED genome (data/truth/<g>.simseq.genome.fa)
#   and will be aligned TO the ORIGINAL reference (data/refs/<g>.fa) in Phase 5.
# This script therefore takes the mutated genome as input, and asserts it.
#
# Baseline condition: 30x coverage, 150 bp paired-end, HS25 profile, qs shift 0.
#
# Outputs, per condition, into work/:
#   <tag>1.fq  <tag>2.fq   the paired reads (R8: never trimmed, all pipelines
#                          receive byte-identical input)
#   <tag>.sam              ART's TRUTH alignment — where each read really came
#                          from. Retained deliberately: Phase 5 uses it to compute
#                          placement accuracy, a pure aligner metric that is
#                          independent of any variant caller.
#
# Usage: bash scripts/simulate_reads.sh <genome> [cov] [len] [qs] [seed]
#        bash scripts/simulate_reads.sh phiX 30 150 0 1

set -euo pipefail

GEN="${1:?usage: simulate_reads.sh <genome> [cov] [len] [qs] [seed]}"
COV="${2:-30}"
LEN="${3:-150}"
QS="${4:-0}"
SEED="${5:-1}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tools.sh"
ART="$(resolve_tool art_illumina sim)"
SAMTOOLS="$(resolve_tool samtools align)"

MUT="$REPO/data/truth/${GEN}.simseq.genome.fa"
REF="$REPO/data/refs/${GEN}.fa"
W="$REPO/work"; mkdir -p "$W" "$REPO/logs"

TAG="${GEN}_cov${COV}_len${LEN}_err${QS}_seed${SEED}"
OUT="$W/${TAG}_"

# --- R1 guard: refuse to run against the original reference -----------------
[[ -s "$MUT" ]] || { echo "FATAL: mutated genome missing: $MUT" >&2; exit 1; }
if cmp -s "$MUT" "$REF"; then
  echo "FATAL: mutated genome is identical to the reference — R1 violated." >&2
  exit 1
fi

echo "[$TAG] simulating from MUTATED genome: $MUT"

# -ss HS25   HiSeq 2500 error profile
# -p         paired-end
# -m 350     mean DNA fragment length; -s 50 its standard deviation
# -rs SEED   random seed (R9: recorded in the filename AND the results table)
# -sam       emit the truth SAM (see header comment)
# -na        skip the .aln files; the SAM carries the same information
# -qs/-qs2   quality-score shift for read 1 / read 2. 0 = the profile's own
#            quality distribution, untouched.
"$ART" -ss HS25 -sam -na \
  -i "$MUT" -p -l "$LEN" -f "$COV" -m 350 -s 50 \
  -qs "$QS" -qs2 "$QS" \
  -rs "$SEED" -o "$OUT" > "$REPO/logs/art_${TAG}.log" 2>&1

for f in "${OUT}1.fq" "${OUT}2.fq" "${OUT}.sam"; do
  [[ -s "$f" ]] || { echo "FATAL: ART did not produce $f" >&2; exit 1; }
done

# --- coverage verification --------------------------------------------------
# Coverage is measured against the MUTATED genome, because that is the template
# the reads were drawn from. Using the reference length instead would be wrong by
# the net indel balance (small here, but wrong on principle).
MUTLEN=$(awk '!/^>/{n+=length($0)} END{print n+0}' "$MUT")
R1=$(( $(wc -l < "${OUT}1.fq") / 4 ))
R2=$(( $(wc -l < "${OUT}2.fq") / 4 ))
TOT=$(( R1 + R2 ))
ACT=$(python3 -c "print(f'{($TOT * $LEN) / $MUTLEN:.2f}')")

printf '%-38s pairs=%-9s reads=%-10s actual_cov=%sx (requested %sx)\n' \
  "$TAG" "$R1" "$TOT" "$ACT" "$COV"

{
  echo "tag=$TAG genome=$GEN mutated_len=$MUTLEN read_len=$LEN seed=$SEED qs=$QS"
  echo "pairs=$R1 reads=$TOT requested_cov=${COV}x actual_cov=${ACT}x"
} >> "$REPO/logs/coverage_verification.txt"
