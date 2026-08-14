#!/usr/bin/env bash
# PHASE 1 — download and index the two reference genomes.
#
# Genomes:
#   phiX  = NC_001422.1              (Escherichia virus phiX174),      5,386 bp
#   ecoli = NC_000913.3 / GCF_000005845.2 (E. coli K-12 MG1655),   4,641,652 bp
#
# RefSeq (GCF_) not GenBank (GCA_): the brief requires it, and RefSeq is the
# curated copy that downstream annotation resources key off.
#
# CONTIG RENAMING (rule R5). NCBI FASTA headers look like
#     >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
# Two problems: (a) the description after the first space is carried into some
# tools' output and dropped by others, and (b) accession-with-version strings are
# easy to mistype inconsistently across a dozen commands. Every downstream file
# — truth VCF, BED, every call set, the RTG SDF — must agree on the contig name
# exactly, or vcfeval silently scores zero overlap. So the name is normalised
# ONCE, here, at download time, to a short token: `phiX` and `ecoli`.
#
# Usage: bash scripts/fetch_references.sh

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
REFS="$REPO/data/refs"
mkdir -p "$REFS" "$REPO/logs"

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# macOS ships bash 3.2, which has no associative arrays (`declare -A`). Rather than
# add a bash-4 dependency just for two lookups, these are plain case statements.
NAMES="phiX ecoli"

acc_of() {
  case "$1" in
    phiX)  echo "NC_001422.1" ;;
    ecoli) echo "NC_000913.3" ;;
    *)     echo "UNKNOWN"     ;;
  esac
}

explen_of() {
  case "$1" in
    phiX)  echo 5386    ;;
    ecoli) echo 4641652 ;;
    *)     echo 0       ;;
  esac
}

fetch() {
  local name="$1"; local acc; acc="$(acc_of "$1")"
  local fa="$REFS/${name}.fa"

  if [[ -s "$fa" ]]; then
    echo "[$name] already present, skipping download"
    return
  fi

  echo "[$name] downloading $acc from NCBI E-utilities"
  # efetch returns the FASTA for a nucleotide accession. Retried because NCBI
  # rate-limits anonymous requests and an empty body would otherwise be written
  # out as a valid-looking but truncated FASTA.
  local url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
  local tmp="$fa.tmp"
  local ok=0
  for attempt in 1 2 3 4 5; do
    if curl -sSL --max-time 300 -o "$tmp" "$url" && [[ -s "$tmp" ]] && head -1 "$tmp" | grep -q '^>'; then
      ok=1; break
    fi
    echo "[$name] attempt $attempt failed, retrying"
    sleep 5
  done
  [[ $ok -eq 1 ]] || { echo "[$name] FATAL: could not download $acc" >&2; exit 1; }

  # Record the original header before we destroy it, so provenance is auditable.
  head -1 "$tmp" > "$REFS/${name}.original_header.txt"

  # Rewrite the header to the short name (R5) and hard-wrap sequence at 60 cols.
  awk -v nm="$name" 'BEGIN{OFS=""}
    /^>/ { if (!seen) { print ">", nm; seen=1 } else { print "ERROR: multi-contig input" > "/dev/stderr"; exit 1 } ; next }
    { print }
  ' "$tmp" > "$fa"
  rm -f "$tmp"
  echo "[$name] wrote $fa"
}

index() {
  local name="$1"
  local fa="$REFS/${name}.fa"

  echo "[$name] samtools faidx"
  conda run -n align samtools faidx "$fa"

  echo "[$name] bwa index"
  conda run -n align bwa index "$fa" 2> "$REPO/logs/bwa_index_${name}.log"

  echo "[$name] bowtie2-build"
  conda run -n align bowtie2-build --threads 4 "$fa" "$REFS/${name}" \
    > "$REPO/logs/bowtie2_build_${name}.log" 2>&1

  echo "[$name] gatk CreateSequenceDictionary"
  rm -f "$REFS/${name}.dict"
  conda run -n callers gatk CreateSequenceDictionary -R "$fa" \
    > "$REPO/logs/gatk_dict_${name}.log" 2>&1

  # RTG's SDF is its own indexed reference format; vcfeval needs it, and it must
  # be built from the SAME renamed FASTA or contig names will not match.
  echo "[$name] rtg format -> SDF"
  rm -rf "$REFS/${name}.sdf"
  "$RTG" format -o "$REFS/${name}.sdf" "$fa" > "$REPO/logs/rtg_format_${name}.log" 2>&1
}

verify() {
  local name="$1"; local exp; exp="$(explen_of "$1")"
  local fai="$REFS/${name}.fa.fai"
  local contig len
  contig=$(cut -f1 "$fai")
  len=$(cut -f2 "$fai")

  printf '%-8s contig=%-8s length=%-10s expected=%-10s ' "$name" "$contig" "$len" "$exp"
  if [[ "$len" == "$exp" ]]; then
    echo "OK"
  else
    echo "MISMATCH"
    echo "FATAL: $name length $len != expected $exp — wrong assembly version." >&2
    exit 1
  fi
}

for n in $NAMES; do fetch  "$n"; done
for n in $NAMES; do index  "$n"; done
echo
echo "=== GATE 1 verification ==="
for n in $NAMES; do verify "$n"; done
