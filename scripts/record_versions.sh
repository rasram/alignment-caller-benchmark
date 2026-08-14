#!/usr/bin/env bash
# Record the exact version of every tool the benchmark depends on (rule R9).
#
# Why this exists: aligner and caller behaviour changes between releases. A result
# table without versions cannot be reproduced or defended. This is run on day one
# and re-run whenever an environment is rebuilt.
#
# Usage:  bash scripts/record_versions.sh
# Writes: logs/versions_<date>.txt  and  logs/versions_latest.txt

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
RTG="$REPO/tools/rtg-tools-3.13/rtg"
SIMUG="$REPO/tools/simuG/simuG.pl"

DATE="$(date +%Y%m%d)"
OUT="$REPO/logs/versions_${DATE}.txt"
mkdir -p "$REPO/logs"

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Run a command inside a conda env and print its version.
#
# Tool version reporting is wildly inconsistent: some use --version, some -v, some
# only print usage with no arguments, several write to stderr, and several exit
# non-zero while still being perfectly healthy (art_illumina and bwa both do).
# So: merge stderr, ignore the exit code, and pull the first line matching an
# optional pattern rather than trusting line 1.
#
# `grep -a` is required — bowtie2 --version emits its own binary path, and without
# -a GNU/BSD grep decides the stream is binary and prints "Binary file matches"
# instead of the version. That is exactly the kind of thing that silently turns a
# provenance record into garbage.
ver() {
  local env="$1" label="$2" pat="$3"; shift 3
  local out
  if [[ -n "$pat" ]]; then
    out="$(conda run -n "$env" "$@" 2>&1 | grep -a -m1 -E "$pat")"
  else
    out="$(conda run -n "$env" "$@" 2>&1 | grep -a -v '^[[:space:]]*$' | head -1)"
  fi
  out="$(echo "$out" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  if [[ -z "$out" ]]; then
    printf '%-14s NOT FOUND / no version string\n' "$label"
  else
    printf '%-14s %s\n' "$label" "$out"
  fi
}

{
  echo "# Tool versions — recorded $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "# Host: $(uname -s) $(uname -r) $(uname -m)"
  echo "# macOS: $(sw_vers -productVersion 2>/dev/null || echo n/a)"
  echo "# conda: $("$CONDA_BASE/bin/conda" --version 2>&1)"
  echo "# NOTE: all tools are native $(uname -m) builds — no Rosetta emulation,"
  echo "#       so Phase 5 runtime comparisons are on equal footing."
  echo

  echo "## env: align"
  ver align "bwa"        '^Version'      bwa
  ver align "bowtie2"    'version'       bowtie2 --version
  ver align "minimap2"   ''              minimap2 --version
  ver align "samtools"   '^samtools'     samtools --version
  echo

  echo "## env: callers"
  ver callers "gatk"      'Genome Analysis Toolkit'  gatk --version
  ver callers "gatk-java" 'openjdk version'          java -version
  ver callers "freebayes" 'version'      freebayes --version
  ver callers "bcftools"  '^bcftools'    bcftools --version
  echo

  echo "## env: sim"
  # art_illumina has no --version flag; it prints a banner on bare invocation.
  ver sim "art_illumina" 'Version'       art_illumina
  ver sim "perl"         'This is perl'  perl --version
  echo

  echo "## env: qc"
  ver qc "fastqc"  ''  fastqc --version
  ver qc "multiqc" ''  multiqc --version
  echo

  echo "## env: ml"
  ver ml "python"     ''  python --version
  ver ml "snakemake"  ''  snakemake --version
  ver ml "pandas"     ''  python -c "import pandas; print(pandas.__version__)"
  ver ml "sklearn"    ''  python -c "import sklearn; print(sklearn.__version__)"
  ver ml "matplotlib" ''  python -c "import matplotlib; print(matplotlib.__version__)"
  echo

  echo "## outside conda"
  printf '%-14s %s\n' "java" "$(java -version 2>&1 | head -1)"
  if [[ -x "$RTG" ]]; then
    printf '%-14s %s\n' "rtg" "$("$RTG" version 2>&1 | head -1)"
    printf '%-14s %s\n' "rtg-core" "$("$RTG" version 2>&1 | sed -n '2p')"
  else
    printf '%-14s NOT FOUND at %s\n' "rtg" "$RTG"
  fi
  if [[ -f "$SIMUG" ]]; then
    printf '%-14s %s\n' "simuG.pl" "present ($(cd "$REPO/tools/simuG" && git rev-parse --short HEAD 2>/dev/null || echo 'unknown rev'))"
  else
    printf '%-14s NOT FOUND at %s\n' "simuG.pl" "$SIMUG"
  fi
} | tee "$OUT"

cp "$OUT" "$REPO/logs/versions_latest.txt"
echo
echo "Wrote $OUT"
