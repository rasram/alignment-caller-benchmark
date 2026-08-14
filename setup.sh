#!/usr/bin/env bash
# One-time setup: install the tools that are NOT conda packages, create the five
# conda environments, and record every version.
#
# Safe to re-run — every step is skipped if already done.
#
# Usage: bash setup.sh

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
RTG_VERSION="3.13"

echo "=== Setup: $REPO ==="

# --- 0. prerequisites --------------------------------------------------------
if [[ ! -x "$CONDA_BASE/bin/conda" ]]; then
  cat >&2 << EOF
FATAL: no conda at $CONDA_BASE

Install Miniforge first:
  curl -L -o miniforge.sh \\
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-\$(uname)-\$(uname -m).sh"
  bash miniforge.sh -b -p "\$HOME/miniforge3"

Or set CONDA_BASE to an existing installation.
EOF
  exit 1
fi

if ! command -v java >/dev/null 2>&1; then
  echo "FATAL: java not found. RTG Tools needs Java 8+." >&2
  exit 1
fi
echo "  java: $(java -version 2>&1 | head -1)"

# --- 1. channels -------------------------------------------------------------
# bioconda ABOVE conda-forge is the ordering bioconda documents and tests.
# strict priority stops conda mixing builds from different channels, which can
# produce a working-looking environment that crashes or computes wrongly.
echo "--- configuring conda channels ---"
"$CONDA_BASE/bin/conda" config --add channels conda-forge  >/dev/null 2>&1
"$CONDA_BASE/bin/conda" config --add channels bioconda     >/dev/null 2>&1
"$CONDA_BASE/bin/conda" config --set channel_priority strict >/dev/null 2>&1
# Anaconda's own channels carry commercial Terms of Service that would make this
# pipeline non-reproducible for anyone who has not accepted them.
for ch in defaults https://repo.anaconda.com/pkgs/main https://repo.anaconda.com/pkgs/r; do
  "$CONDA_BASE/bin/conda" config --remove channels "$ch" >/dev/null 2>&1 || true
done

# --- 2. conda environments ---------------------------------------------------
# Separate environments because GATK4, FreeBayes and hap.py have conflicting
# dependency constraints; one combined environment silently downgrades GATK.
echo "--- creating conda environments (this is the slow part) ---"
SOLVER="$CONDA_BASE/bin/mamba"
[[ -x "$SOLVER" ]] || SOLVER="$CONDA_BASE/bin/conda"

create_env() {
  local name="$1"; shift
  if [[ -d "$CONDA_BASE/envs/$name" ]]; then
    echo "  [$name] already exists — skipping"
    return 0
  fi
  echo "  [$name] creating..."
  if "$SOLVER" create -y -n "$name" "$@" > "$REPO/logs/env_${name}.log" 2>&1; then
    echo "  [$name] OK"
  else
    echo "  [$name] FAILED — see logs/env_${name}.log" >&2
    return 1
  fi
}

mkdir -p "$REPO/logs" "$REPO/tools"
create_env align   bwa bowtie2 minimap2 samtools
create_env callers gatk4 freebayes bcftools
create_env sim     art perl
create_env qc      fastqc multiqc
create_env ml      -c conda-forge python=3.11 scikit-learn pandas matplotlib seaborn \
                   jupyterlab snakemake graphviz

# --- 3. RTG Tools ------------------------------------------------------------
# A Java application distributed as a zip, not a conda-first package. The release
# has no macOS bundle, so the `nojre` build + the system Java is the route.
if [[ -x "$REPO/tools/rtg-tools-${RTG_VERSION}/rtg" ]]; then
  echo "  [rtg] already installed — skipping"
else
  echo "--- installing RTG Tools ${RTG_VERSION} ---"
  ( cd "$REPO/tools" \
    && curl -sL -o rtg.zip \
       "https://github.com/RealTimeGenomics/rtg-tools/releases/download/${RTG_VERSION}/rtg-tools-${RTG_VERSION}-nojre.zip" \
    && unzip -q -o rtg.zip && rm -f rtg.zip )
  # Off by default RTG phones home with usage stats and crash reports. Disabled
  # for reproducibility and so scoring never depends on network availability.
  printf 'RTG_TALKBACK=false\nRTG_USAGE=false\n' \
    > "$REPO/tools/rtg-tools-${RTG_VERSION}/rtg.cfg"
  "$REPO/tools/rtg-tools-${RTG_VERSION}/rtg" version >/dev/null 2>&1 \
    && echo "  [rtg] OK" || { echo "  [rtg] FAILED to run under this Java" >&2; exit 1; }
fi

# --- 4. simuG ----------------------------------------------------------------
# A single Perl script with no packaging, so it is cloned rather than installed.
if [[ -f "$REPO/tools/simuG/simuG.pl" ]]; then
  echo "  [simuG] already cloned — skipping"
else
  echo "--- cloning simuG ---"
  git clone --depth 1 https://github.com/yjx1217/simuG.git "$REPO/tools/simuG" >/dev/null 2>&1 \
    && echo "  [simuG] OK" || { echo "  [simuG] clone FAILED" >&2; exit 1; }
fi

# --- 5. references + truth sets ----------------------------------------------
if [[ -s "$REPO/data/refs/ecoli.fa.fai" && -s "$REPO/data/truth/ecoli.truth.vcf.gz" ]]; then
  echo "  [data] references and truth sets already built — skipping"
else
  echo "--- downloading and indexing references (Phase 1) ---"
  bash "$REPO/scripts/fetch_references.sh" || exit 1
  echo "--- building truth sets (Phase 2) ---"
  for g in phiX ecoli; do
    "$CONDA_BASE/envs/sim/bin/perl" "$REPO/tools/simuG/simuG.pl" \
      -refseq "$REPO/data/refs/${g}.fa" \
      -snp_count "$([[ $g == phiX ]] && echo 50 || echo 5000)" \
      -indel_count "$([[ $g == phiX ]] && echo 10 || echo 1000)" \
      -seed 20260814 -prefix "$REPO/data/truth/${g}" > "$REPO/logs/simug_${g}.log" 2>&1 || exit 1
    bash "$REPO/scripts/build_truth.sh" "$g" || exit 1
    bash "$REPO/scripts/verify_truth.sh" "$g" >/dev/null || exit 1
  done
fi

# --- 6. record versions (R9) -------------------------------------------------
echo "--- recording tool versions ---"
bash "$REPO/scripts/record_versions.sh" >/dev/null 2>&1 \
  && echo "  wrote logs/versions_latest.txt" \
  || echo "  WARNING: version recording failed" >&2

bash "$REPO/scripts/export_envs.sh" >/dev/null 2>&1 \
  && echo "  exported envs/*.yaml" \
  || echo "  WARNING: environment export failed" >&2

cat << EOF

=== Setup complete ===

Reproduce the baseline:
    snakemake --use-conda --cores 8

Verify everything:
    bash scripts/verify_all.sh
EOF
