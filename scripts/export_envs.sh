#!/usr/bin/env bash
# Export each conda environment to envs/<name>.yaml with pinned versions.
#
# Two files are written per environment, because they serve different purposes:
#
#   envs/<name>.yaml         — name + version + build string, no platform-specific
#                              URLs. This is what Snakemake's `conda:` directive uses
#                              and what a collaborator on Linux can actually install.
#
#   envs/<name>.lock.yaml    — `--explicit`-style full URL lock of this exact build
#                              set on osx-arm64. Byte-exact reproduction on an
#                              identical machine; will NOT solve on another platform.
#
# Keeping both is deliberate: the lock file proves what *these* numbers were produced
# with, the portable file lets someone else re-run the benchmark at all.
#
# Usage: bash scripts/export_envs.sh

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
mkdir -p "$REPO/envs"

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

for env in align callers sim qc ml; do
  if ! conda env list | awk '{print $1}' | grep -qx "$env"; then
    echo "SKIP $env — environment does not exist"
    continue
  fi

  # Portable: --from-history keeps only what was explicitly requested, but drops
  # versions. We want the full solved list minus the platform-locked build URLs,
  # so use the normal export and strip the `prefix:` line (it leaks a local path).
  conda env export -n "$env" --no-builds \
    | grep -v '^prefix:' > "$REPO/envs/${env}.yaml"

  # Exact: full URLs including build hashes.
  conda list -n "$env" --explicit --md5 > "$REPO/envs/${env}.lock.yaml"

  n=$(grep -c '^  - ' "$REPO/envs/${env}.yaml" 2>/dev/null || echo 0)
  echo "Exported $env  ($n packages) -> envs/${env}.yaml + envs/${env}.lock.yaml"
done
