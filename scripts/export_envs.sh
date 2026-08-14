#!/usr/bin/env bash
# Export each conda environment to envs/<name>.yaml with pinned versions.
#
# Two files are written per environment, because they serve different purposes:
#
#   envs/<name>.yaml         — the TOP-LEVEL packages only, pinned to the exact
#                              versions in use. This is what Snakemake's `conda:`
#                              directive consumes and what a collaborator on
#                              another platform can actually install.
#
#                              It is deliberately NOT a full `conda env export`.
#                              A full export lists every transitive dependency at
#                              an exact version; re-solving that set on another
#                              machine (or another day) frequently fails, because
#                              it over-constrains packages we never asked for.
#                              Pinning what we chose and letting the solver fill
#                              in the rest is both reproducible and installable.
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

  # Portable: the packages we explicitly asked for, at their installed versions.
  # `--from-history` recovers exactly that request list (but without versions), so
  # the versions are looked up from `conda list`.
  {
    echo "name: $env"
    echo "channels:"
    echo "  - bioconda"
    echo "  - conda-forge"
    echo "dependencies:"
    conda env export -n "$env" --from-history \
      | awk '/^dependencies:/{f=1;next} /^[a-z]/{f=0} f && /^  - /{sub(/^  - /,""); sub(/=.*$/,""); print}' \
      | while read -r pkg; do
          [ -n "$pkg" ] || continue
          v=$(conda list -n "$env" --json "^${pkg}$" 2>/dev/null \
              | python3 -c "import json,sys; d=json.load(sys.stdin); print(d[0]['version'] if d else '')" 2>/dev/null)
          if [ -n "$v" ]; then echo "  - ${pkg}=${v}"; else echo "  - ${pkg}"; fi
        done
  } > "$REPO/envs/${env}.yaml"

  # Exact: full URLs including build hashes.
  conda list -n "$env" --explicit --md5 > "$REPO/envs/${env}.lock.yaml"

  n=$(grep -c '^  - ' "$REPO/envs/${env}.yaml" 2>/dev/null || echo 0)
  echo "Exported $env  ($n packages) -> envs/${env}.yaml + envs/${env}.lock.yaml"
done
