#!/usr/bin/env bash
# Shared tool resolution for every script in this project.
#
# WHY THIS EXISTS
# The scripts were originally written to call tools by absolute path into named
# conda environments ($CONDA_BASE/envs/align/bin/bwa). That is deliberate: it
# avoids `conda run`, which adds ~1s of startup per call and raises
# BrokenPipeError when a consumer like `head` exits early (see NOTES 2).
#
# But Snakemake's per-rule `conda:` directive works by ACTIVATING an environment
# and putting its tools on PATH. A script that ignores PATH would silently keep
# using the developer's local environments, making the `conda:` declarations
# decorative — the workflow would appear portable while not being portable.
#
# So: prefer whatever is already on PATH (the Snakemake-provided environment),
# and fall back to the named conda env only when the tool is not on PATH.
# This makes the same scripts correct under `snakemake --use-conda`, under plain
# `snakemake`, and when run directly from a shell.

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"

resolve_tool() {   # resolve_tool <executable> <fallback-env-name>
  local name="$1" env="$2" p
  p="$(command -v "$name" 2>/dev/null || true)"
  if [[ -n "$p" ]]; then
    printf '%s\n' "$p"
  else
    printf '%s\n' "$CONDA_BASE/envs/${env}/bin/${name}"
  fi
}
