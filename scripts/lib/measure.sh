#!/usr/bin/env bash
# Run ONE command and record its exact wall time, peak RSS and CPU time.
#
# Why not Snakemake's `benchmark:` alone? It SAMPLES memory with psutil on a
# timer. On macOS those samples came back empty (NA) — psutil often cannot read
# another process's memory there — and even where sampling works it misses short
# peaks. The kernel tracks the true peak itself (getrusage), and /usr/bin/time
# reports it exactly. This wrapper also times only the tool, excluding the
# conda-activation overhead that a whole-job benchmark includes (which dominates
# sub-second phiX jobs).
#
# Portable across the two `time` dialects:
#   macOS / BSD : /usr/bin/time -l -o FILE   (peak RSS in BYTES)
#   Linux / GNU : /usr/bin/time -f FMT -o FILE (peak RSS in KILOBYTES)
# Mixing those units up misreports memory by 1024× (NOTES 5.4).
#
# Usage:
#   measure.sh --out M.tsv [--stdout F] [--stderr F] -- cmd arg ...
#   measure.sh --out M.tsv -- sh -c "a | b > out"     # for pipelines
set -uo pipefail
OUT=""; SO="/dev/stdout"; SE="/dev/stderr"
while [ $# -gt 0 ]; do
  case "$1" in
    --out)    OUT="$2"; shift 2 ;;
    --stdout) SO="$2";  shift 2 ;;
    --stderr) SE="$2";  shift 2 ;;
    --)       shift; break ;;
    *) echo "measure.sh: unknown argument $1" >&2; exit 2 ;;
  esac
done
[ -n "$OUT" ] && [ $# -gt 0 ] || { echo "usage: measure.sh --out F -- cmd ..." >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"
TF="$(mktemp "${TMPDIR:-/tmp}/measure.XXXXXX")"

if [ "$(uname)" = "Darwin" ]; then
  /usr/bin/time -l -o "$TF" "$@" > "$SO" 2> "$SE"; rc=$?
  read -r secs usr sys < <(awk '/ real /{print $1, $3, $5; exit}' "$TF")
  rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576; exit}' "$TF")
else
  /usr/bin/time -o "$TF" -f "%e %U %S %M" "$@" > "$SO" 2> "$SE"; rc=$?
  read -r secs usr sys kb < <(tail -1 "$TF")
  rss=$(awk -v k="$kb" 'BEGIN{printf "%.1f", k/1024}')
fi
cpu=$(awk -v u="${usr:-0}" -v s="${sys:-0}" 'BEGIN{printf "%.2f", u+s}')
printf 's\tmax_rss_mb\tcpu_s\texit\n%s\t%s\t%s\t%s\n' "${secs:-NA}" "${rss:-NA}" "$cpu" "$rc" > "$OUT"
rm -f "$TF"
exit "$rc"
