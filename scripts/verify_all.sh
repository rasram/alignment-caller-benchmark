#!/usr/bin/env bash
# END-TO-END VERIFICATION against the Definition of Done (PROJECT_BRIEF.md §6).
#
# This does NOT just check that files exist. Where it can, it re-derives the fact:
# it reads lengths out of the .fai, greps GT fields out of the VCFs, counts rows
# in results.tsv, re-runs the ploidy check, and so on. A check that only tested
# for file presence would pass on a file full of wrong numbers.
#
# Usage: bash scripts/verify_all.sh
# Exit status is non-zero if any REQUIRED check fails.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
# shellcheck disable=SC1091
source "$REPO/scripts/lib/tools.sh"
BCFTOOLS="$(resolve_tool bcftools callers)"
SAMTOOLS="$(resolve_tool samtools align)"
RTG="$REPO/tools/rtg-tools-3.13/rtg"

PASS=0; FAIL=0; WARN=0
ok()   { printf '  \033[32m[PASS]\033[0m %-58s %s\n' "$1" "${2:-}"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31m[FAIL]\033[0m %-58s %s\n' "$1" "${2:-}"; FAIL=$((FAIL+1)); }
warn() { printf '  \033[33m[WARN]\033[0m %-58s %s\n' "$1" "${2:-}"; WARN=$((WARN+1)); }
sec()  { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }

check() { # description expected actual
  if [[ "$2" == "$3" ]]; then ok "$1" "$3"; else bad "$1" "expected=$2 actual=$3"; fi
}

TAG_P=phiX_cov30_len150_err0_seed1
TAG_E=ecoli_cov30_len150_err0_seed1

# =============================================================================
sec "INFRASTRUCTURE"
# =============================================================================
[[ -d .git ]] && ok "git repository initialised" "$(git rev-list --count HEAD 2>/dev/null) commits" \
              || bad "git repository initialised"

for pat in 'work/' '\*.fq' '\*.bam' '\*.sdf/'; do
  grep -qE "^${pat}" .gitignore 2>/dev/null \
    && ok ".gitignore excludes ${pat//\\/}" \
    || bad ".gitignore excludes ${pat//\\/}"
done

# Nothing large should actually be tracked by git.
SMK_TRACKED=$(git ls-files .snakemake 2>/dev/null | wc -l | tr -d ' ')
[[ "$SMK_TRACKED" -eq 0 ]] && ok ".snakemake/ not tracked in git" \
  || bad ".snakemake/ not tracked in git" "$SMK_TRACKED files"
GITSZ=$(du -sm .git 2>/dev/null | cut -f1)
[[ "${GITSZ:-0}" -lt 200 ]] && ok "git repo size sane" "${GITSZ}MB" \
  || warn "git repo size" "${GITSZ}MB — something large may be committed"
BIGTRACKED=$(git ls-files 2>/dev/null | grep -cE '\.(fq|bam|sam|fa)$' || true)
[[ "$BIGTRACKED" -eq 0 ]] && ok "no reads/BAMs/FASTAs tracked in git" \
                          || bad "no reads/BAMs/FASTAs tracked in git" "$BIGTRACKED tracked"

for e in align callers sim qc ml; do
  if [[ -s "envs/${e}.yaml" ]]; then
    npin=$(grep -cE '^\s+- \S+=' "envs/${e}.yaml" 2>/dev/null || echo 0)
    if [[ "$npin" -gt 0 ]]; then ok "envs/${e}.yaml exported with pinned versions" "$npin pinned"
    else bad "envs/${e}.yaml exported with pinned versions" "no versions pinned"; fi
  else bad "envs/${e}.yaml exported"; fi
  [[ -s "envs/${e}.lock.yaml" ]] && ok "envs/${e}.lock.yaml exact lock present" \
                                 || warn "envs/${e}.lock.yaml exact lock present"
done

# Tools actually run, not merely present.
if v=$("$RTG" version 2>/dev/null | head -1); then ok "RTG Tools runs" "$v"; else bad "RTG Tools runs"; fi
if [[ -f tools/simuG/simuG.pl ]]; then
  ok "simuG cloned" "$(cd tools/simuG && git rev-parse --short HEAD 2>/dev/null)"
else bad "simuG cloned"; fi

VERFILE=$(ls -t logs/versions_*.txt 2>/dev/null | head -1)
if [[ -n "$VERFILE" ]]; then
  # grep -c exits 1 on zero matches, so `|| echo 0` would append a SECOND line
  # and make $nf "0\n0". Use `|| true` and let grep's own 0 stand.
  nf=$(grep -c "NOT FOUND" "$VERFILE" 2>/dev/null || true); nf=${nf:-0}
  if [[ "$nf" -eq 0 ]]; then ok "logs/versions_*.txt written, all tools found" "$(basename "$VERFILE")"
  else bad "logs/versions_*.txt has missing tools" "$nf NOT FOUND"; fi
else bad "logs/versions_*.txt written"; fi

# =============================================================================
sec "DATA — references (GATE 1)"
# =============================================================================
check "phiX length from .fai"   5386    "$(cut -f2 data/refs/phiX.fa.fai 2>/dev/null)"
check "ecoli length from .fai"  4641652 "$(cut -f2 data/refs/ecoli.fa.fai 2>/dev/null)"
check "phiX contig renamed"     phiX    "$(cut -f1 data/refs/phiX.fa.fai 2>/dev/null)"
check "ecoli contig renamed"    ecoli   "$(cut -f1 data/refs/ecoli.fa.fai 2>/dev/null)"

for g in phiX ecoli; do
  miss=""
  for f in "${g}.fa.fai" "${g}.fa.bwt" "${g}.fa.sa" "${g}.dict" "${g}.1.bt2"; do
    [[ -s "data/refs/$f" ]] || miss="$miss $f"
  done
  [[ -d "data/refs/${g}.sdf" ]] || miss="$miss ${g}.sdf"
  [[ -z "$miss" ]] && ok "$g fully indexed (fai/bwa/bowtie2/dict/sdf)" \
                   || bad "$g fully indexed" "missing:$miss"
done

# R5: contig names identical across every artefact that vcfeval touches.
for g in phiX ecoli; do
  a=$(cut -f1 "data/refs/${g}.fa.fai" 2>/dev/null)
  b=$("$BCFTOOLS" view -H "data/truth/${g}.truth.vcf.gz" 2>/dev/null | cut -f1 | sort -u | tr -d '\n')
  c=$(cut -f1 "data/truth/${g}.confident.bed" 2>/dev/null)
  d=$(grep -m1 '^@SQ' <("$RTG" sdfstats "data/refs/${g}.sdf" 2>/dev/null) >/dev/null 2>&1; \
      "$RTG" sdfstats --lengths "data/refs/${g}.sdf" 2>/dev/null | awk '/^[[:space:]]*'"$g"'/{print $1; exit}')
  if [[ "$a" == "$b" && "$a" == "$c" && "$a" == "$d" ]]; then
    ok "R5 contig names identical ($g)" "$a"
  else
    bad "R5 contig names identical ($g)" "fai=$a truth=$b bed=$c sdf=$d"
  fi
done

# =============================================================================
sec "DATA — truth sets (GATE 2)"
# =============================================================================
check "phiX truth SNVs"    50   "$("$BCFTOOLS" view -H -v snps   data/truth/phiX.truth.vcf.gz 2>/dev/null | wc -l | tr -d ' ')"
check "phiX truth indels"  10   "$("$BCFTOOLS" view -H -v indels data/truth/phiX.truth.vcf.gz 2>/dev/null | wc -l | tr -d ' ')"
check "ecoli truth SNVs"   5000 "$("$BCFTOOLS" view -H -v snps   data/truth/ecoli.truth.vcf.gz 2>/dev/null | wc -l | tr -d ' ')"
check "ecoli truth indels" 1000 "$("$BCFTOOLS" view -H -v indels data/truth/ecoli.truth.vcf.gz 2>/dev/null | wc -l | tr -d ' ')"

# Genotype column decision: must be haploid GT=1 with a sample column.
for g in phiX ecoli; do
  gts=$("$BCFTOOLS" query -f '[%GT]\n' "data/truth/${g}.truth.vcf.gz" 2>/dev/null | sort -u | tr '\n' ',')
  [[ "$gts" == "1," ]] && ok "truth genotype column is haploid GT=1 ($g)" \
                       || bad "truth genotype column ($g)" "found: $gts"
done

# Length bookkeeping: does the truth VCF describe the actual mutated genome?
for g in phiX ecoli; do
  rl=$(awk '!/^>/{n+=length($0)} END{print n+0}' "data/refs/${g}.fa" 2>/dev/null)
  ml=$(awk '!/^>/{n+=length($0)} END{print n+0}' "data/truth/${g}.simseq.genome.fa" 2>/dev/null)
  nv=$("$BCFTOOLS" view -H -v indels "data/truth/${g}.truth.vcf.gz" 2>/dev/null \
       | awk '{s+=length($5)-length($4)} END{print s+0}')
  check "truth VCF accounts for mutated-genome length ($g)" "$((ml-rl))" "$nv"
done

# Re-run a live REF spot-check rather than trusting the stored log.
for g in phiX ecoli; do
  bads=0
  while IFS=$'\t' read -r c p _ r a _; do
    [[ -n "${p:-}" ]] || continue
    end=$((p + ${#r} - 1))
    fa=$("$SAMTOOLS" faidx "data/refs/${g}.fa" "${c}:${p}-${end}" 2>/dev/null | grep -v '^>' | tr -d '\n')
    [[ "$r" == "$fa" ]] || bads=$((bads+1))
  done <<< "$("$BCFTOOLS" view -H "data/truth/${g}.truth.vcf.gz" 2>/dev/null | head -5 || true)"
  [[ "$bads" -eq 0 ]] && ok "live REF spot-check ($g, 5 variants)" \
                      || bad "live REF spot-check ($g)" "$bads mismatches"
done

for g in phiX ecoli; do
  [[ -s "data/truth/${g}.confident.bed" ]] && ok "confident BED present ($g)" \
      "$(cat "data/truth/${g}.confident.bed")" || bad "confident BED present ($g)"
done

# =============================================================================
# Which run tags are we auditing?
#   VERIFY_SCOPE=all     (default) every tag in config/conditions.tsv — the
#                        project is only complete when the full sweep has run
#   VERIFY_SCOPE=present only tags that have been executed (for partial runs)
# =============================================================================
SCOPE="${VERIFY_SCOPE:-all}"
if [[ "$SCOPE" == "present" && -s results/read_metrics.tsv ]]; then
  TAGS=($(awk -F'\t' 'NR>1{print $1}' results/read_metrics.tsv | sort -u))
else
  TAGS=($(awk -F'\t' 'NR>1{print $9}' config/conditions.tsv))
fi
NT=${#TAGS[@]}
TAGFILE="$(mktemp)"; printf '%s\n' "${TAGS[@]}" > "$TAGFILE"
printf '\n  scope: %s — auditing %d run tags (conditions × seeds × genomes)\n' "$SCOPE" "$NT"

# helper: count rows of a TSV whose column $2 value is one of our tags
count_tags() { awk -F'\t' -v col="$2" 'NR==FNR{t[$1];next} FNR>1 && ($col in t)' "$TAGFILE" "$1" | wc -l | tr -d ' '; }

# =============================================================================
sec "PIPELINE — reads, QC (GATE 3, 4)"
# =============================================================================
if [[ -s results/read_metrics.tsv ]]; then
  check "read_metrics.tsv: one row per tag" "$NT" "$(count_tags results/read_metrics.tsv 1)"
  # actual coverage must match what was requested, for EVERY tag
  offcov=$(awk -F'\t' 'NR==FNR{t[$1];next} FNR>1 && ($1 in t) {
            r=$10/$3; if (r<0.98 || r>1.02) n++ } END{print n+0}' "$TAGFILE" results/read_metrics.tsv)
  [[ "$offcov" -eq 0 ]] && ok "actual coverage within ±2% of requested, every tag" \
                        || bad "coverage off by >2%" "$offcov tags"
  nanq=$(awk -F'\t' 'FNR>1 && ($12=="" || tolower($12)=="nan")' results/read_metrics.tsv | wc -l | tr -d ' ')
  [[ "$nanq" -eq 0 ]] && ok "measured error rate (mean_p) present for every tag" \
                      || bad "mean_p missing" "$nanq tags"
else bad "results/read_metrics.tsv present"; fi
[[ -s results/qc/multiqc_report.html ]] && ok "baseline MultiQC report present" || bad "MultiQC report"

# =============================================================================
sec "PIPELINE — alignment (GATE 5)"
# =============================================================================
nbam=0; nrg=0
for t in "${TAGS[@]}"; do for a in bwa bowtie2 minimap2; do
  b="work/${t}.${a}.md.bam"
  [[ -s "$b" && -s "$b.bai" ]] && nbam=$((nbam+1))
  "$SAMTOOLS" view -H "$b" 2>/dev/null | grep -q '^@RG' && nrg=$((nrg+1))
done; done
check "duplicate-marked, indexed BAMs (3 per tag)" "$((NT*3))" "$nbam"
check "R6 — read group present in every BAM" "$((NT*3))" "$nrg"
[[ -s results/align_metrics.tsv ]] \
  && check "align_metrics.tsv rows (3 per tag)" "$((NT*3))" "$(count_tags results/align_metrics.tsv 1)" \
  || bad "align_metrics.tsv present"
if [[ -s results/placement_accuracy.tsv ]]; then
  check "placement_accuracy.tsv rows (3 per tag)" "$((NT*3))" "$(count_tags results/placement_accuracy.tsv 1)"
  worst=$(awk -F'\t' 'NR>1{print $8}' results/placement_accuracy.tsv | sort -n | head -1)
  ok "placement accuracy computed (minimum across all runs)" "$worst"
else bad "placement_accuracy.tsv present"; fi
if [[ -s results/runtime.tsv ]]; then
  check "runtime.tsv rows (3 align + 9 call per tag)" "$((NT*12))" "$(count_tags results/runtime.tsv 1)"
  nanrt=$(awk -F'\t' 'NR>1 && (tolower($5)=="nan" || tolower($6)=="nan")' results/runtime.tsv | wc -l | tr -d ' ')
  [[ "$nanrt" -eq 0 ]] && ok "every timed job has wall time AND peak RSS" \
                       || bad "timed jobs with missing measurement" "$nanrt"
else bad "results/runtime.tsv present"; fi

# =============================================================================
sec "PIPELINE — variant calling (GATE 6) — R2 PLOIDY"
# =============================================================================
NVCF=0; NDIP=0; NFILT=0
for t in "${TAGS[@]}"; do for a in bwa bowtie2 minimap2; do for c in gatk freebayes bcftools; do
  v="work/${t}.${a}.${c}.raw.vcf.gz"
  [[ -s "$v" ]] || continue
  NVCF=$((NVCF+1))
  [[ -s "work/${t}.${a}.${c}.filt.vcf.gz" ]] && NFILT=$((NFILT+1))
  d=$("$BCFTOOLS" query -f '[%GT]\n' "$v" 2>/dev/null | grep -c '[/|]' || true)
  NDIP=$((NDIP+d))
done; done; done
check "raw VCFs present (9 per tag)" "$((NT*9))" "$NVCF"
if [[ "$NDIP" -eq 0 ]]; then ok "R2 — LIVE re-check: zero diploid genotypes in $NVCF VCFs"
else bad "R2 — diploid genotypes found" "$NDIP"; fi
if [[ -s logs/ploidy_verification.txt ]]; then
  p=$(awk 'NR==FNR{t[$1];next} ($1 in t) && /PASS \(haploid\)/' "$TAGFILE" logs/ploidy_verification.txt | wc -l | tr -d ' ')
  check "logs/ploidy_verification.txt PASS entries (9 per tag)" "$((NT*9))" "$p"
else bad "logs/ploidy_verification.txt"; fi
check "hard-filtered callsets (9 per tag)" "$((NT*9))" "$NFILT"
# NOTE: exclude this file. verify_all.sh lives in scripts/ and contains the
# search pattern in its own source, so an unfiltered grep matches itself.
if grep -rliE "BaseRecalibrator|ApplyBQSR" scripts/ Snakefile 2>/dev/null | grep -qv 'verify_all.sh'; then
  bad "R7 — no BQSR" "found in: $(grep -rliE 'BaseRecalibrator|ApplyBQSR' scripts/ Snakefile | grep -v verify_all.sh | tr '\n' ' ')"
else ok "R7 — no BQSR anywhere in scripts or Snakefile"; fi

# =============================================================================
sec "SCORING (GATE 7a)"
# =============================================================================
NNORM=0
for t in "${TAGS[@]}"; do NNORM=$((NNORM + $(ls work/norm/${t}.*.raw.norm.vcf.gz 2>/dev/null | wc -l))); done
check "call sets normalised (9 per tag)" "$((NT*9))" "$NNORM"
if grep -q -- "-m -any --atomize" scripts/build_truth.sh && grep -q -- "-m -any --atomize" scripts/score_variants.sh \
   && grep -q -- "-m -any --atomize" Snakefile; then
  ok "R3 — identical norm args in truth, scoring and Snakefile"
else bad "R3 — normalisation arguments differ between truth and calls"; fi
NVE=0; NVS=0
for t in "${TAGS[@]}"; do
  NVE=$((NVE + $(ls -d results/vcfeval/${t}__*__raw__all 2>/dev/null | wc -l)))
  NVS=$((NVS + $(ls -d results/vcfeval/${t}__*__raw__snps results/vcfeval/${t}__*__raw__indels 2>/dev/null | wc -l)))
done
check "vcfeval runs, full call set (9 per tag)" "$((NT*9))" "$NVE"
check "vcfeval runs, SNV/indel pre-split (18 per tag)" "$((NT*18))" "$NVS"
if [[ -s results/results.tsv ]]; then
  read -r PRI MISS ZERO < <("$CONDA_BASE/envs/ml/bin/python" - "$TAGFILE" <<'PY'
import csv, sys
want = set(open(sys.argv[1]).read().split())
rows = [r for r in csv.DictReader(open("results/results.tsv"), delimiter="\t")
        if r["scoring_method"] == "single_run" and r["callset"] == "raw"
        and f'{r["genome"]}_cov{r["coverage"]}_len{r["read_length"]}_err{r["qs_shift"]}_seed{r["seed"]}' in want]
feat = ["mapping_rate", "mean_mapq", "mean_depth", "placement_accuracy", "align_seconds",
        "call_seconds", "peak_rss_mb", "actual_coverage", "mean_q", "mean_p"]
miss = sum(1 for r in rows for c in feat if str(r[c]).strip().lower() in ("", "nan", "na"))
zero = sum(1 for r in rows if r["genome"] == "ecoli" and float(r["f1"]) == 0)
print(len(rows), miss, zero)
PY
)
  check "results.tsv raw single_run rows (18 per tag)" "$((NT*18))" "$PRI"
  NFROWS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
           $h["callset"]=="filt" && $h["scoring_method"]=="single_run"' results/results.tsv | wc -l | tr -d ' ')
  check "results.tsv hard-filtered rows scored (18 per tag)" "$((NT*18))" "$NFROWS"
  [[ "$MISS" -eq 0 ]] && ok "no missing feature cells in results.tsv" \
                      || bad "missing feature cells in results.tsv" "$MISS"
  [[ "$ZERO" -eq 0 ]] && ok "no E. coli pipeline scored F1 = 0 (bug signature)" \
                      || bad "E. coli rows with F1 = 0" "$ZERO"
else bad "results/results.tsv present"; fi
[[ -s logs/scoring_negative_control.txt ]] && ok "scoring negative control recorded" \
  || warn "scoring negative control recorded"

# =============================================================================
sec "WORKFLOW (GATE 7b)"
# =============================================================================
if [[ -s config/conditions.tsv ]]; then
  uniq=$(awk -F'\t' 'NR>1 && $2=="phiX"{print $3"_"$4"_"$5}' config/conditions.tsv | sort -u | wc -l | tr -d ' ')
  check "conditions.tsv unique conditions per genome" 11 "$uniq"
  check "conditions.tsv rows (11 × 5 seeds × 2 genomes)" 110 "$(( $(wc -l < config/conditions.tsv) - 1 ))"
else bad "config/conditions.tsv present"; fi
[[ -s results/workflow_dag.svg ]] && ok "workflow DAG exported" || bad "workflow DAG exported"
SM="$CONDA_BASE/envs/ml/bin/snakemake"
export PATH="$CONDA_BASE/bin:$CONDA_BASE/condabin:$PATH"
if out=$("$SM" -n --config run=all 2>&1); then
  if echo "$out" | grep -q "Nothing to be done"; then
    ok "full 990-run sweep is complete and up to date under Snakemake"
  else
    j=$(echo "$out" | awk '/^total/{print $2; exit}')
    [[ "$SCOPE" == "all" ]] && bad "full sweep not complete" "$j jobs outstanding" \
                            || ok "full sweep DAG resolves" "$j jobs outstanding"
  fi
else bad "snakemake -n --config run=all" "$(echo "$out" | grep -iE 'error|exception' | head -1)"; fi
rm -f "$TAGFILE"

# =============================================================================
sec "ANALYSIS, MODEL, FIGURES (STEPS 7–9)"
# =============================================================================
for f in summary_by_condition effects titv_experiment scoring_method_effect \
         runtime_by_condition hard_filter_effect; do
  [[ -s "results/analysis/${f}.tsv" ]] && ok "results/analysis/${f}.tsv" \
                                       || bad "results/analysis/${f}.tsv present"
done
if [[ -s results/analysis/effects.tsv ]]; then
  # every E. coli condition must have been TESTED (>= 2 seeds), not merely summarised
  read -r NT_EC NTESTED < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
       $h["genome"]=="ecoli"{n++; if($h["aligner_q"]!="" && $h["aligner_q"]!="nan") t++}
       END{print n+0, t+0}' results/analysis/effects.tsv)
  check "E. coli condition × type cells with a significance test" "22" "$NTESTED"
fi
for f in validation importances recommendations loco_by_condition; do
  [[ -s "results/model/${f}.tsv" ]] && ok "results/model/${f}.tsv" \
                                    || bad "results/model/${f}.tsv present"
done
[[ -s results/model/tree_rules.txt ]] && ok "decision tree rules exported" || bad "tree_rules.txt"
# a regression tree ranks whole leaves equal; the selection metric must treat a tie
# as a tie (expected regret), never break it by row order — guard the fix
if [[ -s results/model/recommendations.tsv ]] && head -1 results/model/recommendations.tsv | grep -q model_top_n; then
  ok "model recommendations are tie-aware" "top-ranked SETS, not idxmax picks"
else bad "model recommendations are tie-aware (model_top_n column)"; fi
NFIG=$(ls results/figures/F*.png 2>/dev/null | wc -l | tr -d ' ')
check "report figures F1–F7" 7 "$NFIG"
STALE=$(find results/figures -name 'F*.png' ! -newer results/results.tsv | wc -l | tr -d ' ')
check "figures newer than results.tsv (not stale)" 0 "$STALE"

# error-mechanism tests (Step 7b): one row per pipeline x depth / pipeline x {30x,100x}
for f in fp_near_indel:54 fn_repeats:18; do
  n=$(( $(wc -l < "results/analysis/${f%%:*}.tsv" 2>/dev/null || echo 1) - 1 ))
  check "results/analysis/${f%%:*}.tsv rows" "${f##*:}" "$n"
done
[[ -s results/analysis/phix_errors.tsv ]] && ok "results/analysis/phix_errors.tsv" \
  "$(( $(wc -l < results/analysis/phix_errors.tsv) - 1 )) distinct phiX errors outside 5x" \
  || bad "results/analysis/phix_errors.tsv present"

# =============================================================================
sec "DOCUMENTATION"
# =============================================================================
if [[ -s NOTES.md ]]; then
  ok "NOTES.md present" "$(wc -l < NOTES.md | tr -d ' ') lines, $(grep -c '^## ' NOTES.md) sections"
else bad "NOTES.md present"; fi
[[ -s README.md ]]   && ok "README.md present"   || bad "README.md present (one-command reproduction)"
[[ -s HANDOFF.md ]]  && ok "HANDOFF.md present"  || bad "HANDOFF.md present"
[[ -s docs/FINAL_REPORT.md ]] && ok "final report (markdown)" || bad "docs/FINAL_REPORT.md present"
[[ -s docs/FINAL_REPORT.pdf ]] && ok "final report (PDF)" || bad "docs/FINAL_REPORT.pdf present"
[[ -s docs/FINAL_REPORT.docx ]] && ok "final report (DOCX)" || bad "docs/FINAL_REPORT.docx present"
# the report's tables and quoted numbers must be generated, never left as placeholders
if [[ -s docs/FINAL_REPORT.md ]] && grep -q '{{' docs/FINAL_REPORT.md; then
  bad "report contains unfilled placeholders"
else ok "report has no unfilled placeholders"; fi
# ...and must have been generated from the CURRENT results and template
if [[ -s docs/FINAL_REPORT.md && docs/FINAL_REPORT.md -nt results/results.tsv \
      && docs/FINAL_REPORT.md -nt docs/report/FINAL_REPORT.template.md \
      && docs/FINAL_REPORT.pdf -nt docs/FINAL_REPORT.md ]]; then
  ok "report is up to date" "newer than results.tsv and its template; PDF newer than markdown"
else bad "report is up to date (rebuild: snakemake --config run=all)"; fi

# =============================================================================
printf '\n\033[1m== SUMMARY ==\033[0m\n'
printf '  PASS %d   FAIL %d   WARN %d\n' "$PASS" "$FAIL" "$WARN"
if [[ "$FAIL" -eq 0 ]]; then
  printf '  \033[32mAll required checks passed.\033[0m\n'
else
  printf '  \033[31m%d required check(s) failed — see [FAIL] above.\033[0m\n' "$FAIL"
fi
exit $(( FAIL > 0 ? 1 : 0 ))
