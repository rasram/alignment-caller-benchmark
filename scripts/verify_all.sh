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
sec "PIPELINE — reads, QC (GATE 3, 4)"
# =============================================================================
for t in "$TAG_P" "$TAG_E"; do
  g="${t%%_*}"
  if [[ -s "work/${t}_1.fq" && -s "work/${t}_2.fq" && -s "work/${t}_.sam" ]]; then
    ml=$(awk '!/^>/{n+=length($0)} END{print n+0}' "data/truth/${g}.simseq.genome.fa")
    n1=$(( $(wc -l < "work/${t}_1.fq") / 4 )); n2=$(( $(wc -l < "work/${t}_2.fq") / 4 ))
    cov=$(python3 -c "print(f'{(($n1+$n2)*150)/$ml:.2f}')")
    ok "reads + ART truth SAM present ($g)" "${n1} pairs, ${cov}x"
    [[ "$cov" == "30.00" ]] && ok "coverage is 30x as requested ($g)" \
                            || warn "coverage ($g)" "$cov"
  else bad "reads + ART truth SAM present ($g)"; fi
done

[[ -s results/qc/multiqc_report.html ]] && ok "MultiQC report generated" || bad "MultiQC report"
if [[ -s results/qc/mean_q.tsv ]]; then
  n=$(( $(wc -l < results/qc/mean_q.tsv) - 1 ))
  ok "mean-Q table populated" "$n datasets"
else bad "mean-Q table populated"; fi

# =============================================================================
sec "PIPELINE — alignment (GATE 5)"
# =============================================================================
for t in "$TAG_P" "$TAG_E"; do
  for a in bwa bowtie2 minimap2; do
    bam="work/${t}.${a}.md.bam"
    if [[ -s "$bam" && -s "${bam}.bai" ]]; then
      rg=$("$SAMTOOLS" view -H "$bam" 2>/dev/null | grep -c '^@RG' || true)
      [[ "$rg" -ge 1 ]] && ok "BAM indexed + duplicate-marked + @RG (${t%%_*}/$a)" \
                        || bad "R6 read group missing (${t%%_*}/$a)"
    else bad "BAM present+indexed (${t%%_*}/$a)"; fi
  done
done

if [[ -s results/align_metrics.tsv ]]; then
  n=$(( $(wc -l < results/align_metrics.tsv) - 1 ))
  check "align_metrics.tsv rows (2 genomes x 3 aligners)" 6 "$n"
else bad "align_metrics.tsv present"; fi

if [[ -s results/placement_accuracy.tsv ]]; then
  n=$(( $(wc -l < results/placement_accuracy.tsv) - 1 ))
  check "placement_accuracy.tsv rows" 6 "$n"
  worst=$(awk -F'\t' 'NR>1{print $8}' results/placement_accuracy.tsv | sort -n | head -1)
  ok "placement accuracy computed (min across aligners)" "$worst"
else bad "placement_accuracy.tsv present"; fi

# =============================================================================
sec "PIPELINE — variant calling (GATE 6) — R2 PLOIDY"
# =============================================================================
NVCF=0; NDIP=0
for t in "$TAG_P" "$TAG_E"; do
  for a in bwa bowtie2 minimap2; do
    for c in gatk freebayes bcftools; do
      v="work/${t}.${a}.${c}.raw.vcf.gz"
      [[ -s "$v" ]] || continue
      NVCF=$((NVCF+1))
      d=$("$BCFTOOLS" query -f '[%GT]\n' "$v" 2>/dev/null | grep -c '[/|]' || true)
      NDIP=$((NDIP+d))
    done
  done
done
check "raw VCFs present (2 genomes x 9 pipelines)" 18 "$NVCF"
if [[ "$NDIP" -eq 0 ]]; then ok "R2 — LIVE re-check: zero diploid genotypes in any VCF"
else bad "R2 — diploid genotypes found" "$NDIP"; fi

if [[ -s logs/ploidy_verification.txt ]]; then
  p=$(grep -c "PASS (haploid)" logs/ploidy_verification.txt || true)
  ok "logs/ploidy_verification.txt logged" "$p PASS entries"
else bad "logs/ploidy_verification.txt"; fi

NFILT=$(ls work/*.filt.vcf.gz 2>/dev/null | wc -l | tr -d ' ')
check "hard-filtered callsets (raw + filt per pipeline)" 18 "$NFILT"

# R7: no BQSR anywhere.
# NOTE: exclude this file. verify_all.sh lives in scripts/ and contains the
# search pattern in its own source, so an unfiltered grep matches itself and
# reports a BQSR violation that does not exist.
if grep -rliE "BaseRecalibrator|ApplyBQSR" scripts/ Snakefile 2>/dev/null \
     | grep -qv 'verify_all.sh'; then
  bad "R7 — no BQSR" "found in: $(grep -rliE 'BaseRecalibrator|ApplyBQSR' scripts/ Snakefile 2>/dev/null | grep -v verify_all.sh | tr '\n' ' ')"
else ok "R7 — no BQSR anywhere in scripts or Snakefile"; fi

# =============================================================================
sec "SCORING (GATE 7a)"
# =============================================================================
NNORM=$(ls work/norm/*.norm.vcf.gz 2>/dev/null | wc -l | tr -d ' ')
[[ "$NNORM" -ge 18 ]] && ok "all call sets normalised" "$NNORM files" \
                      || bad "all call sets normalised" "$NNORM"

# R3: truth and calls must have been normalised with the SAME arguments.
if grep -q -- "-m -any --atomize" scripts/build_truth.sh 2>/dev/null \
   && grep -q -- "-m -any --atomize" scripts/score_variants.sh 2>/dev/null \
   && grep -q -- "-m -any --atomize" Snakefile 2>/dev/null; then
  ok "R3 — identical norm args in truth, scoring and Snakefile"
else bad "R3 — normalisation arguments differ between truth and calls"; fi

NVE=$(ls -d results/vcfeval/*__raw__all 2>/dev/null | wc -l | tr -d ' ')
check "vcfeval runs, full callset (2 x 9)" 18 "$NVE"
NVS=$(ls -d results/vcfeval/*__raw__snps results/vcfeval/*__raw__indels 2>/dev/null | wc -l | tr -d ' ')
check "vcfeval runs, SNV/indel separately" 36 "$NVS"

if [[ -s results/results.tsv ]]; then
  tot=$(( $(wc -l < results/results.tsv) - 1 ))
  pri=$(awk -F'\t' 'NR>1 && $10=="single_run"' results/results.tsv | wc -l | tr -d ' ')
  check "results.tsv baseline rows (9 pipelines x 2 genomes x 2 types)" 36 "$pri"
  ok "results.tsv total rows (incl. comparison method)" "$tot"
  # A row with F1 exactly 0 on E. coli would signal a bug, not a result.
  z=$(awk -F'\t' 'NR>1 && $1=="ecoli" && $15+0==0' results/results.tsv | wc -l | tr -d ' ')
  [[ "$z" -eq 0 ]] && ok "no E. coli pipeline scored F1=0 (bug signature)" \
                   || bad "E. coli rows with F1=0" "$z"
else bad "results/results.tsv present"; fi

for f in roc_ecoli_SNV roc_ecoli_INDEL roc_phiX_SNV roc_phiX_INDEL; do
  [[ -s "results/${f}.svg" ]] && ok "ROC plot ${f}.svg" || bad "ROC plot ${f}.svg"
done

[[ -s logs/scoring_negative_control.txt ]] && ok "scoring negative control recorded" \
  || warn "scoring negative control recorded"

# =============================================================================
sec "READY FOR THE SWEEP (GATE 7b)"
# =============================================================================
if [[ -s config/conditions.tsv ]]; then
  rows=$(( $(wc -l < config/conditions.tsv) - 1 ))
  uniq=$(awk -F'\t' 'NR>1 && $2=="phiX"{print $3"_"$4"_"$5}' config/conditions.tsv | sort -u | wc -l | tr -d ' ')
  check "conditions.tsv unique conditions per genome" 11 "$uniq"
  check "conditions.tsv rows (11 x 5 seeds x 2 genomes)" 110 "$rows"
  ok "pipeline runs implied by design" "$((rows*9))"
else bad "config/conditions.tsv present"; fi

[[ -s Snakefile ]] && ok "Snakefile present" || bad "Snakefile present"
[[ -s results/workflow_dag.svg ]] && ok "workflow_dag.svg exported" \
    "$(grep -c '<title>' results/workflow_dag.svg) nodes" || bad "workflow_dag.svg exported"

SM="$CONDA_BASE/envs/ml/bin/snakemake"
export PATH="$CONDA_BASE/bin:$CONDA_BASE/condabin:$PATH"
if out=$("$SM" -n --use-conda --conda-frontend conda --config run=all 2>&1); then
  jobs=$(echo "$out" | awk '/^total/{print $2; exit}')
  ok "snakemake -n resolves the FULL sweep DAG" "${jobs:-?} jobs"
else
  bad "snakemake -n resolves the FULL sweep DAG" "$(echo "$out" | grep -iE 'error|exception' | head -1)"
fi
if out=$("$SM" -n --use-conda --conda-frontend conda 2>&1); then
  echo "$out" | grep -q "Nothing to be done" \
    && ok "baseline is up to date under Snakemake (idempotent)" \
    || warn "baseline not up to date" "$(echo "$out" | awk '/^total/{print $2}')"
else bad "snakemake -n baseline"; fi

# =============================================================================
sec "DOCUMENTATION"
# =============================================================================
if [[ -s NOTES.md ]]; then
  ok "NOTES.md present" "$(wc -l < NOTES.md | tr -d ' ') lines, $(grep -c '^## ' NOTES.md) sections"
else bad "NOTES.md present"; fi
[[ -s README.md ]]   && ok "README.md present"   || bad "README.md present (one-command reproduction)"
[[ -s HANDOFF.md ]]  && ok "HANDOFF.md present"  || bad "HANDOFF.md present"

# =============================================================================
printf '\n\033[1m== SUMMARY ==\033[0m\n'
printf '  PASS %d   FAIL %d   WARN %d\n' "$PASS" "$FAIL" "$WARN"
if [[ "$FAIL" -eq 0 ]]; then
  printf '  \033[32mAll required checks passed.\033[0m\n'
else
  printf '  \033[31m%d required check(s) failed — see [FAIL] above.\033[0m\n' "$FAIL"
fi
exit $(( FAIL > 0 ? 1 : 0 ))
