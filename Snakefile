# =============================================================================
# Read alignment & variant-calling benchmark — Snakemake workflow (Phase 7b)
#
# Converts the Phase 3-7a shell pipeline into a dependency graph driven by
# config/conditions.tsv.
#
# WHY FINE-GRAINED RULES RATHER THAN WRAPPING THE SCRIPTS
# The scripts in scripts/ each do a whole phase (all 3 aligners, or all 9
# caller runs) in one invocation. Wrapping them would give a 5-node DAG with no
# real parallelism and no per-pipeline resumability: one failed caller would
# force re-running all nine. So each rule here is one tool on one input, which
# is what lets Snakemake parallelise the sweep and resume it after a failure.
#
# The cost of that choice is duplicated logic — the command lines here must stay
# in step with the scripts. GATE 7b guards exactly this: the baseline is re-run
# through Snakemake and the numbers must match Phase 7a exactly.
#
# CORRECTNESS RULES CARRIED OVER (see NOTES.md):
#   R1  reads simulated FROM the mutated genome, aligned TO the original reference
#   R2  ploidy 1 in every caller, verified from the GT field
#   R3  identical normalisation on truth and every call set, incl. --atomize
#   R5  contig names inherited from the reference everywhere
#   R6  read groups on every BAM (bowtie2 needs different syntax)
#   R7  no BQSR
#   R8  identical thread count across aligners; identical filter across callers
#   R9  seed recorded in every filename
#
# Usage:
#   snakemake -n                      dry run, baseline only (default)
#   snakemake --cores 8               execute the baseline
#   snakemake -n --config run=all     dry run of the FULL 990-run sweep
#   snakemake --dag | dot -Tsvg > results/workflow_dag.svg
# =============================================================================

import csv
import os

configfile: "config/config.yaml"


# --- load the sweep design table --------------------------------------------
CONDITIONS = []
with open(config["conditions"]) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        CONDITIONS.append(row)

RUN_MODE = config.get("run", "baseline")
if RUN_MODE == "baseline":
    # Baseline = the * condition, and only seed 1. This is what Phase 7a ran.
    SELECTED = [c for c in CONDITIONS
                if c["is_baseline"] == "yes" and c["seed"] == str(config["baseline_seed"])]
else:
    SELECTED = CONDITIONS

TAGS = [c["tag"] for c in SELECTED]
COND_BY_TAG = {c["tag"]: c for c in SELECTED}

ALIGNERS = config["aligners"]
CALLERS = config["callers"]
VTYPES = ["snps", "indels"]

THREADS = config["threads"]
RG = config["read_group"]


wildcard_constraints:
    genome   = r"phiX|ecoli",
    aligner  = r"bwa|bowtie2|minimap2",
    caller   = r"gatk|freebayes|bcftools",
    callset  = r"raw|filt",
    vtype    = r"snps|indels",
    tag      = r"[A-Za-z0-9]+_cov\d+_len\d+_err-?\d+_seed\d+",


def genome_of(tag):
    return COND_BY_TAG[tag]["genome"]


def ref(tag):
    return f"data/refs/{genome_of(tag)}.fa"


# --- final targets ----------------------------------------------------------
def all_targets():
    t = ["results/results.tsv"]
    for tag in TAGS:
        for a in ALIGNERS:
            for c in CALLERS:
                t.append(f"results/vcfeval/{tag}__{a}__{c}__raw__all/summary.txt")
                for v in VTYPES:
                    t.append(f"results/vcfeval/{tag}__{a}__{c}__raw__{v}/summary.txt")
    return t


rule all:
    input:
        all_targets()


# =============================================================================
# PHASE 3 — read simulation
# =============================================================================
rule simulate_reads:
    """R1: reads come FROM the mutated genome. R9: seed is in the filename."""
    input:
        mutated = lambda w: f"data/truth/{genome_of(w.tag)}.simseq.genome.fa",
    output:
        r1  = "work/{tag}_1.fq",
        r2  = "work/{tag}_2.fq",
        sam = "work/{tag}_.sam",      # ART truth SAM -> placement accuracy
    params:
        prefix = "work/{tag}_",
        cov    = lambda w: COND_BY_TAG[w.tag]["coverage"],
        length = lambda w: COND_BY_TAG[w.tag]["read_length"],
        qs     = lambda w: COND_BY_TAG[w.tag]["qs_shift"],
        seed   = lambda w: COND_BY_TAG[w.tag]["seed"],
    log:
        "logs/art_{tag}.log",
    conda:
        "envs/sim.yaml"
    shell:
        r"""
        art_illumina -ss HS25 -sam -na \
          -i {input.mutated} -p -l {params.length} -f {params.cov} \
          -m 350 -s 50 -qs {params.qs} -qs2 {params.qs} \
          -rs {params.seed} -o {params.prefix} > {log} 2>&1
        """


# =============================================================================
# PHASE 5 — alignment (R1, R6, R8)
# =============================================================================
rule align_bwa:
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        ref = lambda w: ref(w.tag),
        idx = lambda w: ref(w.tag) + ".bwt",
    output:
        temp("work/{tag}.bwa.sam")
    threads: THREADS
    log: "logs/align_{tag}_bwa.log"
    conda: "envs/align.yaml"
    shell:
        r"""bwa mem -t {threads} -R '{RG}' {input.ref} {input.r1} {input.r2} \
              > {output} 2> {log}"""


rule align_bowtie2:
    """Bowtie2 will not take BWA's tab-delimited @RG string; it needs --rg-id
    plus one --rg per field, or the header is malformed and GATK rejects it."""
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        idx = lambda w: f"data/refs/{genome_of(w.tag)}.1.bt2",
    output:
        temp("work/{tag}.bowtie2.sam")
    params:
        prefix = lambda w: f"data/refs/{genome_of(w.tag)}",
    threads: THREADS
    log: "logs/align_{tag}_bowtie2.log"
    conda: "envs/align.yaml"
    shell:
        r"""bowtie2 -p {threads} \
              --rg-id s1 --rg SM:sim --rg PL:ILLUMINA --rg LB:lib1 \
              -x {params.prefix} -1 {input.r1} -2 {input.r2} \
              > {output} 2> {log}"""


rule align_minimap2:
    """-ax sr is the short-read preset; without it minimap2 uses long-read
    defaults and places 150 bp reads badly."""
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        ref = lambda w: ref(w.tag),
    output:
        temp("work/{tag}.minimap2.sam")
    threads: THREADS
    log: "logs/align_{tag}_minimap2.log"
    conda: "envs/align.yaml"
    shell:
        r"""minimap2 -ax sr -t {threads} -R '{RG}' {input.ref} \
              {input.r1} {input.r2} > {output} 2> {log}"""


rule sort_bam:
    input:  "work/{tag}.{aligner}.sam"
    output: temp("work/{tag}.{aligner}.sorted.bam")
    threads: THREADS
    log: "logs/sort_{tag}_{aligner}.log"
    conda: "envs/align.yaml"
    shell: "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule mark_duplicates:
    """Simulated reads contain no PCR duplicates, so this marks ~0%. Kept for
    pipeline realism; its inertness is reported, not hidden (NOTES 5.8)."""
    input:  "work/{tag}.{aligner}.sorted.bam"
    output:
        bam = "work/{tag}.{aligner}.md.bam",
        met = "logs/{tag}.{aligner}.md.metrics",
    log: "logs/markdup_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""gatk MarkDuplicates -I {input} -O {output.bam} -M {output.met} \
              --VALIDATION_STRINGENCY LENIENT > {log} 2>&1"""


rule index_bam:
    input:  "work/{tag}.{aligner}.md.bam"
    output: "work/{tag}.{aligner}.md.bam.bai"
    conda:  "envs/align.yaml"
    shell:  "samtools index {input}"


# =============================================================================
# PHASE 6 — variant calling (R2 ploidy 1, R7 no BQSR)
# =============================================================================
rule call_gatk:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
        dic = lambda w: f"data/refs/{genome_of(w.tag)}.dict",
    output: "work/{tag}.{aligner}.gatk.raw.vcf.gz"
    threads: THREADS
    log: "logs/gatk_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} \
              --sample-ploidy 1 --native-pair-hmm-threads {threads} > {log} 2>&1"""


rule call_freebayes:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
    output: "work/{tag}.{aligner}.freebayes.raw.vcf.gz"
    log: "logs/freebayes_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""freebayes -f {input.ref} -p 1 {input.bam} 2> {log} \
              | bgzip > {output}"""


rule call_bcftools:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
    output: "work/{tag}.{aligner}.bcftools.raw.vcf.gz"
    log: "logs/bcftools_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""bcftools mpileup -f {input.ref} -a AD,DP -Ou {input.bam} 2> {log} \
              | bcftools call -mv --ploidy 1 -Oz -o {output} 2>> {log}"""


rule index_vcf:
    """One indexing rule for every VCF under work/, including work/norm/.
    Snakemake wildcards match '/' by default, so a separate rule for the norm/
    subdirectory would be AMBIGUOUS with this one rather than more specific."""
    input:  "work/{prefix}.vcf.gz"
    output: "work/{prefix}.vcf.gz.tbi"
    conda:  "envs/callers.yaml"
    shell:  "bcftools index -t -f {input}"


rule verify_ploidy:
    """R2: read the GT field directly. vcfeval will NOT catch a caller emitting
    1/1 (it scores a perfect F1 against haploid truth), so the score cannot be
    used as a ploidy check. Fails the workflow on any diploid genotype."""
    input:
        vcf = "work/{tag}.{aligner}.{caller}.raw.vcf.gz",
        tbi = "work/{tag}.{aligner}.{caller}.raw.vcf.gz.tbi",
    output: "logs/ploidy/{tag}.{aligner}.{caller}.ploidy.txt"
    conda: "envs/callers.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        total=$(bcftools view -H {input.vcf} | wc -l | tr -d ' ')
        dip=$(bcftools query -f '[%GT]\n' {input.vcf} | grep -c '[/|]' || true)
        echo "{wildcards.tag} {wildcards.aligner} {wildcards.caller} records=$total diploid=$dip" > {output}
        if [ "$dip" -ne 0 ]; then
          echo "FATAL: {wildcards.caller} emitted $dip diploid genotypes (R2)" >&2
          exit 1
        fi
        echo "PASS (haploid)" >> {output}
        """


rule hard_filter:
    """R8: identical filter expression for every caller. Only QUAL and DP are
    used because they are the only fields all three emit comparably."""
    input:  "work/{tag}.{aligner}.{caller}.raw.vcf.gz"
    output: "work/{tag}.{aligner}.{caller}.filt.vcf.gz"
    params: expr = config["filter_expr"]
    conda:  "envs/callers.yaml"
    shell:  """bcftools view -i '{params.expr}' -Oz -o {output} {input}"""


# =============================================================================
# PHASE 7a — normalisation (R3) and GA4GH scoring (R4)
# =============================================================================
rule normalise_calls:
    """R3: byte-identical to the normalisation applied to the truth set.
    --atomize decomposes FreeBayes' MNV/complex records (NOTES 6.6)."""
    input:
        vcf = "work/{tag}.{aligner}.{caller}.{callset}.vcf.gz",
        ref = lambda w: ref(w.tag),
        # gate scoring on the ploidy check having passed
        ploidy = "logs/ploidy/{tag}.{aligner}.{caller}.ploidy.txt",
    output: "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz"
    log: "logs/norm_{tag}_{aligner}_{caller}_{callset}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""bcftools norm -f {input.ref} -m -any --atomize \
              -Oz -o {output} {input.vcf} 2> {log}"""


rule split_calls_by_type:
    input:
        vcf = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz",
        tbi = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz.tbi",
    output: "work/norm/{tag}.{aligner}.{caller}.{callset}.{vtype}.vcf.gz"
    conda:  "envs/callers.yaml"
    shell:  "bcftools view -v {wildcards.vtype} -Oz -o {output} {input.vcf}"


rule split_truth_by_type:
    input:  "data/truth/{genome}.truth.vcf.gz"
    output: "work/norm/{genome}.truth.{vtype}.vcf.gz"
    conda:  "envs/callers.yaml"
    shell:  "bcftools view -v {wildcards.vtype} -Oz -o {output} {input}"


rule vcfeval_all:
    """PRIMARY scoring: ONE vcfeval run on the full call set. RTG emits its own
    snp_roc/non_snp_roc breakdown, which is what collect_results.py reads.
    Pre-splitting before scoring degrades vcfeval's haplotype-aware comparison
    into a context-free one (NOTES 7.6)."""
    input:
        calls = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz",
        tbi   = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz.tbi",
        truth = lambda w: f"data/truth/{genome_of(w.tag)}.truth.vcf.gz",
        sdf   = lambda w: directory(f"data/refs/{genome_of(w.tag)}.sdf"),
        bed   = lambda w: f"data/truth/{genome_of(w.tag)}.confident.bed",
    output:
        summary = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__all/summary.txt",
    params:
        outdir = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__all",
        rtg    = config["rtg"],
    log: "logs/vcfeval_{tag}_{aligner}_{caller}_{callset}_all.log"
    shell:
        r"""
        rm -rf {params.outdir}
        {params.rtg} vcfeval -b {input.truth} -c {input.calls} \
          -t {input.sdf} -e {input.bed} --vcf-score-field=QUAL \
          -o {params.outdir} > {log} 2>&1
        """


rule vcfeval_split:
    """The brief's per-type method, retained for comparison only. See NOTES 7.6
    for why its numbers differ from vcfeval_all and are not used as primary."""
    input:
        calls = "work/norm/{tag}.{aligner}.{caller}.{callset}.{vtype}.vcf.gz",
        tbi   = "work/norm/{tag}.{aligner}.{caller}.{callset}.{vtype}.vcf.gz.tbi",
        truth = lambda w: f"work/norm/{genome_of(w.tag)}.truth.{w.vtype}.vcf.gz",
        ttbi  = lambda w: f"work/norm/{genome_of(w.tag)}.truth.{w.vtype}.vcf.gz.tbi",
        sdf   = lambda w: directory(f"data/refs/{genome_of(w.tag)}.sdf"),
        bed   = lambda w: f"data/truth/{genome_of(w.tag)}.confident.bed",
    output:
        summary = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__{vtype}/summary.txt",
    params:
        outdir = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__{vtype}",
        rtg    = config["rtg"],
    log: "logs/vcfeval_{tag}_{aligner}_{caller}_{callset}_{vtype}.log"
    shell:
        r"""
        rm -rf {params.outdir}
        {params.rtg} vcfeval -b {input.truth} -c {input.calls} \
          -t {input.sdf} -e {input.bed} --vcf-score-field=QUAL \
          -o {params.outdir} > {log} 2>&1
        """


rule collect_results:
    input:
        [f"results/vcfeval/{tag}__{a}__{c}__raw__{v}/summary.txt"
         for tag in TAGS for a in ALIGNERS for c in CALLERS
         for v in ["all"] + VTYPES]
    output:
        "results/results.tsv"
    conda: "envs/python.yaml"
    shell:
        "python3 scripts/collect_results.py --set raw --out {output}"
