# =============================================================================
# Read alignment & variant-calling benchmark — Snakemake workflow
#
# One rule per tool per input, driven by config/conditions.tsv. Fine-grained
# rules are what let Snakemake parallelise the sweep and resume it after a
# failure; wrapping the phase scripts would give a 5-node DAG with neither.
#
# CORRECTNESS RULES CARRIED OVER (see NOTES.md):
#   R1  reads simulated FROM the mutated genome, aligned TO the original reference
#   R2  ploidy 1 in every caller, verified from the GT field (not from the score)
#   R3  identical normalisation on truth and every call set, incl. --atomize
#   R5  contig names inherited from the reference everywhere
#   R6  read groups on every BAM (bowtie2 needs different syntax)
#   R7  no BQSR
#   R8  identical threads across aligners, identical filter across callers, and
#       timed jobs run with the machine to themselves (see "Fair timing" below)
#   R9  seed recorded in every filename
#
# FAIR TIMING (R8)
#   Runtime is a reported metric. In a parallel sweep, a tool timed while seven
#   other jobs compete for the CPU is measured under arbitrary contention, so the
#   comparison becomes noise. Every timed rule therefore requests the whole
#   `machine` resource (8 units); every other rule gets 1 by default
#   (profiles/default/config.yaml). A timed job can only start when nothing else
#   is running, and nothing else can start while it runs — reproducing the
#   conditions the Phase 5/6 baseline was timed under.
#
# METRICS (scatter → gather)
#   Each metric is computed per unit into its own file, then aggregated. The
#   original scripts appended to one shared TSV, which corrupts under parallel
#   execution because concurrent appends interleave.
#
# Usage (the profile in profiles/default/ is picked up automatically):
#   snakemake -n                       dry run, baseline only
#   snakemake                          execute the baseline
#   snakemake -n --config run=all      dry run of the full 990-run sweep
#   snakemake --config run=all         execute the sweep
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
    # Baseline = the starred condition at seed 1 — what Phase 7a ran.
    SELECTED = [c for c in CONDITIONS
                if c["is_baseline"] == "yes" and c["seed"] == str(config["baseline_seed"])]
elif RUN_MODE == "tags":
    # An explicit subset, e.g. the pilot:  --config run=tags tags=TAG1,TAG2
    _want = [x for x in str(config.get("tags", "")).split(",") if x]
    _known = {c["tag"] for c in CONDITIONS}
    _bad = [x for x in _want if x not in _known]
    if _bad:
        raise ValueError(f"tags not in {config['conditions']}: {_bad}")
    SELECTED = [c for c in CONDITIONS if c["tag"] in set(_want)]
elif RUN_MODE == "all":
    SELECTED = CONDITIONS
else:
    raise ValueError(f"run must be baseline | tags | all, not {RUN_MODE!r}")

TAGS = [c["tag"] for c in SELECTED]
COND_BY_TAG = {c["tag"]: c for c in CONDITIONS}     # all, so any tag resolves

ALIGNERS = config["aligners"]
CALLERS = config["callers"]
VTYPES = ["snps", "indels"]

THREADS = config["threads"]
RG = config["read_group"]
MACHINE = config["machine_units"]      # a timed job claims all of these
TIMING_SEEDS = {str(s) for s in config.get("timing_seeds", [1])}


def timed_machine(wildcards):
    """Whole machine for timing seeds (clean measurement); one unit otherwise.
    See config.yaml `timing_seeds` for why only some seeds are timed cleanly."""
    return MACHINE if COND_BY_TAG[wildcards.tag]["seed"] in TIMING_SEEDS else 1


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
def _final_targets():
    t = ["results/results.tsv", "logs/ploidy_verification.txt"]
    if RUN_MODE == "all":
        # Statistics need all seeds, so the analysis chain and the report are
        # built only for the full sweep — never for the baseline or a subset.
        t += ["docs/FINAL_REPORT.md", "docs/FINAL_REPORT.pdf", "docs/FINAL_REPORT.docx"]
    return t


rule all:
    input: _final_targets()


# =============================================================================
# PHASE 3 — read simulation
# =============================================================================
rule simulate_reads:
    """R1: reads come FROM the mutated genome. R9: seed is in the filename.

    All three outputs are temp(): Snakemake deletes them once their last consumer
    has run (the three aligners + read_metrics for the FASTQs; the three
    placement jobs for the truth SAM). Measured, this cuts the sweep's disk
    footprint from ~50 GB to ~15 GB. ART is deterministic for a given seed, so
    anything deleted is exactly regenerable."""
    input:
        mutated = lambda w: f"data/truth/{genome_of(w.tag)}.simseq.genome.fa",
    output:
        r1  = temp("work/{tag}_1.fq"),
        r2  = temp("work/{tag}_2.fq"),
        sam = temp("work/{tag}_.sam"),      # ART truth SAM -> placement accuracy
    params:
        prefix = "work/{tag}_",
        cov    = lambda w: COND_BY_TAG[w.tag]["coverage"],
        length = lambda w: COND_BY_TAG[w.tag]["read_length"],
        qs     = lambda w: COND_BY_TAG[w.tag]["qs_shift"],
        seed   = lambda w: COND_BY_TAG[w.tag]["seed"],
    log: "logs/run/art_{tag}.log"
    conda: "envs/sim.yaml"
    shell:
        r"""
        art_illumina -ss HS25 -sam -na \
          -i {input.mutated} -p -l {params.length} -f {params.cov} \
          -m 350 -s 50 -qs {params.qs} -qs2 {params.qs} \
          -rs {params.seed} -o {params.prefix} > {log} 2>&1
        """


rule read_metrics:
    """Actual coverage and measured error rate (mean_p) per run. mean_p, not
    ART's -qs knob, is the physically meaningful error feature (NOTES 4.4)."""
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        mutated = lambda w: f"data/truth/{genome_of(w.tag)}.simseq.genome.fa",
    output: "work/metrics/{tag}.reads.tsv"
    conda: "envs/metrics.yaml"
    shell:
        "python3 scripts/read_metrics.py --tag {wildcards.tag} --r1 {input.r1} "
        "--r2 {input.r2} --genome-fa {input.mutated} --out {output}"


# =============================================================================
# PHASE 5 — alignment (R1, R6, R8). Timed: benchmark + whole-machine resource.
# =============================================================================
rule align_bwa:
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        ref = lambda w: ref(w.tag),
        idx = lambda w: ref(w.tag) + ".bwt",
    output:
        sam  = temp("work/{tag}.bwa.sam"),
        meas = "benchmarks/measure/align/{tag}.bwa.tsv",
    threads: THREADS
    resources: machine = timed_machine
    benchmark: "benchmarks/align/{tag}.bwa.tsv"
    log: "logs/run/align_{tag}_bwa.log"
    conda: "envs/align.yaml"
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} --stdout {output.sam} --stderr {log} -- \
              bwa mem -t {threads} -R '{RG}' {input.ref} {input.r1} {input.r2}"""


rule align_bowtie2:
    """Bowtie2 will not take BWA's tab-delimited @RG string; it needs --rg-id
    plus one --rg per field, or the header is malformed and GATK rejects it."""
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        idx = lambda w: f"data/refs/{genome_of(w.tag)}.1.bt2",
    output:
        sam  = temp("work/{tag}.bowtie2.sam"),
        meas = "benchmarks/measure/align/{tag}.bowtie2.tsv",
    params: prefix = lambda w: f"data/refs/{genome_of(w.tag)}"
    threads: THREADS
    resources: machine = timed_machine
    benchmark: "benchmarks/align/{tag}.bowtie2.tsv"
    log: "logs/run/align_{tag}_bowtie2.log"
    conda: "envs/align.yaml"
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} --stdout {output.sam} --stderr {log} -- \
              bowtie2 -p {threads} \
              --rg-id s1 --rg SM:sim --rg PL:ILLUMINA --rg LB:lib1 \
              -x {params.prefix} -1 {input.r1} -2 {input.r2}"""


rule align_minimap2:
    """-ax sr is the short-read preset; without it minimap2 uses long-read
    defaults and places 150 bp reads badly."""
    input:
        r1 = "work/{tag}_1.fq", r2 = "work/{tag}_2.fq",
        ref = lambda w: ref(w.tag),
    output:
        sam  = temp("work/{tag}.minimap2.sam"),
        meas = "benchmarks/measure/align/{tag}.minimap2.tsv",
    threads: THREADS
    resources: machine = timed_machine
    benchmark: "benchmarks/align/{tag}.minimap2.tsv"
    log: "logs/run/align_{tag}_minimap2.log"
    conda: "envs/align.yaml"
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} --stdout {output.sam} --stderr {log} -- \
              minimap2 -ax sr -t {threads} -R '{RG}' {input.ref} {input.r1} {input.r2}"""


rule sort_bam:
    input:  "work/{tag}.{aligner}.sam"
    output: temp("work/{tag}.{aligner}.sorted.bam")
    threads: THREADS
    log: "logs/run/sort_{tag}_{aligner}.log"
    conda: "envs/align.yaml"
    shell: "samtools sort -@ {threads} -o {output} {input} 2> {log}"


rule mark_duplicates:
    """Simulated reads contain no PCR duplicates, so this marks ~0%. Kept for
    pipeline realism; its inertness is reported, not hidden (NOTES 5.8)."""
    input:  "work/{tag}.{aligner}.sorted.bam"
    output:
        bam = "work/{tag}.{aligner}.md.bam",
        met = "work/metrics/{tag}.{aligner}.markdup.txt",
    log: "logs/run/markdup_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""gatk MarkDuplicates -I {input} -O {output.bam} -M {output.met} \
              --VALIDATION_STRINGENCY LENIENT > {log} 2>&1"""


rule index_bam:
    input:  "work/{tag}.{aligner}.md.bam"
    output: "work/{tag}.{aligner}.md.bam.bai"
    conda:  "envs/align.yaml"
    shell:  "samtools index {input}"


rule align_metrics:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        met = "work/metrics/{tag}.{aligner}.markdup.txt",
    output: "work/metrics/{tag}.{aligner}.align.tsv"
    conda: "envs/metrics.yaml"
    shell:
        "bash scripts/align_metrics_one.sh {wildcards.tag} {wildcards.aligner} "
        "{input.bam} {input.met} > {output}"


rule placement_accuracy:
    """Pure aligner metric. Converts the ART truth SAM from MUTATED to REFERENCE
    coordinates before comparing (NOTES 3.4) — without that, 99% reads as 9%."""
    input:
        truth  = "work/{tag}_.sam",
        bam    = "work/{tag}.{aligner}.md.bam",
        bai    = "work/{tag}.{aligner}.md.bam.bai",
        indels = lambda w: f"data/truth/{genome_of(w.tag)}.refseq2simseq.INDEL.vcf",
    output: "work/metrics/{tag}.{aligner}.placement.tsv"
    log: "logs/run/placement_{tag}_{aligner}.log"
    conda: "envs/metrics.yaml"
    shell:
        "python3 scripts/placement_accuracy.py --truth-sam {input.truth} "
        "--indel-vcf {input.indels} --bam {input.bam} --aligner {wildcards.aligner} "
        "--tag {wildcards.tag} --out {output} --single > {log} 2>&1"


# =============================================================================
# PHASE 6 — variant calling (R2 ploidy 1, R7 no BQSR). Timed.
# =============================================================================
rule call_gatk:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
        dic = lambda w: f"data/refs/{genome_of(w.tag)}.dict",
    output:
        vcf  = "work/{tag}.{aligner}.gatk.raw.vcf.gz",
        meas = "benchmarks/measure/call/{tag}.{aligner}.gatk.tsv",
    threads: THREADS
    resources: machine = timed_machine
    benchmark: "benchmarks/call/{tag}.{aligner}.gatk.tsv"
    log: "logs/run/gatk_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} --stdout {log} --stderr {log}.err -- \
              gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} \
              --sample-ploidy 1 --native-pair-hmm-threads {threads}"""


rule call_freebayes:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
    output:
        vcf  = "work/{tag}.{aligner}.freebayes.raw.vcf.gz",
        meas = "benchmarks/measure/call/{tag}.{aligner}.freebayes.tsv",
    resources: machine = timed_machine
    benchmark: "benchmarks/call/{tag}.{aligner}.freebayes.tsv"
    log: "logs/run/freebayes_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    # `bash -o pipefail`: under plain sh a pipeline's status is the LAST
    # command's, so a FreeBayes crash would let bgzip "succeed" and leave an
    # empty, valid-looking VCF behind.
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} -- bash -o pipefail -c \
              "freebayes -f {input.ref} -p 1 {input.bam} 2> {log} | bgzip > {output.vcf}" """


rule call_bcftools:
    input:
        bam = "work/{tag}.{aligner}.md.bam",
        bai = "work/{tag}.{aligner}.md.bam.bai",
        ref = lambda w: ref(w.tag),
    output:
        vcf  = "work/{tag}.{aligner}.bcftools.raw.vcf.gz",
        meas = "benchmarks/measure/call/{tag}.{aligner}.bcftools.tsv",
    resources: machine = timed_machine
    benchmark: "benchmarks/call/{tag}.{aligner}.bcftools.tsv"
    log: "logs/run/bcftools_{tag}_{aligner}.log"
    conda: "envs/callers.yaml"
    shell:
        r"""bash scripts/lib/measure.sh --out {output.meas} -- bash -o pipefail -c \
              "bcftools mpileup -f {input.ref} -a AD,DP -Ou {input.bam} 2> {log} | bcftools call -mv --ploidy 1 -Oz -o {output.vcf} 2>> {log}" """


rule index_vcf:
    """One indexing rule for every VCF under work/, including work/norm/.
    Snakemake wildcards match '/' by default, so a separate rule for norm/ would
    be AMBIGUOUS with this one rather than more specific."""
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
        total=$(bcftools view -H {input.vcf} | wc -l | tr -d ' ')
        dip=$(bcftools query -f '[%GT]\n' {input.vcf} | grep -c '[/|]' || true)
        if [ "$dip" -ne 0 ]; then
          echo "FATAL: {wildcards.caller} emitted $dip diploid genotypes (R2)" >&2
          exit 1
        fi
        echo "{wildcards.tag} {wildcards.aligner} {wildcards.caller} records=$total diploid=0 PASS (haploid)" > {output}
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
        ploidy = "logs/ploidy/{tag}.{aligner}.{caller}.ploidy.txt",   # gate on R2
    output: "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz"
    log: "logs/run/norm_{tag}_{aligner}_{caller}_{callset}.log"
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


# RTG defaults to every core and 90% of RAM per JVM. Under a parallel sweep that
# oversubscribes the machine badly, so threads and heap are set explicitly.
RTG_ENV = "RTG_MEM=3g"

rule vcfeval_all:
    """PRIMARY scoring: ONE vcfeval run on the full call set. RTG emits its own
    snp_roc/non_snp_roc breakdown, which is what collect_results.py reads.
    Pre-splitting before scoring degrades vcfeval's haplotype-aware comparison
    into a context-free one (NOTES 7.6)."""
    input:
        calls = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz",
        tbi   = "work/norm/{tag}.{aligner}.{caller}.{callset}.norm.vcf.gz.tbi",
        truth = lambda w: f"data/truth/{genome_of(w.tag)}.truth.vcf.gz",
        sdf   = lambda w: f"data/refs/{genome_of(w.tag)}.sdf",
        bed   = lambda w: f"data/truth/{genome_of(w.tag)}.confident.bed",
    output:
        summary = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__all/summary.txt",
    params:
        outdir = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__all",
        rtg    = config["rtg"],
    threads: 2
    log: "logs/run/vcfeval_{tag}_{aligner}_{caller}_{callset}_all.log"
    shell:
        r"""
        rm -rf {params.outdir}
        {RTG_ENV} {params.rtg} vcfeval -b {input.truth} -c {input.calls} \
          -t {input.sdf} -e {input.bed} --vcf-score-field=QUAL \
          --threads {threads} -o {params.outdir} > {log} 2>&1
        """


rule vcfeval_split:
    """The brief's per-type method, retained for comparison only (NOTES 7.6)."""
    input:
        calls = "work/norm/{tag}.{aligner}.{caller}.{callset}.{vtype}.vcf.gz",
        tbi   = "work/norm/{tag}.{aligner}.{caller}.{callset}.{vtype}.vcf.gz.tbi",
        truth = lambda w: f"work/norm/{genome_of(w.tag)}.truth.{w.vtype}.vcf.gz",
        ttbi  = lambda w: f"work/norm/{genome_of(w.tag)}.truth.{w.vtype}.vcf.gz.tbi",
        sdf   = lambda w: f"data/refs/{genome_of(w.tag)}.sdf",
        bed   = lambda w: f"data/truth/{genome_of(w.tag)}.confident.bed",
    output:
        summary = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__{vtype}/summary.txt",
    params:
        outdir = "results/vcfeval/{tag}__{aligner}__{caller}__{callset}__{vtype}",
        rtg    = config["rtg"],
    threads: 2
    log: "logs/run/vcfeval_{tag}_{aligner}_{caller}_{callset}_{vtype}.log"
    shell:
        r"""
        rm -rf {params.outdir}
        {RTG_ENV} {params.rtg} vcfeval -b {input.truth} -c {input.calls} \
          -t {input.sdf} -e {input.bed} --vcf-score-field=QUAL \
          --threads {threads} -o {params.outdir} > {log} 2>&1
        """


# =============================================================================
# GATHER — one table per metric, then the master results table
# =============================================================================
def _cat_tsv(inputs, out):
    """Concatenate single-row TSVs that share a header."""
    with open(out, "w") as w:
        for i, path in enumerate(sorted(inputs)):
            with open(path) as r:
                lines = r.read().splitlines()
            w.write("\n".join(lines if i == 0 else lines[1:]) + "\n")


rule aggregate_read_metrics:
    input:  expand("work/metrics/{tag}.reads.tsv", tag=TAGS)
    output: "results/read_metrics.tsv"
    run:    _cat_tsv(input, output[0])


rule aggregate_align_metrics:
    input:  expand("work/metrics/{tag}.{aligner}.align.tsv", tag=TAGS, aligner=ALIGNERS)
    output: "results/align_metrics.tsv"
    run:    _cat_tsv(input, output[0])


rule aggregate_placement:
    input:  expand("work/metrics/{tag}.{aligner}.placement.tsv", tag=TAGS, aligner=ALIGNERS)
    output: "results/placement_accuracy.tsv"
    run:    _cat_tsv(input, output[0])


rule aggregate_runtime:
    """Primary timing = scripts/lib/measure.sh (tool only, exact peak RSS).
    Snakemake's own benchmark wall time is carried alongside as job_seconds, so
    the conda-activation overhead it includes is visible rather than hidden."""
    input:
        expand("benchmarks/measure/align/{tag}.{aligner}.tsv", tag=TAGS, aligner=ALIGNERS),
        expand("benchmarks/measure/call/{tag}.{aligner}.{caller}.tsv",
               tag=TAGS, aligner=ALIGNERS, caller=CALLERS),
    output: "results/runtime.tsv"
    conda: "envs/python.yaml"
    params: seeds = ",".join(sorted(TIMING_SEEDS))
    shell: "python3 scripts/aggregate_runtime.py --timing-seeds {params.seeds} --out {output} {input}"


rule aggregate_ploidy:
    input:
        expand("logs/ploidy/{tag}.{aligner}.{caller}.ploidy.txt",
               tag=TAGS, aligner=ALIGNERS, caller=CALLERS)
    output: "logs/ploidy_verification.txt"
    run:
        with open(output[0], "w") as w:
            w.write(f"# R2 ploidy verification — {len(input)} call sets, all must be haploid\n")
            for p in sorted(input):
                w.write(open(p).read())


rule collect_results:
    """--strict: fail if any metric column is empty, so a missing join can never
    again produce silently blank feature columns."""
    input:
        summaries = [f"results/vcfeval/{tag}__{a}__{c}__raw__{v}/summary.txt"
                     for tag in TAGS for a in ALIGNERS for c in CALLERS
                     for v in ["all"] + VTYPES],
        # The hard-filtered call sets (QUAL>=20 && DP>=5, identical for every
        # caller — R8). Scored with the primary single-run method only. Until
        # this was added, the hard_filter rule existed but no target requested
        # its output, so the workflow never built a filtered call set at all.
        filtered = [f"results/vcfeval/{tag}__{a}__{c}__filt__all/summary.txt"
                    for tag in TAGS for a in ALIGNERS for c in CALLERS],
        reads   = "results/read_metrics.tsv",
        align   = "results/align_metrics.tsv",
        place   = "results/placement_accuracy.tsv",
        runtime = "results/runtime.tsv",
        ploidy  = "logs/ploidy_verification.txt",
    output: "results/results.tsv"
    params: tags = " ".join(TAGS)
    conda: "envs/python.yaml"
    shell:
        "python3 scripts/collect_results.py --set both --strict "
        "--tags {params.tags} --out {output}"


# =============================================================================
# STEPS 7-10 — analysis, model, figures, report (full sweep only)
# =============================================================================
rule analyse_sweep:
    """Blocked two-way ANOVA + Friedman per condition, BH-corrected (Step 7)."""
    input:
        results = "results/results.tsv",
        runtime = "results/runtime.tsv",
        conds   = config["conditions"],
    output:
        "results/analysis/summary_by_condition.tsv",
        "results/analysis/effects.tsv",
        "results/analysis/titv_experiment.tsv",
        "results/analysis/scoring_method_effect.tsv",
        "results/analysis/runtime_by_condition.tsv",
        "results/analysis/hard_filter_effect.tsv",
        "results/analysis/headline.txt",
    conda: "envs/analysis.yaml"
    shell: "python3 scripts/analyse_sweep.py > /dev/null"


rule diagnose_errors:
    """Mechanism tests for the error patterns the sweep revealed (Step 7b):
    Bowtie2 false SNVs vs true indels, the missed-variant floor vs mappability,
    and every phiX error listed by position."""
    input:
        results = "results/results.tsv",
        bams    = expand("work/ecoli_cov{cov}_len150_err0_seed1.{aligner}.md.bam",
                         cov=[30, 100], aligner=ALIGNERS),
    output:
        "results/analysis/fp_near_indel.tsv",
        "results/analysis/fn_repeats.tsv",
        "results/analysis/phix_errors.tsv",
    conda: "envs/metrics.yaml"
    shell: "python3 scripts/diagnose_errors.py > /dev/null"


rule fit_model:
    """Decision tree + forest importances; held-out-seed and LOCO validation (Step 8)."""
    input: "results/results.tsv"
    output:
        "results/model/validation.tsv", "results/model/importances.tsv",
        "results/model/recommendations.tsv", "results/model/loco_by_condition.tsv",
        "results/model/tree_rules.txt", "results/model/tree.pkl",
    conda: "envs/analysis.yaml"
    shell: "python3 scripts/fit_model.py > /dev/null"


FIGURES = ["F1_f1_vs_coverage", "F2_f1_vs_read_length", "F3_f1_vs_error_rate",
           "F4_heatmaps_by_coverage", "F5_decision_tree", "F6_runtime_vs_accuracy",
           "F7_variance_shares"]

rule make_figures:
    input:
        "results/analysis/summary_by_condition.tsv", "results/analysis/effects.tsv",
        "results/analysis/runtime_by_condition.tsv", "results/read_metrics.tsv",
        "results/model/tree.pkl",
    output: expand("results/figures/{f}.png", f=FIGURES)
    conda: "envs/analysis.yaml"
    shell: "python3 scripts/make_figures.py"


rule build_report:
    """Every table in the report is GENERATED here from results/ (Step 10)."""
    input:
        template = "docs/report/FINAL_REPORT.template.md",
        analysis = rules.analyse_sweep.output,
        diagnose = rules.diagnose_errors.output,
        model    = rules.fit_model.output,
        figures  = rules.make_figures.output,
        metrics  = ["results/read_metrics.tsv", "results/align_metrics.tsv",
                    "results/placement_accuracy.tsv"],
    output: "docs/FINAL_REPORT.md"
    conda: "envs/analysis.yaml"
    shell: "python3 scripts/build_report.py"


rule render_report:
    """PDF via the system XeLaTeX (TeX is not a conda package worth pinning here);
    DOCX via pandoc alone. Figures are referenced relative to docs/."""
    input:  "docs/FINAL_REPORT.md"
    output:
        pdf  = "docs/FINAL_REPORT.pdf",
        docx = "docs/FINAL_REPORT.docx",
    conda: "envs/report.yaml"
    shell:
        r"""
        cd docs
        pandoc FINAL_REPORT.md -o FINAL_REPORT.pdf --pdf-engine=xelatex \
            -V geometry:margin=2.2cm -V mainfont="Helvetica" -V monofont="Menlo" \
            -V fontsize=10pt -V colorlinks=true --resource-path=.:.. --toc
        pandoc FINAL_REPORT.md -o FINAL_REPORT.docx --resource-path=.:.. --toc
        """
