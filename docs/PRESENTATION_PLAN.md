# Mid-semester presentation — build specification

A slide-by-slide blueprint for the mid-sem talk. Every number here is taken from
`results/results.tsv`, `results/align_metrics.tsv` and `logs/` — nothing is invented. Where a
result is provisional it is marked, and the wording to use is given.

**Target: 12 minutes speaking + 3 minutes questions, 14 slides.** A cut-list for an 8-minute
slot is at the end.

---

## 0. Before you build anything: the one decision that shapes the talk

You have two possible stories, and they lead to different presentations.

| Story A — "I ran a benchmark" | Story B — "I built a benchmark that cannot lie to me" |
|---|---|
| Leads with the F1 tables | Leads with the design, then the results as evidence it works |
| Assessor's question: *is 0.9953 vs 0.9906 meaningful?* | Assessor's question: *how do you know your numbers are right?* — which you can answer |
| Weak, because the baseline differences are tiny and single-seed | Strong, because the correctness work is the real intellectual content |

**Build Story B.** At mid-sem you have a working, verified pipeline and only *one* baseline
condition. The differences between pipelines are real but small and not yet
statistically supported. If you front-load the tables, you invite exactly the question you
cannot yet answer. If you front-load the method, the tables become supporting evidence and the
missing statistics become "the next phase", which is where they belong.

The spine of the talk is therefore: **the benchmark is only worth as much as its truth set and
its scoring — here is how I made both trustworthy, and here is what they show so far.**

---

## 1. Slide-by-slide

Timings are cumulative targets. Speaker notes are what to *say*, not what to put on the slide.

---

### Slide 1 — Title (0:00–0:20)

**On the slide**
- Benchmarking Read Aligners and Variant Callers on Known-Truth Genomes
- Your name, course, date
- One-line subtitle: *3 aligners × 3 callers = 9 pipelines, scored against variants we injected ourselves*

**Visual:** none, or a faint background of the workflow rulegraph.

**Say:** Nothing clever. Name the problem and move on within 20 seconds.

---

### Slide 2 — The question (0:20–1:20)

**On the slide**
- To find genetic differences you need two tools: an **aligner** (where did this read come from?) and a **caller** (is this position genuinely different?)
- There are many of each, all widely used, all claiming to work
- **Which combination should you actually use — and does it depend on your data?**

**Visual:** a simple 3×3 grid, aligners down the side, callers across the top, all nine cells
empty. You will fill this grid in later — reusing the same graphic makes the talk feel
designed rather than assembled.

**Say:** Frame it as a practical decision a lab actually faces. Avoid tool trivia here. The
point is that the choice is usually made by habit or by what a tutorial used, not by evidence.

---

### Slide 3 — Why you cannot benchmark on real data (1:20–2:20)

**On the slide**
- To grade a pipeline you must already know the right answer
- For a real biological sample, **nobody does**
- So: manufacture a sample where the answer is known **by construction**

**Visual:** THE key diagram of the talk. Four boxes, left to right, with arrows:

```
  Reference genome ──inject known SNPs+indels──▶ Mutated genome  +  TRUTH VCF
                                                        │
                                                        │ simulate reads FROM here
                                                        ▼
                                                     Reads
                                                        │
                                                        │ align TO the original reference
                                                        ▼
                                                  Called variants ──compare──▶ TRUTH VCF
```

Make the two directional arrows a different colour and label them **"FROM mutated"** and
**"TO original"**.

**Say:** This is the conceptual core. Spend a full minute. Explain that because *we* wrote the
mutation list, every call is gradeable as correct, spurious, or missed — and that this is the
only way to get a denominator.

**Then land the hook for slide 11:** "That FROM/TO distinction is the single easiest thing in
this project to get backwards — and if you do, every pipeline reports zero variants and they
all look identical. I'll come back to that."

---

### Slide 4 — Experimental design (2:20–3:20)

**On the slide**

| | |
|---|---|
| Aligners | BWA-MEM, Bowtie2, minimap2 |
| Callers | GATK4 HaplotypeCaller, FreeBayes, BCFtools |
| Genomes | phiX174 (5,386 bp) · *E. coli* K-12 (4,641,652 bp) |
| Truth | 50 SNV + 10 indel · 5,000 SNV + 1,000 indel |
| Baseline | 30× coverage, 150 bp paired-end, HiSeq 2500 error profile |
| Scoring | GA4GH standard (`rtg vcfeval`) |

- **Held constant so the comparison is fair:** identical reads to all nine pipelines, identical thread count, identical filter logic, no read trimming, no BQSR

**Visual:** the table above. Keep it dense but readable — this is a reference slide.

**Say:** Explain the two genomes have different jobs. phiX is the **debugging** organism —
runs end to end in seconds, so every stage is proven there before touching *E. coli*. *E. coli*
is ~860× larger and has real repetitive sequence, so it can actually distinguish aligners.

Mention **no BQSR** and say why in one sentence: it needs a database of known variants that
doesn't exist for these organisms, and bootstrapping one would give GATK a preprocessing step
its competitors don't get. **This is a fairness decision, not an omission** — say that
explicitly, because an examiner who knows GATK Best Practices will notice its absence.

---

### Slide 5 — The pipeline (3:20–4:00)

**On the slide**
- simulate → align (×3) → mark duplicates → call (×3) → normalise → score
- Implemented twice: as shell scripts, then as a Snakemake workflow
- 990-run parameter sweep already wired up; **only the baseline executed so far**

**Visual:** `results/workflow_rulegraph.svg` — the *rulegraph*, not the full DAG. The rulegraph
has ~15 nodes and reads cleanly on a projector. The full DAG has 822 and will be an unreadable
hairball. Do not use it.

**Say:** Keep this short — it's orientation, not content. One sentence worth making: the
workflow is reproducible from a single command, and rebuilding the baseline from scratch
through Snakemake reproduced the shell-script results **byte for byte**. Two independent
implementations agreeing is the strongest evidence you have that neither has a silly bug.

---

### Slide 6 — Where the project stands (4:00–4:40)

**On the slide** — a progress table, honest about what is not done:

| Phase | Status |
|---|---|
| Environment, tools, references | ✅ |
| Truth sets generated + verified | ✅ |
| Read simulation, QC | ✅ |
| Alignment + alignment metrics | ✅ |
| Variant calling, all 9 pipelines | ✅ |
| GA4GH scoring, ROC curves | ✅ |
| Reproducible workflow + full sweep design | ✅ |
| **Parameter sweep executed** | ⏳ next |
| **Predictive model** | ⏳ later |

**Say:** Be direct that the sweep has not run. Then give the number that makes it credible:
*measured* cost is ~2.25 hours for all 990 runs; the constraint is disk (~80 GB), not time.
Knowing your own compute budget precisely reads as competence.

---

### Slide 7 — Results 1: the alignment layer (4:40–5:40)

**On the slide** — *E. coli*, baseline:

| Aligner | Mapping | Placement acc. (±10 bp) | Mean MAPQ | Runtime | Peak RAM |
|---|---|---|---|---|---|
| BWA-MEM | 100% | 99.03% | 59.1 | 5.2 s | 382 MB |
| Bowtie2 | 99.98% | 99.01% | 41.2 | 17.6 s | 67 MB |
| minimap2 | 100% | 99.02% | 59.1 | 1.9 s | 491 MB |

**Say — three points, in this order:**

1. **Placement accuracy is a pure aligner metric.** Because the simulator records where every
   read really came from, you can measure the fraction the aligner put back within ±10 bp
   *without involving a caller at all*. That cleanly separates "the aligner placed reads well"
   from "the caller reasoned well".
2. **All three are effectively tied at ~99.0%.** The residual 1% is reads from repetitive
   sequence where the true location is not recoverable from 150 bp — all three hit the same
   information limit. **Say plainly: at this baseline, alignment is an easy problem.**
3. **Runtime and memory are inverted** — minimap2 is ~9× faster than Bowtie2 but uses ~7× the
   memory. A genuine engineering trade-off, not a defect.

**Do not** compare mean MAPQ across tools as if it were a quality score — BWA caps at 60,
Bowtie2 at 42. Flag that on the slide as a footnote; it sets up slide 9.

---

### Slide 8 — Results 2: the 3×3 matrices (5:40–6:40)

**On the slide** — the same 3×3 grid from slide 2, now filled, twice:

***E. coli* — SNV F1**

| | GATK | FreeBayes | BCFtools |
|---|---|---|---|
| BWA-MEM | **0.9953** | **0.9953** | 0.9948 |
| Bowtie2 | 0.9909 | 0.9906 | 0.9914 |
| minimap2 | 0.9946 | **0.9953** | 0.9944 |

***E. coli* — indel F1**

| | GATK | FreeBayes | BCFtools |
|---|---|---|---|
| BWA-MEM | **0.9975** | 0.9970 | 0.9970 |
| Bowtie2 | 0.9935 | 0.9890 | 0.9815 |
| minimap2 | 0.9970 | **0.9975** | 0.9970 |

Colour-scale the cells (best green → worst red). The Bowtie2 row will visibly be the weak one,
which does the argument for you.

**Say:** Define F1 in one line — harmonic mean of precision and recall, dominated by whichever
is worse. Give the one-sentence reason it's the headline: *a caller that reports only its single
most confident variant gets precision 1.0 and recall 0.0002 — F1 correctly calls that useless.*

Then note **phiX scored 1.0000 for all nine pipelines** and explain it in one breath: 5,386 bp
with no repeats at 30× is unambiguously solvable, so phiX is a smoke test, not a comparison.
Pre-empt the "is that a bug?" question — see slide 12.

---

### Slide 9 — The headline finding (6:40–7:40)

**On the slide, large:**

> ### At this baseline, the **aligner** matters more than the **caller**

| Effect (marginal mean F1 spread) | SNV | Indel |
|---|---|---|
| Across aligners | 0.0041 | 0.0092 |
| Across callers | 0.0002 | 0.0042 |
| **Aligner effect ÷ caller effect** | **20.5×** | **2.2×** |

- Reframed as errors: worst pipeline makes **2.0× more SNV errors** and **7.4× more indel errors** than the best
- Every Bowtie2 row is the weakest, in both variant types

**Say:** This is your one memorable claim — spend time here. F1 differences in the third
decimal place sound trivial; **converting to error counts is what makes them land.** "0.9975
versus 0.9815" is forgettable; "7.4× more indel errors" is not.

Then give the mechanism as a *hypothesis*, and label it as one:

> Bowtie2's MAPQ scale caps at 42 where the others cap at 60. Callers filter on MAPQ using
> tool-agnostic thresholds, so the same cut-off is a stricter filter on Bowtie2's output.
> Consistent with the data: Bowtie2+GATK has the worst SNV false-negative count (90 vs 47) with
> **zero** false positives — the signature of a caller discarding evidence, not of an aligner
> misplacing reads.

**Critically — say the limitation out loud before anyone asks:** this is **one seed, one
condition**. The SNV spread is 0.0047 and nothing yet establishes that it exceeds seed-to-seed
noise. The sweep's 5 seeds provide the first variance estimate. Volunteering this converts your
weakest point into evidence of judgement.

---

### Slide 10 — ROC curves (7:40–8:20)

**On the slide:** `results/roc_ecoli_INDEL.svg` (indel separates the pipelines more clearly
than SNV — use the one that shows the effect).

**Say:** Explain why the ROC matters more than the table: a single F1 is one operating point at
one quality threshold, and **QUAL is not calibrated the same way across the three callers** —
so comparing single F1 values partly compares the callers' quality scales rather than their
ability to find variants. The ROC sweeps the threshold and removes that arbitrariness.

Point at one concrete case: **Bowtie2+FreeBayes has the highest SNV recall of any pipeline
(0.9912) but also by far the most false positives (50, where six pipelines have zero).** It's
the one aggressive pipeline in the set — invisible in the F1 table, obvious on the curve.

---

### Slide 11 — **The money slide: three silent failures** (8:20–10:00)

This is the intellectual core. Budget the most time here. Title it something like
*"Three ways this benchmark could have produced confident, wrong numbers"*.

**On the slide** — three rows, each: what could go wrong → what it would have looked like → how
it was caught.

| # | The trap | What it would have shown | Caught by |
|---|---|---|---|
| 1 | **Ploidy** — these organisms are haploid; all three callers default to diploid | Nothing. No error message. And scoring can't detect it — a caller emitting `1/1` scores a **perfect F1** against haploid truth | Grepping the genotype field directly on all 18 call sets |
| 2 | **Coordinate systems** — the simulator's truth file uses *mutated-genome* coordinates; aligners report *reference* coordinates. Same contig name, no error | Placement accuracy **8.9%** instead of **99.0%** — all three aligners looking catastrophically broken | Converting coordinates explicitly, then self-validating the conversion against the simulator's own records |
| 3 | **Variant representation** — FreeBayes merges nearby variants into one record (`215 AAA>TAT` = two SNPs) | Splitting SNVs from indels before scoring silently discards those records *and the variants inside them* — FreeBayes' recall understated as a pure artefact | Decomposing complex records before scoring, applied identically to truth and calls |

**Visual:** for trap 2, a small before/after bar — 8.9% vs 99.0% — is extremely effective.
For trap 3, show the concrete example: truth has `215 A>T` and `217 A>T`; FreeBayes writes
`215 AAA>TAT`; both describe identical sequence.

**Say:** The unifying lesson, and the sentence to end the slide on:

> **Every one of these fails silently.** None produces an error message. Each produces
> plausible-looking numbers that are wrong. That's why the project spends more effort on
> verification than on running the tools.

This slide is what separates "I ran some bioinformatics software" from "I understand what
makes a benchmark valid". If you have to cut elsewhere to protect this slide, do.

---

### Slide 12 — How the results are verified (10:00–10:50)

**On the slide**
- **Negative control:** deliberately corrupt the call set and confirm the scoring punishes it

| Call set | F1 |
|---|---|
| Unmodified | 1.0000 |
| All positions shifted +5 bp | 0.0000 |
| Every SNV allele changed | 0.1667 |
| Half the calls removed | 0.6667 (recall exactly 0.500) |

- **Truth-set integrity:** the net length change of all injected indels matches the mutated genome's actual size difference exactly (+14 bp phiX, +187 bp *E. coli*)
- **Automated audit:** 83 checks, re-derived from the data rather than checking files exist
- **Two independent implementations** (shell + Snakemake) produce identical results

**Say:** Lead with the negative control, because it answers the phiX = 1.0 question directly:
*a perfect score is only believable if the same machinery can be shown to fail when it should.*
Shifting positions by 5 bp collapses F1 to zero; removing half the calls gives recall of exactly
0.500. The scoring discriminates, so 1.0 on phiX is a real result about an easy genome.

The length-bookkeeping check is worth 15 seconds: it ties the truth file to the actual FASTA
the reads came from. Counts and spot-checks both pass even if those two disagree — and if they
disagreed, every number in the project would be wrong with no other symptom.

---

### Slide 13 — What's next (10:50–11:30)

**On the slide**
- **Run the sweep:** 11 conditions × 5 seeds × 2 genomes × 9 pipelines = **990 runs**
  - Coverage 5 → 100× · read length 75 → 150 bp · error rate 0 → −10 quality shift
  - Measured cost: ~2.25 h compute; the binding constraint is ~80 GB disk
- **Then model:** which pipeline wins under which conditions?
- **Expected payoff:** all three aligners are tied at 30× — differences should emerge at low coverage and short reads, which is exactly where the sweep looks

**Visual:** a small schematic of the one-factor-at-a-time design — baseline in the centre, three
axes radiating out.

**Say:** Be explicit that the design is one-factor-at-a-time and that its known limitation is
**it cannot detect interactions** — if a pipeline only degrades when coverage is low *and* reads
are short, this design misses it. A full grid would be 72 conditions per genome instead of 11.
Naming your own design's weakness before the examiner does is worth a lot.

---

### Slide 14 — Summary (11:30–12:00)

**On the slide — three lines only:**
1. Built a benchmark where the correct answer is known by construction, and **verified that it can detect its own failure**
2. At 30× coverage, the **aligner choice matters more than the caller** — up to 7.4× fewer indel errors — though this is one seed and needs the sweep to confirm
3. Three classes of silent failure identified and controlled; workflow reproduces byte-for-byte

**Say:** Don't read the slide. Say the one thing you want remembered: *the hard part of a
benchmark isn't running the tools, it's making sure the numbers mean what you think they mean.*

---

## 2. Anticipated questions, with answers

Prepare these. The first four are near-certain.

**Q: The differences are in the third decimal place. Are they real?**
> Honest answer: not yet established. One seed per condition, so there's no variance estimate.
> The spread across aligners is 20× the spread across callers for SNVs, which is suggestive,
> and it's consistent across both variant types and both directions — but the sweep's 5 seeds
> are what will tell us whether it exceeds noise. That's the first thing the next phase produces.

**Q: Why is phiX 1.0 everywhere? Isn't that a bug?**
> That's exactly the right question, and I tested it rather than assuming. Negative control:
> shifting every call by 5 bp drops F1 to 0.0000, removing half the calls gives recall of
> exactly 0.500. The scoring discriminates. phiX is 5,386 bp with no repetitive sequence — at
> 30× every variant is unambiguously recoverable. It's a smoke test, not a comparison.

**Q: Why no BQSR? That's part of GATK Best Practices.**
> It needs a database of known variants to separate real variation from systematic error, and
> none exists for phiX or *E. coli*. I could bootstrap one, but that gives GATK a data-driven
> preprocessing step FreeBayes and BCFtools don't get — any GATK advantage afterwards couldn't
> be attributed to the caller. It's a deliberate fairness trade-off, and it does mean the GATK
> arm here isn't full Best Practices.

**Q: Simulated reads are easier than real data. Does this transfer?**
> Partly. Real data has PCR duplicates, adapter contamination, coverage bias and mapping
> artefacts that ART doesn't simulate — so absolute numbers here are optimistic. What does
> transfer is the *ranking* under controlled conditions, because all nine pipelines see
> byte-identical input. One specific caveat: my mutation spectrum uses the simulator's default
> transition/transversion ratio of ~0.5, where real bacteria run 1–2, so absolute recall isn't
> directly comparable to a real resequencing project.

**Q: Why these three aligners / callers?**
> They're the most widely used in each category and they represent genuinely different
> algorithms — BWA and Bowtie2 are both FM-index based but tuned differently, minimap2 uses
> minimizer seeding; GATK does local reassembly, FreeBayes is Bayesian haplotype-based,
> BCFtools is a position-wise pileup model. So the comparison spans real algorithmic diversity,
> not three variations of one idea.

**Q: How do you know your truth set is correct?**
> Three independent checks. Counts match what was requested. Reference bases in the truth file
> match the actual genome at those coordinates. And the strongest one: the net length change of
> all injected indels equals the mutated genome's actual size difference exactly — +14 bp for
> phiX, +187 for *E. coli*. That last one ties the truth file to the exact sequence the reads
> were generated from.

**Q: What's the hardest thing you hit?**
> The coordinate-system mismatch. The simulator's truth file and the aligner output use the same
> contig name and differ in length by only the net indel balance, so nothing errors — placement
> accuracy just silently reads 8.9% instead of 99.0%, and worse toward one end of the genome.
> It's the kind of bug you only find if you're suspicious of a number that looks plausible.

---

## 3. Design and production notes

**Slide count and density.** 14 slides in 12 minutes is ~50 seconds each — tight. Every slide
should carry one idea. If a slide needs two ideas, split it.

**The 3×3 grid motif.** Introduce it empty on slide 2, fill it on slide 8, colour it on slide 9.
Reusing one graphic three times makes the talk cohere and saves the audience re-orienting.

**Numbers on slides.** Show F1 to four decimals in tables (that's the precision that exists),
but **speak in error ratios** — "7.4× more indel errors", not "0.9815 versus 0.9975". Nobody
retains four decimal places from a projector.

**Figures to use directly from the repo** (all already generated):
- `results/workflow_rulegraph.svg` — slide 5 (**not** `workflow_dag.svg`, which has 822 nodes)
- `results/roc_ecoli_INDEL.svg` — slide 10
- `results/qc/multiqc_report.html` — screenshot the per-base quality panel only if you have a spare slide; it's supporting material, not a result

**Colour.** Use one colour per aligner consistently across every slide and figure. The ROC
curves already distinguish nine pipelines; if the projector washes them out, consider
regenerating with only the three best-of-each-aligner curves for legibility.

**What to put in backup slides** (after the summary, for questions only):
- phiX 3×3 matrices (all 1.0000)
- The full per-pipeline TP/FP/FN table
- Truth-set indel size distribution
- The scoring-method comparison (single-run vs pre-split) — deep methodological detail, only if asked
- Read QC / mean quality profile

**Rehearsal check:** time slide 11 alone. It should take ~100 seconds. If it takes 60, you're
rushing the most valuable content in the talk.

---

## 4. Claims to avoid

Stating any of these would be indefensible under questioning:

- ❌ "BWA-MEM is the best aligner." → ✅ "BWA and minimap2 outperform Bowtie2 at this baseline; one seed, not yet a general claim."
- ❌ "The pipelines differ significantly." → Nothing here has been tested for significance. Say "differ", never "significantly".
- ❌ Presenting phiX's nine 1.0000 scores as a comparison.
- ❌ Quoting mean MAPQ across tools as a quality comparison — different scales (60 vs 42 cap).
- ❌ Presenting caller runtimes as a controlled speed benchmark — GATK is threaded, FreeBayes isn't.
- ❌ "MarkDuplicates removed duplicates" — it marked ~0%; no PCR was simulated. Say it's run for pipeline realism and is inert here.

---

## 5. Cut-list for an 8-minute slot

Drop, in this order:
1. Slide 10 (ROC) — fold its one insight into slide 8 as a sentence
2. Slide 7 (alignment metrics) — keep only the placement-accuracy row, move into slide 8
3. Slide 6 (progress table) — say it in one sentence over slide 5
4. Trim slide 11 from three traps to **two** (keep ploidy and coordinates; drop representation)

**Never cut:** slide 3 (the FROM/TO diagram), slide 9 (the headline), slide 11 (silent
failures), slide 12 (negative control). Those four are the talk.
