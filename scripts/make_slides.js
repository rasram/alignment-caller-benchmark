// Build the mid-semester presentation from docs/PRESENTATION_PLAN.md.
//
// Every number in this deck is taken from results/results.tsv,
// results/align_metrics.tsv and logs/ — see docs/PRESENTATION_PLAN.md for the
// provenance of each figure.
//
// Usage:
//   npm install pptxgenjs            # once
//   node scripts/make_slides.js      # writes docs/midsem_presentation.pptx

const path = require("path");
const fs = require("fs");
const PptxGenJS = require("pptxgenjs");

const REPO = path.resolve(__dirname, "..");
const IMG = process.env.SLIDE_IMG_DIR || path.join(REPO, "docs", "img");
const OUT = path.join(REPO, "docs", "midsem_presentation.pptx");

// ---------------------------------------------------------------------------
// Palette — deep teal / mint with a coral alert accent. Chosen so the
// "silent failures" content can be visually distinct from the results content.
// ---------------------------------------------------------------------------
const INK = "0B2E36"; // near-black teal, dark slide backgrounds
const INK2 = "12414C"; // slightly lifted, for cards on dark
const TEAL = "0E7C86";
const MINT = "3BC4A8";
const CORAL = "E06C4F";
const AMBER = "E8A33D";
const PAPER = "FFFFFF";
const MIST = "EDF4F5"; // light card fill
const SLATE = "5C7278"; // muted body text
const HEAD = "Cambria"; // safe-list serif for headers
const BODY = "Calibri"; // safe-list sans for body

const W = 13.33,
  H = 7.5;
const M = 0.62; // page margin

const pres = new PptxGenJS();
pres.layout = "LAYOUT_WIDE"; // MUST be set before any slide is added
pres.author = "DNA benchmark project";
pres.title = "Benchmarking Read Aligners and Variant Callers";

// --- small helpers ---------------------------------------------------------
// pptxgenjs mutates option objects in place, so every shadow must be a fresh
// object rather than a shared constant.
const shadow = (o = {}) =>
  Object.assign(
    { type: "outer", color: "0B2E36", blur: 10, offset: 2, angle: 90, opacity: 0.12 },
    o
  );

function darkSlide() {
  const s = pres.addSlide();
  s.background = { color: INK };
  return s;
}
function lightSlide(title, kicker) {
  const s = pres.addSlide();
  s.background = { color: PAPER };
  if (kicker) {
    s.addText(kicker.toUpperCase(), {
      x: M, y: 0.36, w: 9, h: 0.26,
      fontFace: BODY, fontSize: 12, bold: true, color: TEAL, charSpacing: 2, margin: 0,
    });
  }
  s.addText(title, {
    x: M, y: kicker ? 0.66 : 0.5, w: W - 2 * M, h: 0.7,
    fontFace: HEAD, fontSize: 34, bold: true, color: INK, margin: 0,
  });
  return s;
}
// Numbered circular badge — the repeated visual motif across the deck.
function badge(s, x, y, label, fill, textColor) {
  s.addShape(pres.ShapeType.ellipse, {
    x, y, w: 0.46, h: 0.46, fill: { color: fill }, line: { color: fill },
  });
  s.addText(String(label), {
    x, y, w: 0.46, h: 0.46, align: "center", valign: "middle",
    fontFace: BODY, fontSize: 15, bold: true, color: textColor || PAPER, margin: 0,
  });
}
function card(s, x, y, w, h, fill) {
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.08,
    fill: { color: fill || MIST }, line: { color: fill || MIST }, shadow: shadow(),
  });
}

// 3x3 pipeline grid. Reused three times: empty (slide 2), filled (slide 8).
function grid3x3(s, x, y, cellW, cellH, values, opts = {}) {
  const aligners = ["BWA-MEM", "Bowtie2", "minimap2"];
  const callers = ["GATK", "FreeBayes", "BCFtools"];
  const labelW = 1.15;
  callers.forEach((c, j) => {
    s.addText(c, {
      x: x + labelW + j * cellW, y: y - 0.34, w: cellW, h: 0.3,
      align: "center", fontFace: BODY, fontSize: 12, bold: true,
      color: opts.onDark ? MINT : TEAL, margin: 0,
    });
  });
  aligners.forEach((a, i) => {
    s.addText(a, {
      x, y: y + i * cellH, w: labelW - 0.1, h: cellH,
      align: "right", valign: "middle", fontFace: BODY, fontSize: 12, bold: true,
      color: opts.onDark ? MINT : TEAL, margin: 0,
    });
    callers.forEach((c, j) => {
      const cx = x + labelW + j * cellW;
      const v = values ? values[i][j] : null;
      const fill = v ? v.fill : opts.onDark ? INK2 : MIST;
      s.addShape(pres.ShapeType.roundRect, {
        x: cx + 0.045, y: y + i * cellH + 0.045,
        w: cellW - 0.09, h: cellH - 0.09, rectRadius: 0.06,
        fill: { color: fill }, line: { color: fill },
      });
      if (v) {
        s.addText(v.text, {
          x: cx + 0.045, y: y + i * cellH + 0.045, w: cellW - 0.09, h: cellH - 0.09,
          align: "center", valign: "middle", fontFace: BODY,
          fontSize: v.size || 14, bold: true, color: v.color || PAPER, margin: 0,
        });
      }
    });
  });
}

// Colour a 3x3 of F1 values best->worst.
function heat(vals) {
  const flat = vals.flat();
  const lo = Math.min(...flat),
    hi = Math.max(...flat);
  const ramp = (v) => {
    const t = hi === lo ? 1 : (v - lo) / (hi - lo);
    if (t > 0.82) return MINT;
    if (t > 0.55) return "8FD9C6";
    if (t > 0.3) return "F3D9A4";
    return CORAL;
  };
  return vals.map((row) =>
    row.map((v) => ({
      text: v.toFixed(4),
      fill: ramp(v),
      color: ramp(v) === MINT || ramp(v) === CORAL ? PAPER : INK,
      size: 14,
    }))
  );
}

/* =========================================================================
   SLIDE 1 — Title
   ========================================================================= */
{
  const s = darkSlide();
  s.addShape(pres.ShapeType.ellipse, {
    x: 9.7, y: -1.7, w: 5.6, h: 5.6,
    fill: { color: TEAL, transparency: 82 }, line: { color: TEAL, transparency: 82 },
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: 11.1, y: 4.2, w: 3.4, h: 3.4,
    fill: { color: MINT, transparency: 88 }, line: { color: MINT, transparency: 88 },
  });
  s.addText("MID-SEMESTER REVIEW", {
    x: M, y: 1.55, w: 8, h: 0.3,
    fontFace: BODY, fontSize: 13, bold: true, color: MINT, charSpacing: 3, margin: 0,
  });
  s.addText("Benchmarking Read Aligners\nand Variant Callers", {
    x: M, y: 2.0, w: 9.4, h: 1.9,
    fontFace: HEAD, fontSize: 44, bold: true, color: PAPER, lineSpacing: 50, margin: 0,
  });
  s.addText(
    "3 aligners × 3 callers = 9 pipelines, scored against variants we injected ourselves",
    { x: M, y: 4.05, w: 9.2, h: 0.5, fontFace: BODY, fontSize: 16, color: "BBD3D6", margin: 0 }
  );
  s.addShape(pres.ShapeType.rect, {
    x: M, y: 4.95, w: 2.2, h: 0.035, fill: { color: MINT }, line: { color: MINT },
  });
  s.addText(
    [
      { text: "Rashwanth", options: { bold: true, color: PAPER } },
      { text: "   ·   Semester 7 · DNA Sequencing Project", options: { color: "8FAEB3" } },
    ],
    { x: M, y: 5.25, w: 9, h: 0.35, fontFace: BODY, fontSize: 14, margin: 0 }
  );
  s.addNotes(
    "Name the problem and move on within 20 seconds. Do not explain the subtitle - slide 3 does that."
  );
}

/* =========================================================================
   SLIDE 2 — The question
   ========================================================================= */
{
  const s = lightSlide("Which pipeline should you actually use?", "The question");
  const bullets = [
    { t: "Finding genetic differences takes two tools in sequence", b: false },
    { t: "An ALIGNER decides where each short read came from in the genome", b: true },
    { t: "A CALLER decides whether a position genuinely differs from the reference", b: true },
    { t: "There are many of each — all widely used, all claiming to work", b: false },
  ];
  s.addText(
    bullets.map((x, i) => ({
      text: x.t,
      options: {
        bullet: x.b ? { indent: 18 } : false,
        breakLine: i !== bullets.length - 1,
        bold: !x.b,
        color: x.b ? SLATE : INK,
        paraSpaceAfter: 10,
      },
    })),
    { x: M, y: 1.85, w: 6.2, h: 2.2, fontFace: BODY, fontSize: 15, margin: 0 }
  );
  card(s, M, 4.3, 6.2, 1.35, MIST);
  s.addText("The choice is usually made by habit, or by whatever a tutorial used — not by evidence.", {
    x: M + 0.3, y: 4.5, w: 5.6, h: 0.95, fontFace: BODY, fontSize: 15, italic: true,
    color: INK, valign: "middle", margin: 0,
  });
  s.addText("9 combinations. Which one wins — and does it depend on the data?", {
    x: 7.35, y: 1.75, w: 5.35, h: 0.55, fontFace: BODY, fontSize: 14, bold: true, color: INK, margin: 0,
  });
  grid3x3(s, 7.35, 2.75, 1.55, 0.95, null);
  s.addNotes(
    "Frame as a practical decision a lab faces. Avoid tool trivia. The grid stays empty here - it gets filled on slide 8, which is why it is introduced now."
  );
}

/* =========================================================================
   SLIDE 3 — Known truth by construction (THE key diagram)
   ========================================================================= */
{
  const s = lightSlide("You cannot benchmark on real data", "The method");
  s.addText(
    "To grade a pipeline you must already know the right answer. For a real sample, nobody does — so we manufacture one where the answer is known by construction.",
    { x: M, y: 1.62, w: 11.4, h: 0.5, fontFace: BODY, fontSize: 15, color: SLATE, margin: 0 }
  );

  const boxes = [
    { x: M, label: "Reference\ngenome", fill: MIST, color: INK },
    { x: M + 3.05, label: "Mutated genome\n+ TRUTH list", fill: TEAL, color: PAPER },
    { x: M + 6.1, label: "Simulated\nreads", fill: MIST, color: INK },
    { x: M + 9.15, label: "Called\nvariants", fill: MIST, color: INK },
  ];
  boxes.forEach((b) => {
    s.addShape(pres.ShapeType.roundRect, {
      x: b.x, y: 2.5, w: 2.4, h: 1.15, rectRadius: 0.1,
      fill: { color: b.fill }, line: { color: b.fill }, shadow: shadow(),
    });
    s.addText(b.label, {
      x: b.x, y: 2.5, w: 2.4, h: 1.15, align: "center", valign: "middle",
      fontFace: BODY, fontSize: 14, bold: true, color: b.color, margin: 0,
    });
  });
  // arrows between boxes
  const arrows = [
    { x: M + 2.42, label: "inject known\nSNPs + indels", color: SLATE },
    { x: M + 5.47, label: "simulate reads\nFROM mutated", color: CORAL },
    { x: M + 8.52, label: "align TO the\nORIGINAL reference", color: CORAL },
  ];
  arrows.forEach((a) => {
    s.addShape(pres.ShapeType.rightArrow, {
      x: a.x + 0.05, y: 2.93, w: 0.55, h: 0.3,
      fill: { color: a.color }, line: { color: a.color },
    });
    // Captions must clear the boxes, which end at y=3.65. At y=3.35 they were
    // drawn on top of the box fills and became unreadable.
    s.addText(a.label, {
      x: a.x - 0.45, y: 3.78, w: 1.55, h: 0.6, align: "center",
      fontFace: BODY, fontSize: 10, bold: a.color === CORAL,
      color: a.color, margin: 0,
    });
  });
  // feedback comparison
  s.addShape(pres.ShapeType.roundRect, {
    x: M + 3.05, y: 4.68, w: 8.5, h: 0.05, rectRadius: 0,
    fill: { color: MINT }, line: { color: MINT },
  });
  s.addText("compare  →  true positives · false positives · false negatives", {
    x: M + 3.05, y: 4.81, w: 8.5, h: 0.4, align: "center",
    fontFace: BODY, fontSize: 13, bold: true, color: TEAL, margin: 0,
  });

  card(s, M, 5.5, 11.4, 1.2, "FDF0EC");
  s.addText(
    [
      { text: "The one thing to get right:  ", options: { bold: true, color: CORAL } },
      {
        text: "reads are simulated FROM the mutated genome and aligned TO the original reference. Invert it and every pipeline reports zero variants — and all nine look identical.",
        options: { color: INK },
      },
    ],
    { x: M + 0.35, y: 5.7, w: 10.7, h: 0.8, fontFace: BODY, fontSize: 14, valign: "middle", margin: 0 }
  );
  s.addNotes(
    "Spend a full minute here - this is the conceptual core. Because WE wrote the mutation list, every call is gradeable. Then plant the hook: that FROM/TO distinction is the easiest thing to get backwards, and I will come back to it on the silent-failures slide."
  );
}

/* =========================================================================
   SLIDE 4 — Experimental design
   ========================================================================= */
{
  const s = lightSlide("Experimental design", "What is being compared");
  const cells = [
    { t: "Aligners", v: "BWA-MEM\nBowtie2\nminimap2", c: TEAL },
    { t: "Callers", v: "GATK4 HaplotypeCaller\nFreeBayes\nBCFtools", c: TEAL },
    { t: "Genomes", v: "phiX174 — 5,386 bp\nE. coli K-12 — 4,641,652 bp\nboth haploid", c: MINT },
    { t: "Truth set", v: "phiX: 50 SNV + 10 indel\nE. coli: 5,000 + 1,000\ninjected by simuG", c: MINT },
    { t: "Baseline", v: "30× coverage\n150 bp paired-end\nHiSeq 2500 profile", c: AMBER },
    { t: "Scoring", v: "GA4GH standard\nrtg vcfeval\nSNV and indel separately", c: AMBER },
  ];
  cells.forEach((c, i) => {
    const col = i % 3,
      row = Math.floor(i / 3);
    const x = M + col * 3.95,
      y = 1.72 + row * 1.72;
    card(s, x, y, 3.65, 1.5, MIST);
    s.addShape(pres.ShapeType.ellipse, {
      x: x + 0.25, y: y + 0.28, w: 0.16, h: 0.16,
      fill: { color: c.c }, line: { color: c.c },
    });
    s.addText(c.t, {
      x: x + 0.52, y: y + 0.18, w: 2.9, h: 0.36,
      fontFace: BODY, fontSize: 13, bold: true, color: INK, margin: 0,
    });
    s.addText(c.v, {
      x: x + 0.52, y: y + 0.56, w: 2.95, h: 0.85,
      fontFace: BODY, fontSize: 11.5, color: SLATE, lineSpacing: 15, margin: 0,
    });
  });
  card(s, M, 5.28, 12.09, 1.42, INK);
  s.addText("HELD CONSTANT SO THE COMPARISON IS FAIR", {
    x: M + 0.35, y: 5.45, w: 6, h: 0.3,
    fontFace: BODY, fontSize: 11, bold: true, color: MINT, charSpacing: 1.5, margin: 0,
  });
  s.addText(
    "Byte-identical reads to all nine pipelines  ·  identical thread count  ·  identical filter logic  ·  no read trimming  ·  no BQSR",
    { x: M + 0.35, y: 5.82, w: 11.4, h: 0.7, fontFace: BODY, fontSize: 14, color: PAPER, margin: 0 }
  );
  s.addNotes(
    "Two genomes have different jobs: phiX is the debugging organism (runs in seconds, proves each stage), E. coli is 860x larger with real repeats so it can distinguish aligners. Explain no-BQSR in one sentence: it needs a known-variant database that does not exist for these organisms, and bootstrapping one would give GATK a step its competitors do not get. Say explicitly this is a fairness decision, not an omission."
  );
}

/* =========================================================================
   SLIDE 5 — The pipeline
   ========================================================================= */
{
  const s = lightSlide("One command, end to end", "Implementation");
  const rg = path.join(IMG, "rulegraph.png");
  if (fs.existsSync(rg)) {
    s.addImage({ path: rg, x: M, y: 1.62, w: 12.09, h: 2.58 });
  }
  const pts = [
    { n: "1", t: "Implemented twice", d: "Shell scripts first, then a Snakemake workflow — independently written" },
    { n: "2", t: "Reproduces byte-for-byte", d: "Rebuilt from scratch, the workflow reproduced the script results exactly" },
    { n: "3", t: "Sweep already wired", d: "990-run parameter sweep resolves; only the baseline has been executed" },
  ];
  pts.forEach((p, i) => {
    const x = M + i * 4.03;
    card(s, x, 4.5, 3.75, 1.9, MIST);
    badge(s, x + 0.28, 4.75, p.n, TEAL);
    s.addText(p.t, {
      x: x + 0.85, y: 4.78, w: 2.8, h: 0.35,
      fontFace: BODY, fontSize: 14, bold: true, color: INK, margin: 0,
    });
    s.addText(p.d, {
      x: x + 0.28, y: 5.3, w: 3.2, h: 0.95,
      fontFace: BODY, fontSize: 12, color: SLATE, lineSpacing: 15, margin: 0,
    });
  });
  s.addNotes(
    "Keep this short - it is orientation, not content. The one sentence worth making: two independent implementations producing identical results is the strongest evidence that neither has a silly bug."
  );
}

/* =========================================================================
   SLIDE 6 — Where the project stands
   ========================================================================= */
{
  const s = lightSlide("Where the project stands", "Progress");
  const done = [
    "Environment, tools, reference genomes",
    "Truth sets generated and verified",
    "Read simulation and quality control",
    "Alignment + alignment metrics",
  ];
  const done2 = [
    "Variant calling, all 9 pipelines",
    "GA4GH scoring, ROC analysis",
    "Reproducible workflow + sweep design",
  ];
  const todo = ["Parameter sweep executed", "Predictive model fitted"];

  const drawList = (items, x, y, mark, colr) => {
    items.forEach((t, i) => {
      s.addShape(pres.ShapeType.ellipse, {
        x, y: y + i * 0.52, w: 0.3, h: 0.3, fill: { color: colr }, line: { color: colr },
      });
      s.addText(mark, {
        x, y: y + i * 0.52, w: 0.3, h: 0.3, align: "center", valign: "middle",
        fontFace: BODY, fontSize: 13, bold: true, color: PAPER, margin: 0,
      });
      s.addText(t, {
        x: x + 0.45, y: y + i * 0.52, w: 5.1, h: 0.3, valign: "middle",
        fontFace: BODY, fontSize: 14, color: INK, margin: 0,
      });
    });
  };
  drawList(done, M, 1.72, "✓", MINT);
  drawList(done2, M + 6.1, 1.72, "✓", MINT);
  s.addText("STILL TO COME", {
    x: M + 6.1, y: 3.3, w: 5, h: 0.3,
    fontFace: BODY, fontSize: 11, bold: true, color: AMBER, charSpacing: 1.5, margin: 0,
  });
  drawList(todo, M + 6.1, 3.62, "→", AMBER);

  // Card must start below the last todo row (3.62 + 0.52 + 0.30 = 4.44).
  card(s, M, 4.72, 12.09, 1.72, INK);
  s.addText("The sweep has not run yet — but its cost is measured, not guessed", {
    x: M + 0.4, y: 4.9, w: 11.3, h: 0.35,
    fontFace: BODY, fontSize: 15, bold: true, color: PAPER, margin: 0,
  });
  const stats = [
    { v: "990", l: "pipeline runs in the full sweep" },
    { v: "~2.25 h", l: "measured compute at 8 cores" },
    { v: "~80 GB", l: "disk — the real constraint" },
  ];
  stats.forEach((st, i) => {
    const x = M + 0.4 + i * 3.9;
    s.addText(st.v, {
      x, y: 5.36, w: 3.6, h: 0.6,
      fontFace: HEAD, fontSize: 30, bold: true, color: MINT, margin: 0,
    });
    s.addText(st.l, {
      x, y: 5.97, w: 3.6, h: 0.4, fontFace: BODY, fontSize: 12, color: "9FBCC0", margin: 0,
    });
  });
  s.addNotes(
    "Be direct that the sweep has not run. Then give the measured cost - knowing your own compute budget precisely reads as competence. Note the constraint is disk, not time."
  );
}

/* =========================================================================
   SLIDE 7 — Results 1: the alignment layer
   ========================================================================= */
{
  const s = lightSlide("The aligners are effectively tied", "Results · alignment");
  const rows = [
    ["Aligner", "Mapping", "Placement ±10 bp", "Mean MAPQ", "Runtime", "Peak RAM"],
    ["BWA-MEM", "100%", "99.03%", "59.1", "5.2 s", "382 MB"],
    ["Bowtie2", "99.98%", "99.01%", "41.2", "17.6 s", "67 MB"],
    ["minimap2", "100%", "99.02%", "59.1", "1.9 s", "491 MB"],
  ];
  s.addTable(
    rows.map((r, ri) =>
      r.map((c, ci) => ({
        text: c,
        options: {
          bold: ri === 0 || ci === 0,
          color: ri === 0 ? PAPER : INK,
          fill: { color: ri === 0 ? TEAL : ri % 2 ? PAPER : MIST },
          fontSize: ri === 0 ? 12 : 13,
          align: ci === 0 ? "left" : "center",
        },
      }))
    ),
    {
      x: M, y: 1.72, w: 12.09, colW: [2.4, 1.7, 2.6, 1.9, 1.7, 1.79],
      rowH: 0.42, fontFace: BODY, valign: "middle", border: { pt: 0, color: PAPER },
    }
  );
  const takeaways = [
    { t: "A pure aligner metric", d: "The simulator records where every read really came from, so placement accuracy is measurable with no caller involved.", c: TEAL },
    { t: "All three tie at ~99%", d: "The missing 1% is repetitive sequence where 150 bp cannot resolve the location. Same information limit for all.", c: MINT },
    { t: "Speed and memory invert", d: "minimap2 is ~9× faster than Bowtie2 but uses ~7× the memory. A trade-off, not a defect.", c: AMBER },
  ];
  takeaways.forEach((tk, i) => {
    const x = M + i * 4.03;
    card(s, x, 3.72, 3.75, 2.05, MIST);
    s.addShape(pres.ShapeType.ellipse, {
      x: x + 0.28, y: 3.98, w: 0.34, h: 0.34, fill: { color: tk.c }, line: { color: tk.c },
    });
    s.addText(tk.t, {
      x: x + 0.75, y: 3.96, w: 2.9, h: 0.38,
      fontFace: BODY, fontSize: 13.5, bold: true, color: INK, margin: 0,
    });
    s.addText(tk.d, {
      x: x + 0.28, y: 4.46, w: 3.2, h: 1.2,
      fontFace: BODY, fontSize: 11.5, color: SLATE, lineSpacing: 14, margin: 0,
    });
  });
  s.addText(
    "Mean MAPQ is not comparable across tools — BWA caps at 60, Bowtie2 at 42. That difference matters later.",
    { x: M, y: 6.0, w: 12.09, h: 0.35, fontFace: BODY, fontSize: 11.5, italic: true, color: SLATE, margin: 0 }
  );
  s.addNotes(
    "Three points in order: (1) placement accuracy isolates the aligner completely - no caller involved; (2) all three tied at ~99%, so at this baseline alignment is an EASY problem; (3) runtime and memory are inverted. Do NOT compare mean MAPQ across tools as a quality score - it sets up slide 9."
  );
}

/* =========================================================================
   SLIDE 8 — The 3x3 matrices
   ========================================================================= */
{
  const s = lightSlide("Nine pipelines, scored", "Results · F1");
  const snv = [
    [0.9953, 0.9953, 0.9948],
    [0.9909, 0.9906, 0.9914],
    [0.9946, 0.9953, 0.9944],
  ];
  const ind = [
    [0.9975, 0.997, 0.997],
    [0.9935, 0.989, 0.9815],
    [0.997, 0.9975, 0.997],
  ];
  s.addText("SNV — F1", {
    x: M + 1.15, y: 1.72, w: 4.6, h: 0.35, align: "center",
    fontFace: BODY, fontSize: 15, bold: true, color: INK, margin: 0,
  });
  grid3x3(s, M, 2.48, 1.55, 0.86, heat(snv));
  s.addText("Indel — F1", {
    x: M + 7.05, y: 1.72, w: 4.6, h: 0.35, align: "center",
    fontFace: BODY, fontSize: 15, bold: true, color: INK, margin: 0,
  });
  grid3x3(s, M + 5.9, 2.48, 1.55, 0.86, heat(ind));

  card(s, M, 5.35, 5.85, 1.5, MIST);
  s.addText("Why F1 is the headline", {
    x: M + 0.3, y: 5.5, w: 5.2, h: 0.3,
    fontFace: BODY, fontSize: 13, bold: true, color: INK, margin: 0,
  });
  s.addText(
    "A caller reporting only its single most confident variant scores precision 1.0 and recall 0.0002. F1 is the harmonic mean — it correctly calls that useless.",
    { x: M + 0.3, y: 5.84, w: 5.3, h: 0.9, fontFace: BODY, fontSize: 11.5, color: SLATE, lineSpacing: 14, margin: 0 }
  );
  card(s, M + 6.24, 5.35, 5.85, 1.5, "FDF0EC");
  s.addText("phiX scored 1.0000 for all nine", {
    x: M + 6.54, y: 5.5, w: 5.2, h: 0.3,
    fontFace: BODY, fontSize: 13, bold: true, color: CORAL, margin: 0,
  });
  s.addText(
    "5,386 bp with no repetitive sequence — at 30× every variant is unambiguously recoverable. phiX is a smoke test, not a comparison.",
    { x: M + 6.54, y: 5.84, w: 5.3, h: 0.9, fontFace: BODY, fontSize: 11.5, color: SLATE, lineSpacing: 14, margin: 0 }
  );
  s.addNotes(
    "Define F1 in one line, then give the worked example. Note phiX = 1.0000 for all nine and explain it in one breath, pre-empting the 'is that a bug' question - slide 12 answers it properly. The Bowtie2 row is visibly the weak one; let the colour do the work."
  );
}

/* =========================================================================
   SLIDE 9 — The headline finding
   ========================================================================= */
{
  const s = darkSlide();
  s.addText("THE HEADLINE", {
    x: M, y: 0.62, w: 8, h: 0.3,
    fontFace: BODY, fontSize: 12, bold: true, color: MINT, charSpacing: 3, margin: 0,
  });
  s.addText("At this baseline, the aligner matters\nmore than the caller", {
    x: M, y: 1.0, w: 11.5, h: 1.4,
    fontFace: HEAD, fontSize: 36, bold: true, color: PAPER, lineSpacing: 42, margin: 0,
  });

  const big = [
    { v: "20.5×", l: "aligner effect ÷ caller effect", s: "SNV — marginal mean F1 spread" },
    { v: "2.2×", l: "aligner effect ÷ caller effect", s: "Indel — marginal mean F1 spread" },
    { v: "7.4×", l: "more indel errors", s: "worst pipeline vs best" },
  ];
  big.forEach((b, i) => {
    const x = M + i * 4.03;
    s.addShape(pres.ShapeType.roundRect, {
      x, y: 2.75, w: 3.75, h: 1.85, rectRadius: 0.1,
      fill: { color: INK2 }, line: { color: INK2 },
    });
    s.addText(b.v, {
      x: x + 0.3, y: 2.95, w: 3.2, h: 0.75,
      fontFace: HEAD, fontSize: 40, bold: true, color: MINT, margin: 0,
    });
    s.addText(b.l, {
      x: x + 0.3, y: 3.72, w: 3.2, h: 0.32,
      fontFace: BODY, fontSize: 13, bold: true, color: PAPER, margin: 0,
    });
    s.addText(b.s, {
      x: x + 0.3, y: 4.05, w: 3.2, h: 0.4,
      fontFace: BODY, fontSize: 11, color: "8FAEB3", margin: 0,
    });
  });

  s.addText("Likely mechanism — a hypothesis, not yet proven", {
    x: M, y: 4.85, w: 7.2, h: 0.32,
    fontFace: BODY, fontSize: 13, bold: true, color: AMBER, margin: 0,
  });
  s.addText(
    "Bowtie2's MAPQ caps at 42 where the others cap at 60. Callers filter on MAPQ with tool-agnostic thresholds, so the same cut-off is a stricter filter on Bowtie2. Consistent with the data: Bowtie2+GATK has the worst SNV false-negative count (90 vs 47) with zero false positives — the signature of a caller discarding evidence.",
    { x: M, y: 5.2, w: 7.2, h: 1.3, fontFace: BODY, fontSize: 12, color: "BBD3D6", lineSpacing: 16, margin: 0 }
  );
  s.addShape(pres.ShapeType.roundRect, {
    x: 8.3, y: 4.85, w: 4.42, h: 1.75, rectRadius: 0.1,
    fill: { color: "3A2420" }, line: { color: CORAL },
  });
  s.addText("Limitation, stated up front", {
    x: 8.6, y: 5.02, w: 3.9, h: 0.3,
    fontFace: BODY, fontSize: 12.5, bold: true, color: CORAL, margin: 0,
  });
  s.addText(
    "One seed, one condition. The SNV spread is 0.0047 and nothing yet shows it exceeds seed-to-seed noise. The sweep's 5 seeds give the first variance estimate.",
    { x: 8.6, y: 5.38, w: 3.85, h: 1.1, fontFace: BODY, fontSize: 11.5, color: "E8C4BB", lineSpacing: 14, margin: 0 }
  );
  s.addNotes(
    "This is your one memorable claim - spend time here. Convert F1 differences into error ratios; '7.4x more indel errors' lands where '0.9815 vs 0.9975' does not. Label the mechanism explicitly as a hypothesis. Then say the limitation OUT LOUD before anyone asks - volunteering it converts your weakest point into evidence of judgement."
  );
}

/* =========================================================================
   SLIDE 10 — Precision / recall, native chart
   ========================================================================= */
{
  const s = lightSlide("Precision is free — the difference is recall", "Results · indel detail");
  // A precision-vs-recall SCATTER was tried first and rejected: BWA-MEM and
  // minimap2 land on almost identical coordinates, so one series was drawn
  // completely underneath the other and simply vanished. A grouped bar chart
  // shows the same three-cluster structure with no occlusion — equal values
  // just produce equal bars.
  s.addChart(
    pres.ChartType.bar,
    [
      { name: "GATK",      labels: ["BWA-MEM", "Bowtie2", "minimap2"], values: [99.5, 98.7, 99.4] },
      { name: "FreeBayes", labels: ["BWA-MEM", "Bowtie2", "minimap2"], values: [99.4, 98.7, 99.5] },
      { name: "BCFtools",  labels: ["BWA-MEM", "Bowtie2", "minimap2"], values: [99.4, 98.2, 99.4] },
    ],
    {
      x: M, y: 1.8, w: 7.3, h: 4.3,
      barDir: "col", barGapWidthPct: 45,
      chartColors: [TEAL, MINT, AMBER],
      showTitle: true, title: "Indel recall (%)  —  E. coli, 30× baseline",
      titleColor: INK, titleFontFace: BODY, titleFontSize: 13,
      showValue: true, dataLabelPosition: "outEnd",
      dataLabelColor: INK, dataLabelFontFace: BODY, dataLabelFontSize: 10,
      dataLabelFormatCode: "0.0",
      valAxisMinVal: 97.5, valAxisMaxVal: 100,
      catAxisLabelColor: INK, valAxisLabelColor: SLATE,
      catAxisLabelFontFace: BODY, valAxisLabelFontFace: BODY,
      catAxisLabelFontSize: 12, valAxisLabelFontSize: 10,
      valGridLine: { color: "E4EDEE", size: 1 },
      catGridLine: { style: "none" },
      showLegend: true, legendPos: "b", legendFontSize: 11, legendFontFace: BODY,
    }
  );
  const notes = [
    { t: "Seven of nine make zero false positives", d: "There is no precision/recall trade-off to exploit — the pipelines differ almost entirely in what they miss.", c: TEAL },
    { t: "The Bowtie2 group is lower, whichever caller", d: "All three Bowtie2 bars sit below every BWA and minimap2 bar. The aligner sets the ceiling.", c: CORAL },
    { t: "Why a single F1 is not enough", d: "F1 is one threshold, and QUAL is not calibrated the same way across the three callers — so comparing single F1 values partly compares quality scales.", c: AMBER },
  ];
  notes.forEach((n, i) => {
    const y = 1.85 + i * 1.55;
    card(s, 8.2, y, 4.5, 1.35, MIST);
    s.addShape(pres.ShapeType.ellipse, {
      x: 8.45, y: y + 0.22, w: 0.16, h: 0.16, fill: { color: n.c }, line: { color: n.c },
    });
    s.addText(n.t, {
      x: 8.72, y: y + 0.12, w: 3.75, h: 0.36,
      fontFace: BODY, fontSize: 12.5, bold: true, color: INK, margin: 0,
    });
    s.addText(n.d, {
      x: 8.45, y: y + 0.5, w: 4.0, h: 0.78,
      fontFace: BODY, fontSize: 11, color: SLATE, lineSpacing: 13, margin: 0,
    });
  });
  s.addNotes(
    "The three Bowtie2 bars are all lower than every BWA/minimap2 bar - the aligner sets the ceiling regardless of caller. Mention Bowtie2+FreeBayes has the highest SNV recall of any pipeline (0.9912) but by far the most false positives (50, where six pipelines have zero) - the one aggressive pipeline, invisible in the F1 table."
  );
}

/* =========================================================================
   SLIDE 11 — THE MONEY SLIDE: three silent failures
   ========================================================================= */
{
  const s = darkSlide();
  s.addText("WHAT NEARLY WENT WRONG", {
    x: M, y: 0.5, w: 8, h: 0.3,
    fontFace: BODY, fontSize: 12, bold: true, color: CORAL, charSpacing: 3, margin: 0,
  });
  s.addText("Three ways this benchmark could have produced\nconfident, wrong numbers", {
    x: M, y: 0.85, w: 11.6, h: 1.0,
    fontFace: HEAD, fontSize: 29, bold: true, color: PAPER, lineSpacing: 34, margin: 0,
  });

  const traps = [
    {
      n: "1", title: "Ploidy",
      what: "These organisms are haploid. All three callers default to diploid.",
      show: "Nothing. No error message — and scoring cannot detect it: a caller emitting 1/1 scores a PERFECT F1 against haploid truth.",
      caught: "Grepping the genotype field directly on all 18 call sets",
    },
    {
      n: "2", title: "Coordinate systems",
      what: "The simulator's truth file uses mutated-genome coordinates; aligners report reference coordinates. Same contig name.",
      show: "Placement accuracy 8.9% instead of 99.0% — all three aligners looking catastrophically broken.",
      caught: "Converting coordinates, then self-validating against the simulator's own records",
    },
    {
      n: "3", title: "Variant representation",
      what: "FreeBayes merges nearby variants into one record: 215 AAA>TAT is two SNPs.",
      show: "Splitting SNVs from indels before scoring silently discards those records and the variants inside them.",
      caught: "Decomposing complex records first, applied identically to truth and calls",
    },
  ];
  traps.forEach((t, i) => {
    const y = 2.12 + i * 1.42;
    s.addShape(pres.ShapeType.roundRect, {
      x: M, y, w: 12.09, h: 1.26, rectRadius: 0.08,
      fill: { color: INK2 }, line: { color: INK2 },
    });
    badge(s, M + 0.28, y + 0.4, t.n, CORAL);
    s.addText(t.title, {
      x: M + 0.92, y: y + 0.16, w: 2.5, h: 0.34,
      fontFace: BODY, fontSize: 15, bold: true, color: PAPER, margin: 0,
    });
    s.addText(t.what, {
      x: M + 0.92, y: y + 0.52, w: 3.1, h: 0.66,
      fontFace: BODY, fontSize: 10.5, color: "8FAEB3", lineSpacing: 12.5, margin: 0,
    });
    s.addText("WOULD HAVE SHOWN", {
      x: M + 4.25, y: y + 0.16, w: 4.2, h: 0.26,
      fontFace: BODY, fontSize: 9, bold: true, color: CORAL, charSpacing: 1, margin: 0,
    });
    s.addText(t.show, {
      x: M + 4.25, y: y + 0.44, w: 4.3, h: 0.74,
      fontFace: BODY, fontSize: 10.5, color: PAPER, lineSpacing: 12.5, margin: 0,
    });
    s.addText("CAUGHT BY", {
      x: M + 8.75, y: y + 0.16, w: 3.2, h: 0.26,
      fontFace: BODY, fontSize: 9, bold: true, color: MINT, charSpacing: 1, margin: 0,
    });
    s.addText(t.caught, {
      x: M + 8.75, y: y + 0.44, w: 3.3, h: 0.74,
      fontFace: BODY, fontSize: 10.5, color: "BBD3D6", lineSpacing: 12.5, margin: 0,
    });
  });
  s.addText(
    "Every one of these fails silently. None produces an error message. Each produces plausible-looking numbers that are wrong.",
    { x: M, y: 6.5, w: 12.09, h: 0.5, fontFace: BODY, fontSize: 14.5, bold: true, italic: true, color: MINT, margin: 0 }
  );
  s.addNotes(
    "The intellectual core of the talk - budget the most time here, about 100 seconds. End on the unifying line: every one of these fails silently, which is why the project spends more effort on verification than on running the tools. If you must cut elsewhere to protect this slide, do."
  );
}

/* =========================================================================
   SLIDE 12 — Verification
   ========================================================================= */
{
  const s = lightSlide("A perfect score is only believable if the scoring can fail", "Verification");
  s.addText("Negative control — deliberately corrupt the call set", {
    x: M, y: 1.7, w: 6.5, h: 0.32,
    fontFace: BODY, fontSize: 14, bold: true, color: INK, margin: 0,
  });
  const nc = [
    ["Call set", "F1"],
    ["Unmodified", "1.0000"],
    ["All positions shifted +5 bp", "0.0000"],
    ["Every SNV allele changed", "0.1667"],
    ["Half the calls removed", "0.6667"],
  ];
  s.addTable(
    nc.map((r, ri) =>
      r.map((c, ci) => ({
        text: c,
        options: {
          bold: ri === 0,
          color: ri === 0 ? PAPER : ri === 1 ? INK : CORAL,
          fill: { color: ri === 0 ? INK : ri % 2 ? MIST : PAPER },
          fontSize: ri === 0 ? 12 : 13,
          align: ci === 1 ? "center" : "left",
        },
      }))
    ),
    { x: M, y: 2.1, w: 6.4, colW: [4.4, 2.0], rowH: 0.44, fontFace: BODY, valign: "middle", border: { pt: 0, color: PAPER } }
  );
  s.addText(
    "Shifting every call by 5 bp collapses F1 to zero. Removing half the calls gives recall of exactly 0.500. The scoring discriminates — so 1.0 on phiX is a real result about an easy genome.",
    { x: M, y: 4.5, w: 6.4, h: 0.95, fontFace: BODY, fontSize: 12.5, color: SLATE, lineSpacing: 16, margin: 0 }
  );

  const checks = [
    { v: "+14 / +187 bp", t: "Truth-set integrity", d: "Net length change of all injected indels matches the mutated genome's actual size difference exactly." },
    { v: "83 checks", t: "Automated audit", d: "Re-derived from the data — lengths from the index, genotypes re-grepped — not just checking files exist." },
    { v: "2 implementations", t: "Independent agreement", d: "Shell scripts and the Snakemake workflow produce identical results." },
  ];
  checks.forEach((c, i) => {
    const y = 1.7 + i * 1.62;
    card(s, 7.4, y, 5.32, 1.42, MIST);
    s.addText(c.v, {
      x: 7.7, y: y + 0.14, w: 4.7, h: 0.42,
      fontFace: HEAD, fontSize: 20, bold: true, color: TEAL, margin: 0,
    });
    s.addText(c.t, {
      x: 7.7, y: y + 0.56, w: 4.7, h: 0.28,
      fontFace: BODY, fontSize: 12.5, bold: true, color: INK, margin: 0,
    });
    s.addText(c.d, {
      x: 7.7, y: y + 0.84, w: 4.75, h: 0.5,
      fontFace: BODY, fontSize: 10.5, color: SLATE, lineSpacing: 12.5, margin: 0,
    });
  });
  s.addNotes(
    "Lead with the negative control - it answers the phiX=1.0 question directly. The length-bookkeeping check deserves 15 seconds: it ties the truth file to the exact FASTA the reads came from. Counts and spot-checks both pass even if those two disagree, and if they disagreed every number in the project would be wrong with no other symptom."
  );
}

/* =========================================================================
   SLIDE 13 — What's next
   ========================================================================= */
{
  const s = lightSlide("The sweep, then the model", "Next phase");
  // One-factor-at-a-time schematic.
  // ORDER MATTERS: the connector lines are drawn FIRST so the centre circle and
  // the axis boxes paint over their endpoints. Drawn afterwards, the lines cut
  // straight through the circle and its label.
  // Geometry is also constrained so the left box clears the 0.62" page margin
  // and the right box stops short of the text column at x = 6.55.
  const cx = 3.5, cy = 3.8;
  const axes = [
    { dx: 0, dy: -1.8, label: "Coverage", vals: "5 · 10 · 20 · 50 · 100×", c: MINT },
    { dx: -1.85, dy: 1.4, label: "Read length", vals: "75 · 100 bp", c: AMBER },
    { dx: 1.85, dy: 1.4, label: "Error rate", vals: "−2 · −5 · −10 shift", c: CORAL },
  ];
  const BW = 1.75, BH = 0.68;

  axes.forEach((a) => {
    s.addShape(pres.ShapeType.line, {
      x: Math.min(cx, cx + a.dx), y: Math.min(cy, cy + a.dy),
      w: Math.abs(a.dx) || 0.01, h: Math.abs(a.dy) || 0.01,
      line: { color: a.c, width: 2 },
      flipH: a.dx < 0, flipV: a.dy < 0,
    });
  });

  s.addShape(pres.ShapeType.ellipse, {
    x: cx - 0.72, y: cy - 0.72, w: 1.44, h: 1.44,
    fill: { color: TEAL }, line: { color: TEAL }, shadow: shadow(),
  });
  s.addText("BASELINE\n30× · 150 bp", {
    x: cx - 0.72, y: cy - 0.72, w: 1.44, h: 1.44, align: "center", valign: "middle",
    fontFace: BODY, fontSize: 10.5, bold: true, color: PAPER, margin: 0,
  });

  axes.forEach((a) => {
    const bx = cx + a.dx - BW / 2, by = cy + a.dy - BH / 2;
    s.addShape(pres.ShapeType.roundRect, {
      x: bx, y: by, w: BW, h: BH, rectRadius: 0.08,
      fill: { color: a.c }, line: { color: a.c },
    });
    s.addText(a.label, {
      x: bx, y: by + 0.06, w: BW, h: 0.28, align: "center",
      fontFace: BODY, fontSize: 12, bold: true, color: PAPER, margin: 0,
    });
    s.addText(a.vals, {
      x: bx, y: by + 0.34, w: BW, h: 0.28, align: "center",
      fontFace: BODY, fontSize: 9.5, color: PAPER, margin: 0,
    });
  });

  const right = [
    { t: "11 conditions × 5 seeds × 2 genomes × 9 pipelines", d: "990 pipeline runs — the DAG already resolves for all of them" },
    { t: "Where the differences should appear", d: "All three aligners tie at 30×. Low coverage and short reads are exactly where they should separate." },
    { t: "Known weakness of this design", d: "One-factor-at-a-time cannot detect interactions — if a pipeline only fails when coverage is low AND reads are short, this misses it. A full grid would be 72 conditions, not 11." },
  ];
  right.forEach((r, i) => {
    const y = 1.78 + i * 1.6;
    card(s, 6.55, y, 6.17, 1.4, i === 2 ? "FDF0EC" : MIST);
    s.addText(r.t, {
      x: 6.85, y: y + 0.16, w: 5.6, h: 0.34,
      fontFace: BODY, fontSize: 12.5, bold: true, color: i === 2 ? CORAL : INK, margin: 0,
    });
    s.addText(r.d, {
      x: 6.85, y: y + 0.52, w: 5.6, h: 0.78,
      fontFace: BODY, fontSize: 11, color: SLATE, lineSpacing: 13.5, margin: 0,
    });
  });
  s.addNotes(
    "Be explicit that the design is one-factor-at-a-time and name its limitation before the examiner does - it cannot detect interactions. Naming your own design's weakness is worth a lot."
  );
}

/* =========================================================================
   SLIDE 14 — Summary
   ========================================================================= */
{
  const s = darkSlide();
  s.addShape(pres.ShapeType.ellipse, {
    x: 10.4, y: 3.9, w: 4.6, h: 4.6,
    fill: { color: TEAL, transparency: 86 }, line: { color: TEAL, transparency: 86 },
  });
  s.addText("IN SUMMARY", {
    x: M, y: 0.85, w: 8, h: 0.3,
    fontFace: BODY, fontSize: 12, bold: true, color: MINT, charSpacing: 3, margin: 0,
  });
  const pts = [
    { n: "1", t: "A benchmark that can detect its own failure", d: "Truth known by construction, and verified — deliberately corrupted call sets are correctly punished." },
    { n: "2", t: "The aligner choice matters more than the caller", d: "Up to 7.4× fewer indel errors. One seed so far — the sweep confirms or refutes it." },
    { n: "3", t: "Three classes of silent failure identified and controlled", d: "Ploidy, coordinate systems, variant representation. The workflow reproduces byte-for-byte." },
  ];
  pts.forEach((p, i) => {
    const y = 1.6 + i * 1.35;
    badge(s, M, y + 0.06, p.n, MINT, INK);
    s.addText(p.t, {
      x: M + 0.75, y, w: 9.3, h: 0.42,
      fontFace: HEAD, fontSize: 20, bold: true, color: PAPER, margin: 0,
    });
    s.addText(p.d, {
      x: M + 0.75, y: y + 0.46, w: 9.3, h: 0.6,
      fontFace: BODY, fontSize: 12.5, color: "8FAEB3", lineSpacing: 15, margin: 0,
    });
  });
  s.addShape(pres.ShapeType.roundRect, {
    x: M, y: 5.72, w: 9.7, h: 0.95, rectRadius: 0.1,
    fill: { color: INK2 }, line: { color: INK2 },
  });
  s.addText(
    "The hard part of a benchmark isn't running the tools — it's making sure the numbers mean what you think they mean.",
    { x: M + 0.35, y: 5.72, w: 9.1, h: 0.95, valign: "middle",
      fontFace: HEAD, fontSize: 15, italic: true, color: MINT, margin: 0 }
  );
  s.addNotes(
    "Do not read the slide. Say the one thing you want remembered - the closing line. Then stop and take questions."
  );
}

pres.writeFile({ fileName: OUT }).then(() => console.log("wrote " + OUT));
