# Experimental Design Simulation for Gene Expression Analysis

**Genetics Assistant Internship — Week 2** | Rabi Narayan Patra | 27 September 2026

A complete, executable experimental design for a gene expression study — **a controlled transcriptomic analysis of hypoxia and reoxygenation injury in human endothelial cells**, with power-based sample sizing, verified randomisation and orthogonal qRT-PCR validation.

## The design in one line

Four arms × six biological replicates (24 RNA-seq libraries) testing three pre-specified contrasts — **B vs A** (hypoxia effect), **C vs B** (reoxygenation effect), **D vs B** (HIF1A dependence) — with four categories of control, a computed sample size, and five independent validation lines.

| Arm | Condition | Purpose |
|-----|-----------|---------|
| A | Normoxia 21% O₂, 24 h, control siRNA | Baseline control |
| B | Hypoxia 1% O₂, 24 h, control siRNA | Hypoxia effect |
| C | Hypoxia 1% O₂ 24 h → reoxygenation 21% O₂ 6 h | Reoxygenation effect |
| D | Hypoxia 1% O₂, 24 h, HIF1A siRNA | Mechanism arm |

## Computed design decisions (not asserted)

**Sample size** — negative-binomial power analysis (Hart et al. 2013 formula), reproduced by `scripts/02_power_and_randomisation.py`:

| Biological CV | n for 2-fold (80% power) | n for 1.5-fold | n for 4-fold | Power at n=6 (2-fold) | Smallest detectable FC at n=6 |
|---|---|---|---|---|---|
| 0.20 | 2.9 | 8.6 | 0.7 | 0.979 | 1.63 |
| 0.30 | 4.6 | 13.4 | 1.1 | 0.894 | 1.83 |
| 0.40 | 6.9 | 20.1 | 1.7 | 0.745 | 2.10 |
| 0.60 | 13.4 | 39.1 | 3.3 | 0.466 | 2.82 |

n = 6 per arm matches the benchmark recommendation of Schurch et al. (2016) *and* is justified quantitatively here. The design is honest about its limit: effects below ~1.8-fold are not reliably detectable.

**Randomisation** — 24 samples balanced across 3 processing batches, 2 passages and 2 operators, verified programmatically to be unconfounded (every arm contributes exactly 2 samples per batch, 3 per passage, 3 per operator). Output: `results/randomisation_plan.csv`.

## Expected results (falsifiable)

- **B vs A:** ~300–1,500 genes at FDR < 0.05, |log2FC| ≥ 1; VEGFA, SLC2A1, LDHA, PGK1, BNIP3 up 2–8×
- **HIF1A mRNA falls 20–80% during hypoxia** while HIF-1α protein accumulates — the counter-intuitive prediction a naive design would misread as failure
- **C vs B:** partial reversal of the hypoxic programme plus inflammatory induction (IL6, CXCL8, ICAM1, VCAM1)
- **D vs B:** ≥50% of hypoxia genes lose the response (HIF1A-dependent subset)

## Validation (five independent lines)

Orthogonal qRT-PCR (Pearson r ≥ 0.85, ≥90% directional agreement; 3 reference genes selected by geNorm M < 0.5) · independent replication (≥70% of candidates reproduce) · dose–response · mechanism/rescue · protein-level check. A gene is reported as validated only on the statistical threshold **plus** two validation lines.

## Repository contents

```
├── README.md
├── scripts/
│   └── 02_power_and_randomisation.py   # power analysis + balanced randomisation (computes both)
└── results/
    ├── power_analysis.csv              # required n across CVs and effect sizes
    ├── power_summary.json
    └── randomisation_plan.csv          # the 24-sample allocation
```

The full report, figures and reproducibility bundle are attached to the [v1.0 release](https://github.com/rabinarayanpatra350-biotech/genetics-week2-experimental-design/releases/tag/v1.0).

## Literature anchors

Schurch et al. 2016 (replication benchmark) · Hart et al. 2013 (sample size) · Love et al. 2014 (DESeq2) · Dobin et al. 2013 (STAR) · Conesa et al. 2016 (best practices) · Bustin et al. 2009 & MIQE 2.0 2025 (qPCR standards) · Vandesompele et al. 2002 (geNorm) · Livak & Schmittgen 2001 (ΔΔCq) · Semenza 2012 (hypoxia biology) · King & Wilson 1975, Brawand et al. 2011 (expression evolution) · Stark et al. 2019 (RNA-seq review)

## Author

**Rabi Narayan Patra** — Genetics Assistant Internship (2026). GitHub: [@rabinarayanpatra350-biotech](https://github.com/rabinarayanpatra350-biotech)
