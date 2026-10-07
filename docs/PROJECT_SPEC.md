# Project Specification

**Gut Microbiome Dysbiosis in Rheumatoid Arthritis**

*Last updated: 7 October 2026 · Phase 1 complete*

---

## 1. Purpose

Use real 16S rRNA sequencing data to compare gut bacterial composition between
rheumatoid arthritis (RA) patients and healthy controls (HC). Identify which genera
differ, whether those differences relate to disease severity, and what biological
mechanisms might link gut bacteria to RA.

---

## 2. Dataset

**Source:** Li et al., *Scientific Data* (2025)
**Hosted:** [Figshare](https://figshare.com/articles/dataset/Data_for_publication_in_Scientific_Data/27603876/2)

| Property | Value |
| --- | --- |
| Total samples | 2,238 |
| RA patients | 1,034 |
| Healthy controls | 1,204 |
| Sequencing | 16S rRNA, V3–V4 region |
| Processing | QIIME2, pre-rarefied to 10,000 reads/sample |
| Raw ASVs | 39,868 |
| Recruitment sites | Beijing (1,928), Shanxi (129), Henan (105), Hubei (47), Guangdong (29) |

**Files used**

| File | Contents |
| --- | --- |
| `1.16S.ASV.profile.original.rds` | ASV count table |
| `1.sample.info.original.rds` | Sample metadata |
| `1.taxonomy.info.rds` | ASV → taxonomy mapping |

---

## 3. Guiding Principles

These apply across every stage and are not optional add-ons.

**3.1 Compositional awareness.** Counts sum to a fixed total per sample, so they are
proportions in disguise. CLR transform before any statistical operation. Raw counts
for descriptive plots only.

**3.2 Effect size alongside significance.** Report how much things differ, not just
whether the difference clears a p-value threshold. Cliff's delta for non-parametric
tests, log fold-change for differential abundance.

**3.3 Confounder control.** Location, stool type, and where available medication can
all masquerade as disease signal. Handle explicitly, not as an afterthought.

**3.4 Replication.** A finding present in one subgroup or under one parameter choice
is a hypothesis, not a result.

---

## 4. Phase 1 — MVP

**Goal:** a complete, working, defensible analysis end to end.

**Status: complete.**

| # | Notebook | Contents |
| --- | --- | --- |
| 01 | `01_data_acquisition.ipynb` | Load `.rds` via pyreadr, transpose ASV table, parse taxonomy strings, clean metadata, align indices, save to `data/processed/` |
| 02 | `02_preprocessing.ipynb` | Genus rollup, prevalence filtering, CLR transform |
| 03 | `03_eda.ipynb` | Phylum composition, genus heatmap, confounder visualisation, read-depth check |
| 04 | `04_diversity.ipynb` | Alpha diversity (Shannon, Simpson, observed) with Cliff's delta. Beta diversity (Aitchison, PCoA, PERMANOVA, PERMDISP) |
| 05 | `05_diff_abundance.ipynb` | Per-genus testing, FDR correction, effect sizes, volcano plot, replication check |
| 06 | `06_ml.ipynb` | Random Forest, stratified 5-fold CV, ROC/AUC, Gini importances |
| 06b | `06_ml.ipynb` (extended) | SHAP values, direction-aware attribution, resolution of the Gini/SHAP ranking disagreement |

**Phase 1 output:** a ranked list of RA-associated genera, a statistical test of
whether RA and HC communities separate, and a classifier quantifying how much RA
signal microbiome composition carries.

---

## 5. Phase 2 — Rigour and Depth

**Goal:** make the Phase 1 findings defensible and biologically meaningful.

> **Sequencing note.** Item 2.3 (SHAP) was promoted out of Phase 2 and completed as
> part of Phase 1, because the notebook 06 Random Forest produced a feature ranking
> that could not be honestly interpreted without it. Item 2.9 was also completed
> early, in notebook 05. Numbering is retained so journal cross-references stay valid.

### 2.1 Sensitivity analysis — prevalence threshold

Rerun the pipeline at 5%, 10%, and 20% prevalence filters. Compare which genera are
called significant under each. Report findings as robust (significant at all
thresholds) or threshold-dependent.

Given n=2,238, 5% is well-powered and 10% may be unnecessarily conservative.

**Status:** pending. **Estimated:** 2–3 hours.

### 2.2 Confounder baseline model

Train a classifier using only non-bacterial variables (location, stool type,
collection date). If it performs comparably to the microbiome model, the microbiome
is not adding signal beyond demographics.

This is the key check against the model learning recruitment site rather than
disease, and the first question any reader will ask about the 0.773 AUC.

**Status:** pending. **Estimated:** 1 hour. **Priority: do first.**

### 2.3 SHAP values — ✅ COMPLETE (promoted to Phase 1)

**Original rationale.** Replace raw Gini importances with SHAP. Gini is biased
toward high-cardinality features and gives no direction; SHAP gives per-sample
attribution and tells you whether high abundance of a genus pushes prediction
toward RA or away from it.

**Why it was promoted.** Notebook 06 produced only 7/20 overlap between the top
Random Forest features and the top differential abundance hits. Two explanations
were consistent with that result and could not be separated using Gini importance:

1. **Genuine multivariate signal** — genera uninformative alone but useful in
   combination. This is the stated justification for running ML at all.
2. **Gini cardinality bias** — several RF-only genera (*Klebsiella*,
   *Escherichia-Shigella*, *Enterobacter*) are Proteobacteria: high-variance and
   zero-inflated, exactly the profile that inflates Gini.

Reporting the RF ranking as settled while that ambiguity stood would have been
misleading, making SHAP a Phase 1 completion requirement rather than a refinement.

**Outcome:** both explanations proved partly right. See journal, Notebook 06.

### 2.4 XGBoost comparison

Second model alongside Random Forest. Tree-based methods consistently outperform
linear models on sparse microbiome data; comparing two guards against
model-specific artefacts.

**Status:** pending. **Estimated:** 1–2 hours.

### 2.5 Disease activity regression

If DAS28 / CRP / ESR scores can be obtained from the remaining Figshare files,
reframe as regression: predict disease severity from microbiome composition rather
than binary RA/HC classification.

Moves the headline claim from "RA patients differ" to "dysbiosis scales with
disease severity."

**Priority raised.** The notebook 05 replication found effect sizes *larger* in the
non-Beijing cohort, hinting at severity-dependent signal. This item would test that
directly.

**Status:** blocked — requires clinical variables not present in the loaded
metadata. **Next action:** audit remaining Figshare files.

### 2.6 Compositional differential abundance methods

Add ANCOM-BC and ALDEx2 alongside the Phase 1 Mann-Whitney approach. Both handle
compositionality and zero-inflation with purpose-built models. Treat a genus as a
high-confidence hit only if significant under more than one method.

Given that 96/128 genera reached significance in Phase 1, a method with different
assumptions about zeros is a meaningful robustness test rather than a formality.

**Status:** pending. **Estimated:** 4–6 hours (rpy2 interop is fiddly).

### 2.7 PICRUSt2 functional prediction

Infer metabolic pathways from taxonomy. Pathways of interest: SCFA synthesis
(butyrate, propionate), LPS biosynthesis, tryptophan metabolism, folate
biosynthesis (methotrexate interaction).

State clearly that this is *predicted*, not measured.

**Status:** pending. **Estimated:** 4–8 hours (installation is the main cost).

### 2.8 Location as covariate

Test how much variance location explains relative to disease status.

**Partially complete.** PERMANOVA with location as grouping gave R² = 0.0082 (0.82%)
versus disease status at 0.91% — comparable, not dominant. See journal, Notebook 04.

### 2.9 Non-Beijing replication — ✅ COMPLETE

Re-run the pipeline on the 310 non-Beijing samples and compare against the Beijing
primary analysis.

**Outcome:** 45/46 Beijing hits same direction, 38/46 also formally significant.
Findings replicate. See journal, Notebook 05.

### 2.10 Location-as-covariate sensitivity

Conditional on 2.9. If findings had failed to replicate, rerun on the full
2,238-sample cohort with location as a covariate.

**Status:** optional. The pre-registered rule resolved in favour of the stratified
design, so this is no longer required. Would remain informative as a robustness
check. Note that Hubei (n=47, 100% RA) cannot be adjusted for and would need
excluding under this approach.

---

## 6. Phase 3 — Extensions

### 3.1 Cross-cohort replication

Obtain a second independent RA cohort. Train on one, test on the other with no
retraining. If AUC holds, the signal generalises across populations and sequencing
batches. If it collapses, within-cohort performance was partly batch effect.

### 3.2 Multi-omics integration

Correlate microbiome with metabolomics or serology where matched data exists.
Pairwise correlation networks are achievable; MOFA+ is the stretch goal.

### 3.3 Streamlit dashboard

Interactive PCoA explorer, volcano plot, classifier ROC curve, replication table.

### 3.4 Longitudinal analysis

If multi-timepoint data becomes available: mixed-effects models with patient as
random effect, treatment response trajectories.

---

## 7. Known Limitations

To be stated explicitly in the write-up, not buried.

**Correlational, not causal.** Cannot distinguish bacteria driving disease from
bacteria responding to inflammation, medication effects, or recruitment artefacts.

**Genus-level resolution only.** 16S V3–V4 cannot reliably resolve species or
strain. Claims about *Prevotella* are supportable; claims about *Prevotella copri*
are not.

**Within-genus cancellation.** Rolling up sums ASVs within a genus. If one strain is
enriched and another depleted, they cancel and the effect disappears.

**16S under-detects rare taxa** relative to shotgun metagenomics, and over-weights
dominant organisms.

**Unassigned taxa discarded.** 5,702 ASVs (~14%) with no genus label were dropped.
These are real bacteria excluded because they cannot be interpreted.

**CLR cannot distinguish zero types.** A structural zero ("not present") and a
sampling zero ("present but missed") are treated identically by zero replacement.
No clean solution exists.

**Single cohort, single country.** All samples are Chinese, 86% from Beijing.
Generalisability to other populations is untested.

**Stratification cost.** 310 samples held out of the primary analysis. Hubei
(n=47) carries no separable information at all.

**PICRUSt2 (Phase 2) predicts function from taxonomy** rather than measuring gene
content directly.

---

## 8. Field Context

Of 42 autoimmune gut microbiome studies in one published meta-analysis, 30 used 16S
rRNA, 9 used shotgun metagenomics, and 3 used both — roughly 70% of the field uses
the approach taken here.

A pediatric ulcerative colitis study that sequenced the same samples with both
methods found 16S produced comparable results to shotgun for alpha diversity, beta
diversity, and prediction accuracy, which covers the whole of Phase 1.

Credibility rests on analytical rigour, not sequencing method.
