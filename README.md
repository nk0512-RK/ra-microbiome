# Gut Microbiome Dysbiosis in Rheumatoid Arthritis

16S rRNA analysis comparing gut bacterial composition in 1,034 RA patients against
1,204 healthy controls.

Data from Li et al., *Scientific Data* 2025
([Figshare](https://figshare.com/articles/dataset/Data_for_publication_in_Scientific_Data/27603876/2)).
QIIME2-processed ASV tables, V3–V4 region, pre-rarefied to 10,000 reads per sample.

---

## What's here

| Notebook | Contents |
| --- | --- |
| `01_data_acquisition` | Load `.rds` files, parse taxonomy, align tables |
| `02_preprocessing` | Genus rollup, prevalence filtering, CLR transform |
| `03_eda` | Composition plots, confounder checks |
| `04_diversity` | Alpha and beta diversity, PERMANOVA, PERMDISP |
| `05_diff_abundance` | Per-genus testing, FDR correction, replication |
| `06_ml` | Random Forest, SHAP |

`docs/PROJECT_JOURNAL.md` has the full record of decisions and findings, including
the ones that went the wrong way. `docs/PROJECT_SPEC.md` has the plan and what's
still outstanding.

---

## Design decision worth flagging

Recruitment site is confounded with disease status in this dataset:

| Location | HC | RA | RA % |
| --- | --- | --- | --- |
| Beijing | 1,109 | 819 | 42.5% |
| Shanxi | 46 | 83 | 64.3% |
| Henan | 40 | 65 | 61.9% |
| Guangdong | 9 | 20 | 69.0% |
| Hubei | 0 | 47 | 100% |

Hubei has no controls at all, so "from Hubei" and "has RA" are the same statement for
those 47 samples. No covariate adjustment fixes that.

Primary analysis therefore runs on Beijing only (n=1,928). The other 310 samples are
held back as a replication set. The rule for what counts as replication was set
before any results were looked at.

---

## Results

**Diversity.** RA patients have modestly fewer genera (Cliff's d = −0.252 on observed
features). Shannon shows a smaller effect, Simpson essentially none — so the loss is
in rare taxa, not the dominant community.

**Community structure.** PERMANOVA p = 0.001, R² = 0.0091. RA and HC are reliably
distinguishable as groups, and disease status explains under 1% of total variation.
Both of those are true at once. PERMDISP was non-significant, so this is a genuine
compositional shift rather than one group being more scattered.

**Differential abundance.** 96 of 128 genera clear FDR correction. 46 have a
non-negligible effect size. None reach medium. That pattern is the result: broad
shallow restructuring, not a few dramatic shifts.

**The main biological finding** is five oral-associated genera enriched in stool —
*Veillonella*, *Rothia*, *Haemophilus*, *Granulicatella*, *Staphylococcus*. *Rothia*
and *Granulicatella* are near-exclusively oral organisms, so finding them in the gut
means translocation. This lines up with the periodontal disease association in RA.

*[Ruminococcus] gnavus* had the largest raw shift (+1.31 CLR), consistent with its
reported expansion in IBD.

**Classification.** Random Forest, AUC 0.773 ± 0.022 across 5 folds. SHAP used for
feature attribution after Gini importances turned out to be biased toward
high-variance Proteobacteria.

**Replication.** 45 of 46 Beijing findings go the same direction in the held-out
cohort, 38 also reach significance there despite n=310. Effect sizes were
consistently *larger* in the replication set, which is the opposite of what noise
would do.

---

## Things that didn't fit the expected story

*Faecalibacterium* shows no depletion here (4.33 HC vs 4.29 RA) despite being the
most frequently reported RA-depleted genus in the literature.

*Klebsiella* and *Bilophila* are both depleted in RA, not enriched. Both are usually
framed as pro-inflammatory.

The SCFA-depletion narrative is weakly supported at best — butyrate producers appear
in the depleted list but almost all at negligible effect sizes.

---

## Limitations

- Correlational. Can't separate bacteria driving disease from bacteria responding to
  it, or from medication effects.
- Genus-level only. 16S V3–V4 can't reliably resolve species, so claims about
  *Prevotella* hold and claims about *Prevotella copri* don't.
- 5,702 ASVs (~14%) with no genus assignment were dropped. Real bacteria, excluded
  because they can't be interpreted.
- Single cohort, single country, 86% from one city.
- 310 samples held out of the primary analysis by the stratification.
- A confounders-only baseline model hasn't been run yet, so the claim that the
  microbiome adds signal beyond demographics is still untested.

---

## Reproducing

```bash
conda env create -f environment.yml
conda activate ra-microbiome
jupyter lab
```

Download the three `.rds` files from the Figshare link into `data/raw/`, then run the
notebooks in order. Data isn't committed — notebooks 01 and 02 regenerate everything
in `data/processed/`.

`src/` is scaffolding for refactoring the notebook code into modules. It's stubbed,
not implemented — the working code lives in the notebooks.

---

## Still to do

Confounders-only baseline model, prevalence threshold sensitivity analysis,
ANCOM-BC/ALDEx2 cross-check, PICRUSt2 functional prediction. Disease activity
regression is blocked on clinical variables that aren't in the metadata file used
here.

Full list in `docs/PROJECT_SPEC.md`.
