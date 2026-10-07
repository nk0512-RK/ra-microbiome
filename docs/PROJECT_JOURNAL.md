# Project Journal

**Gut Microbiome Dysbiosis in Rheumatoid Arthritis**

A running record of decisions, findings, and things that caught us out. Newest
entries at the bottom.

---

## Contents

- [Setup](#setup)
- [Notebook 01 — Data Acquisition](#notebook-01--data-acquisition)
- [Notebook 02 — Preprocessing](#notebook-02--preprocessing)
- [Notebook 03 — EDA](#notebook-03--eda)
- [Notebook 04 — Diversity](#notebook-04--diversity)
- [Notebook 05 — Differential Abundance](#notebook-05--differential-abundance)
- [Notebook 06 — Machine Learning](#notebook-06--machine-learning)
- [Open Questions](#open-questions)

---

## Setup

**Environment.** conda env `ra-microbiome`, Python 3.10, macOS. Core stack: pandas,
numpy, scipy, scikit-learn, scikit-bio, matplotlib, seaborn, plotly, biom-format.
Added later via pip: `pyreadr` (reads R `.rds` files without needing R installed),
`shap`.

**Decision — build the full structure up front.** Created the complete folder tree,
including directories for stages not yet started, and stubbed every `src/` module.
Reasoning: when PICRUSt2 output needs a home in a month, the folder already exists.
Nothing gets reorganised later.

**Decision — central config.** All paths and parameters live in `config.py`.
Notebooks import from it and never hardcode. Adding a second cohort becomes a
one-line change rather than a find-and-replace across six notebooks.

---

## Notebook 01 — Data Acquisition

### Decision — real data over simulated

Initially drafted a synthetic dataset to move fast. Rejected: a portfolio project
built on fabricated data proves nothing.

Found the Li et al. 2025 *Scientific Data* release instead — 2,238 samples, already
QIIME2-processed, openly hosted on Figshare. Larger and better documented than the
older Scher 2013 cohort that was the first candidate.

### Findings

**Data arrives pre-rarefied.** Every sample has exactly 10,000 reads; the authors
already subsampled to even depth. Saves a preprocessing step and means samples are
directly comparable for diversity metrics.

**Orientation.** The `.rds` ASV table loads as ASVs-in-rows, samples-in-columns.
Transposed to samples-in-rows, which is what every downstream tool expects.

**Taxonomy parsing.** Taxa arrive as single Greengenes-format strings
(`k__Bacteria;p__Firmicutes;...;g__Faecalibacterium`). Split on `;`, stripped the
`k__`/`p__`/`g__` prefixes, expanded into one column per rank.

**Alignment.** Perfect 2,238-sample overlap between ASV table and metadata, no
orphans either direction. Taxonomy had 729 extra ASVs not in the feature table
(likely filtered during the authors' QC but retained in the reference); dropped to
the 39,868-ASV intersection.

### Gotcha #1 — the `object` dtype bug

`disease_status` came back entirely NaN. The metadata's `group_numeric` column held
`1` and `2`, and `.map({1: "HC", 2: "RA"})` matched nothing.

Initial diagnosis was that the values were floats (`1.0 ≠ 1` in a dict lookup).
Wrong. Printing `.dtype` showed `object` — they were **strings** (`"1"`, `"2"`). The
fix (`.astype(float).astype(int)`) worked anyway because it parses strings fine, but
the reason it worked was not the reason given.

> **Lesson.** `object` is pandas' catch-all for "not a clean number," usually
> strings. It silently breaks numeric operations rather than erroring. When
> something numeric behaves oddly, check `.dtype` first.

After the fix: 1,204 HC / 1,034 RA, zero NaN — exactly matching the published
sample counts.

### Confounder flagged

Location is heavily skewed. Beijing accounts for 1,928 of 2,238 samples (86%).
Flagged as the primary batch/site variable to control for.

---

## Notebook 02 — Preprocessing

### Decision — roll up to genus *before* filtering

Originally specced the other order. Changed because a genus can be represented by
dozens of individually rare ASVs that are collectively common; an ASV-level
prevalence filter would delete all of them and lose a widespread genus.

> **General principle.** Filter at the level you intend to analyse at, not upstream
> of it.

### Decision — drop unassigned taxa

ASVs labelled `Unknown` or `uncultured` were removed. These are real bacteria and
discarding them is a real cost, but an unnameable feature cannot be interpreted,
compared against literature, or turned into a probiotic candidate.

Recorded in limitations rather than done silently.

**Rollup result:** 39,868 ASVs → 445 genera. 5,702 ASVs dropped (~14%), but mean
reads per sample only fell from 10,000 to 9,612 (~4% loss). The unassigned ASVs were
overwhelmingly rare — many distinct sequences, little actual abundance. Cut the
unnameable tail without losing signal.

### Gotcha #2 — quotes vs variable names

`NameError: name 'Genus' is not defined`. Writing `tax[Genus]` makes Python look for
a *variable* called Genus. Column labels are literal text and need quotes:
`tax["Genus"]`.

> **Rule.** Quotes = literal text. No quotes = the name of something defined
> earlier. `tax` is a variable (no quotes); `"Genus"` is a label (quotes).

Also nearly shipped a silent bug — `"Unkown"` misspelled. Would not have errored,
just quietly matched nothing. Worse than a crash.

### Dtype check — clean

`.describe()` reported `float64`, which looked suspicious after Gotcha #1. Checked
`genus_counts.dtypes.value_counts()`: all 445 columns `int64`. The float was just
`.describe()` computing means and standard deviations, which are inherently
non-integer. Data was fine.

Worth keeping: the reflex to check rather than assume.

### Sanity check passed

Most prevalent genera after rollup: *Bacteroides* (99.4%), *Blautia* (97.9%),
*Faecalibacterium* (95.8%), *Streptococcus* (95.2%), *Escherichia-Shigella* (94.4%),
*Bifidobacterium* (93.2%).

Exactly the genera that should dominate a healthy human gut, including several known
butyrate producers. Strong evidence the taxonomy parsing and rollup are correct.

### Decision — prevalence threshold lowered 10% → 5%

At 10%: 445 → 109 genera, 336 dropped, accounting for only ~1.5% of reads (mean
9,612 → 9,467). The filter was removing genuinely near-empty features.

But 10% is the aggressive end for a dataset this large. Thresholds like 10% partly
exist to protect underpowered studies from noise, and n=2,238 is not underpowered —
at 5%, a genus still needs to appear in ~112 samples to survive, which is plenty for
a group comparison. 336 dropped genera is a lot of chances to have discarded a
rare-but-real RA-associated taxon before ever testing it.

Changed `min_prevalence` to 0.05. The 10% and 20% comparisons move to Phase 2 as a
formal sensitivity analysis (spec 2.1).

**Effect of the change:** 109 → 128 genera, only 19 recovered, worth ~0.6% of reads.
Indicates a gap in the prevalence distribution — genera are mostly either widespread
or near-absent, with little in between. Suggests the Phase 2 sensitivity analysis
will likely find results stable across thresholds.

### Outlier noted — do not drop

One sample retained only 4,572 of 10,000 reads after filtering (46%), against a
cohort mean of 9,467. Its gut is unusually dominated by bacteria that are either
unnameable or rare in this cohort — a real biological property, not a data error.

> **Explicitly not excluded.** Dropping samples for looking strange has no
> principled stopping point, and outliers in disease cohorts are often the
> interesting cases rather than the broken ones. Exclusion criteria must be set on
> outcome-blind grounds and applied before looking at results.

**Action:** check whether this sample appears as a beta-diversity outlier in the
PCoA. *(Resolved in Notebook 04 — it does not.)*

### Gotcha #3 — scikit-bio function rename

`from skbio.stats.composition import clr, multiplicative_replacement` failed with
`ImportError: cannot import name 'multiplicative_replacement'`.

The function was renamed in newer scikit-bio versions. Correct import:

```python
from skbio.stats.composition import clr, multi_replace
```

Diagnosed by listing the module's contents:
`[f for f in dir(skbio.stats.composition) if not f.startswith("_")]`

> **Lesson.** An ImportError for a function that definitely exists usually means a
> version rename, not a broken install. Listing what the module actually exposes
> settles it in one line. Worth pinning the scikit-bio version in
> `environment.yml`.

### Decision — CLR over the alternatives

| Method | Why not |
| --- | --- |
| **ALR** (additive log-ratio) | Results depend on which feature is chosen as denominator — arbitrary and hard to defend |
| **ILR** (isometric log-ratio) | Mathematically cleanest, produces genuinely independent coordinates. But coordinates are *combinations* of taxa, not taxa, so you cannot say "*Prevotella* was enriched." Fatal for our purposes |
| **PhILR** | ILR structured along a phylogenetic tree — more interpretable than plain ILR, but requires a tree we do not have |
| **Rarefaction + relative abundance** | Still compositional; all the same problems remain. Fine for plots, not for stats |
| **TSS / CSS / TMM / DESeq2** | Borrowed from RNA-seq; ongoing debate whether those assumptions transfer to sparser, more variable microbiome data |

**Chose CLR because** genus names stay interpretable, it is compositionally valid,
ALDEx2 (Phase 2) uses it internally so Phase 1 and Phase 2 stay methodologically
consistent, and it is established enough that it will not be questioned.

The trade — slight mathematical compromise (CLR's residual row-sum dependency) for
full interpretability — is the standard one in this field.

> **Known weakness for limitations.** CLR cannot distinguish a structural zero
> ("this genus is not present") from a sampling zero ("present but sequencing missed
> it"). Zero replacement treats both identically. No clean solution exists.

### Notebook 02 complete

**Outputs**

| File | Shape | Purpose |
| --- | --- | --- |
| `genus_counts.csv` | 2,238 × 128 | Filtered raw counts, for descriptive plots |
| `genus_clr.csv` | 2,238 × 128 | CLR-transformed, for all statistics and ML |

**Validation passed**

- Row sums exactly 0.0 across all 2,238 samples — a mathematical property of CLR,
  confirming the transform applied correctly
- Zero infinities, zero NaNs — zero replacement caught everything before the log
- Value range −2.74 to 9.50. Asymmetry is expected: CLR compresses the low end while
  leaving the top open, so dominant genera push well above zero

---

## Notebook 03 — EDA

### Why look before testing

Everything to this point was mechanical — load, clean, transform. EDA is where you
form expectations before running any test. The order matters: run statistics first
and you will unconsciously read the plots to fit the result you already have.

### Major finding — recruitment site is confounded with disease status

| Location | HC | RA | RA % |
| --- | --- | --- | --- |
| Beijing | 1,109 | 819 | 42.5% |
| Shanxi | 46 | 83 | 64.3% |
| Henan | 40 | 65 | 61.9% |
| Guangdong | 9 | 20 | 69.0% |
| **Hubei** | **0** | **47** | **100%** |

**Hubei is completely confounded** — 47 RA, zero HC. For those samples "from Hubei"
and "has RA" are the same statement. No statistical adjustment can separate them;
this is a logical impossibility, not a power problem.

The gradient runs the same direction everywhere: Beijing is 42.5% RA, every other
site is 62–100% RA. Since gut microbiome composition varies substantially with
geography and diet, a classifier could learn "not from Beijing → probably RA" and
achieve respectable AUC while learning nothing about bacteria.

Likely cause: controls recruited at the authors' Beijing base, RA patients pulled
from collaborating hospitals elsewhere. Standard recruitment practice with this exact
side effect.

### Decision — stratified analysis, Beijing as primary cohort

Primary analysis runs on Beijing only: 1,928 samples, 1,109 HC / 819 RA. Well
powered, reasonable group balance, internally consistent for site.

The remaining 310 non-Beijing samples become a replication set.

> **Pre-registered decision rule** (set before seeing any results, deliberately):
> if findings replicate in the non-Beijing samples, report as robust. If they do
> not, rerun on the full cohort with location as a covariate and evaluate which
> approach is better justified at that point.

**Alternatives considered**

- *Location as covariate on full cohort* — standard, retains all samples, but cannot
  rescue Hubei. A covariate adjusts for variables that overlap between groups; Hubei
  has zero overlap.
- *Drop Hubei only* — defensible on outcome-blind grounds, but does not address the
  broader site-RA gradient.

**Cost:** 310 samples. **Benefit:** a result that survives the first objection any
reviewer or domain expert will raise.

### Step 1 — Phylum composition: sanity check, no signal

Mean relative abundance (%), Beijing cohort:

| Phylum | HC | RA |
| --- | --- | --- |
| Firmicutes | 65.11 | 66.56 |
| Bacteroidetes | 23.77 | 21.86 |
| Proteobacteria | 7.51 | 7.03 |
| Actinobacteria | 2.96 | 4.05 |
| Fusobacteria | 0.29 | 0.26 |
| Tenericutes | 0.27 | 0.18 |

Firmicutes + Bacteroidetes ≈ 89%, textbook human faecal 16S. **This is the main
value of the plot** — it confirms the whole pipeline is sound end to end.

F:B ratio 2.74 (HC) → 3.05 (RA), direction consistent with published RA findings but
small. Actinobacteria is the largest proportional shift (+37% relative), notable
because it contains *Collinsella*, repeatedly implicated in RA.

**Concluded: not a finding.** Means of proportions, no dispersion, no test. More
importantly, phylum aggregation *masks* opposing genus-level shifts: Firmicutes
contains both *Faecalibacterium* (anti-inflammatory) and *Streptococcus* (often
elevated in inflammation). If one rises and the other falls, the phylum total barely
moves and the signal disappears entirely.

**→ Decision: go down to genus level.**

### Step 2 — Genus heatmap: no visible structure

Top 30 genera by mean CLR, all 1,928 Beijing samples, ordered HC-then-RA so groups
form contiguous blocks. Used CLR rather than counts (raw counts span orders of
magnitude — *Bacteroides* would saturate the colour scale).

**Result: no visible difference between the HC and RA blocks.**

This is an honest negative, not a failed plot. With 1,928 samples at a few pixels
each, a genus would have to differ dramatically *and* consistently to be visible.
Real microbiome effects are shifts in distribution, not on/off switches. The eye
cannot do this job — which is precisely why differential abundance testing exists.

Figure retained as a legitimate "no gross structure" result.

**→ Decision: the visual approach is exhausted. Go numerical.**

### Step 3 — Mean CLR differences: direction-finding

Computed mean CLR per genus per group and took the difference. **Explicitly not a
statistical test** — no p-values, no FDR correction, no effect sizes. Purpose is to
know where to look in notebook 05.

| Enriched in RA | Δ CLR | | Depleted in RA | Δ CLR |
| --- | --- | --- | --- | --- |
| *[Ruminococcus] gnavus* | +1.31 | | *Klebsiella* | −0.76 |
| *Erysipelatoclostridium* | +0.86 | | *Bilophila* | −0.59 |
| *Veillonella* | +0.84 | | *Ruminococcaceae UCG-002* | −0.58 |
| *Tyzzerella 4* | +0.81 | | *Alistipes* | −0.57 |
| *Bifidobacterium* | +0.62 | | *Parasutterella* | −0.51 |
| *Lactobacillus* | +0.53 | | *Dorea* | −0.44 |

**What this surfaced that the plots could not**

*R. gnavus* is the largest single effect by a wide margin. Well documented as
expanded in IBD; produces an inflammatory polysaccharide driving TNF-α secretion.
Emerged without being looked for.

Four oral-associated genera enriched independently — *Veillonella*, *Streptococcus*,
*Granulicatella*, *Lactobacillus*. Oral taxa colonising the gut is documented in RA
and connects to the periodontal disease / citrullination hypothesis. Four separate
genera pointing the same way is a pattern, not noise.

**Counter-observations — keep, do not quietly drop**

- *Faecalibacterium* is essentially flat (4.33 HC vs 4.29 RA) despite being the
  single most-reported RA-depleted genus in the literature.
- *Klebsiella* and *Bilophila*, both usually framed as pro-inflammatory, are
  **depleted** here. Cuts against a simple dysbiosis narrative.

### Step 4 — Abundance plot: informal significance check

Top 20 genera, RA vs HC with standard error bars. Added because the difference table
gives point estimates with no sense of spread — a large mean difference means little
if within-group variance is larger still.

Reading the error bars is an eyeball version of a significance test:

- *Bacteroides* — bars overlap almost entirely → no real difference
- *Blautia*, *Bifidobacterium* — clean separation → likely real

### Step 5 — Read depth confounder check

Retention after dropping unassigned ASVs and filtering rare genera varies by sample
(4,572–10,000). If RA samples systematically retained fewer reads, depth itself
becomes a confounder.

| Group | n | mean | std | min | median |
| --- | --- | --- | --- | --- | --- |
| HC | 1,109 | 9,502 | 423 | 4,572 | 9,609 |
| RA | 819 | 9,574 | 336 | 5,580 | 9,645 |

Mann-Whitney U: **p = 1.52e-05**

> **This is a textbook large-n significance trap.** The p-value is highly
> significant. The actual difference in medians is 36 reads out of ~9,600 — 0.4%.
>
> Had we looked only at the p-value, the conclusion would have been "read depth is a
> confounder, adjust for it." Looking at magnitude: note it, move on.
>
> This is exactly why the project's guiding principles require effect size alongside
> significance.

**Conclusion:** read depth is not a meaningful confounder.

**Outlier identified: HC0745**, 4,572 reads retained (46%). A healthy control, not
RA. Not excluded, per the earlier outcome-blind exclusion decision.

Incidental: 8 of the 10 lowest-depth samples are HC, consistent with HC showing
higher variance (std 423 vs 336). Plausibly because RA patients were recruited
through clinics with tighter criteria.

### EDA summary

Progression was: phylum (too coarse, masks opposing shifts) → genus heatmap (right
resolution, wrong tool) → numerical differences (found structure) → abundance plot
with SE (rough check on whether differences exceed variance).

Each step was chosen because the previous one hit a specific limit, not because it
was next on a list.

---

## Notebook 04 — Diversity

### Alpha diversity — three metrics, deliberately

Ran observed features, Shannon, and Simpson rather than one. They weight different
parts of the abundance distribution: observed counts all genera equally, Simpson
weights dominant taxa heavily, Shannon sits between. Agreement across all three
means a robust result; disagreement localises where the signal is.

Computed on **count** data, not CLR. These metrics expect abundances; CLR values are
centred and include negatives, which breaks the maths. The only place in the project
that returns to the count matrix.

| Metric | HC mean | RA mean | p | Cliff's d | Effect |
| --- | --- | --- | --- | --- | --- |
| Observed features | 53.9 | 47.8 | <0.0001 | −0.252 | small |
| Shannon | 2.695 | 2.597 | <0.0001 | −0.172 | small |
| Simpson | 0.859 | 0.858 | 0.0001 | −0.105 | **negligible** |

### Interpretation

**All three significant. Only one is interpretable.** Simpson at p=0.0001 with
d=−0.105 is the clearest demonstration yet of why significance alone fails at
n=1,928.

The effect-size ordering (observed > Shannon > Simpson) is itself the finding:

- Simpson weights dominant taxa and shows nothing → the abundant core community is
  unchanged in RA
- Observed counts everything equally and shows the largest effect → genera are
  missing
- Since Simpson is flat, the missing genera must be low-abundance
- Shannon lands between, exactly where it should if loss is confined to rare taxa

**Claim supported:** RA patients show modestly reduced genus richness, concentrated
in low-abundance taxa, with the dominant community largely intact. Cliff's d =
−0.252 means a random control has more genera than a random RA patient ~63% of the
time.

**Claim NOT supported:** "RA patients have depleted gut diversity." Overstates it.

This precision came from running three metrics instead of one. A Shannon-only
analysis would have given "significantly reduced diversity, small effect" and missed
where the loss actually sits.

### Beta diversity — Aitchison distance

Used Aitchison (Euclidean on CLR) rather than Bray-Curtis. Bray-Curtis is the field
convention but operates on relative abundances, inheriting every compositional
problem CLR was adopted to solve. Deliberate departure from convention.

1,928 samples → 1.86M pairwise distances. Mean 27.28, range 6.08–42.68.

### PCoA

PC1 16.1%, PC2 6.3%, first five PCs 33.7% cumulative. Typical for microbiome data.

> **Caveat for reading the plot.** The 2D view shows 22.4% of total variation. Two
> thirds of the structure is in unseen dimensions, so "the plot looks mixed" is weak
> evidence either way. PERMANOVA uses the full distance matrix and is the real test.

Visible gradient along PC1: left edge below PC2=0 heavily HC, right edge beyond
PC1=10 heavily RA, centre thoroughly mixed. A distributional shift, not two clusters.

**HC0745 resolved.** PC1 77th percentile, PC2 11th percentile. Unusual but well
inside the main cloud, not a beta-diversity outlier. Atypical in read retention and
alpha diversity but normal in overall community structure — which supports the
earlier decision to retain it.

### PERMANOVA — the headline result

pseudo-F = 17.74, **p = 0.001** (floor for 999 permutations — no random relabelling
produced a statistic this extreme), **R² = 0.0091**.

> **Both numbers are needed for an honest statement.** RA and HC communities are
> reliably distinguishable as groups, *and* disease status accounts for 0.91% of the
> variation in gut composition. The remaining 99% is diet, age, host genetics,
> medication, stool consistency, and noise.
>
> Not a contradiction — at n=1,928 a consistent small shift is detectable with high
> confidence.

Typical microbiome disease-status R² is 1–5%; studies reporting much higher usually
have a confounder doing the work. 0.91% is at the low end and consistent with the
modest alpha effects and 0.4–1.3 CLR genus differences.

### PERMDISP — ruling out the dispersion artefact

F = 2.61, **p = 0.101, non-significant.**

PERMANOVA cannot distinguish a centroid shift from a dispersion difference. This
mattered here: alpha diversity showed RA with visibly tighter distributions, making
dispersion a live alternative explanation. It is ruled out — the PERMANOVA result is
a genuine compositional shift.

Many published microbiome papers omit PERMDISP. Running it is the correct practice.

### Location comparison — partially revises an earlier framing

PERMANOVA with location as grouping, full 2,238-sample cohort:
pseudo-F = 4.63, p = 0.001, **R² = 0.0082 (0.82%)**.

> **This undercuts how the stratification decision was originally framed.** The
> implication at the time was that recruitment site might be a dominant source of
> structure swamping the disease signal. It is not — location (0.82%) and disease
> (0.91%) are comparable, with disease slightly ahead.

Caveats on direct comparison: different sample sets (2,238 vs 1,928), different group
counts (5 vs 2), and location's effect is diluted by Beijing being 86% of samples.

**Does the decision still hold? Yes, on narrower grounds.** The real argument was
always structural, not magnitude: Hubei is 47 RA / 0 HC, which no covariate
adjustment can fix. That reasoning is untouched. What changes is that the magnitude
concern is downgraded.

---

## Notebook 05 — Differential Abundance

### Method

Mann-Whitney U per genus across all 128, Benjamini-Hochberg FDR correction, Cliff's
delta for effect size. Beijing primary (n=1,928), non-Beijing replication (n=310).

**Why FDR rather than Bonferroni.** 128 tests at p<0.05 would yield ~6 false
positives by chance. Bonferroni controls the probability of *any* false positive,
overly harsh across 128 correlated features. BH controls the expected *proportion* of
false positives among significant results.

**Why effect size governs interpretation.** At n=1,928 significance is cheap. The
q-value answers "is this real"; Cliff's delta answers "does it matter." Only the
second question is interesting at this sample size.

### The headline numbers

| | Count |
| --- | --- |
| Genera tested | 128 |
| Significant (q<0.05) | 96 |
| …with non-negligible effect (\|d\| ≥ 0.147) | 46 |
| …with medium+ effect (\|d\| ≥ 0.330) | **0** |

**75% of genera clear FDR correction.** Reported bare, that looks dramatic. It is
not. Half of those hits have effects too small to interpret, and nothing in the
entire dataset reaches medium effect size. Largest is *Veillonella* at d = 0.274.

> **Conclusion: RA involves broad, shallow restructuring of the gut community rather
> than a few dramatic shifts.** Third independent line of evidence pointing the same
> way — PERMANOVA R² = 0.91%, alpha diversity effects all small or negligible, and
> now zero medium-effect genera.

### Top enriched in RA

| Genus | Δ CLR | Cliff's d |
| --- | --- | --- |
| *Veillonella* | +0.83 | 0.274 |
| *Terrisporobacter* | +0.49 | 0.268 |
| *Erysipelatoclostridium* | +0.86 | 0.257 |
| *[Ruminococcus] gnavus* | **+1.31** | 0.251 |
| *Rothia* | +0.27 | 0.246 |
| *Lachnospiraceae NC2004* | +0.34 | 0.243 |
| *Haemophilus* | +0.18 | 0.235 |
| *Tyzzerella 4* | +0.81 | 0.234 |
| *Staphylococcus* | +0.21 | 0.228 |
| *Granulicatella* | +0.39 | 0.225 |
| *Blautia* | +0.43 | 0.222 |

All q < 0.0001.

### Top depleted in RA

| Genus | Δ CLR | Cliff's d | Effect |
| --- | --- | --- | --- |
| *Bilophila* | −0.59 | −0.221 | small |
| *Alistipes* | −0.57 | −0.159 | small |
| *Lachnospiraceae UCG-010* | −0.46 | −0.150 | small |
| *Ruminococcaceae UCG-002* | −0.58 | −0.137 | negligible |
| *Parasutterella* | −0.50 | −0.127 | negligible |
| *Dorea* | −0.44 | −0.106 | negligible |
| *Klebsiella* | −0.76 | −0.103 | negligible |

### Primary finding — oral taxa translocation

Five of the top 15 enriched genera are oral-associated: ***Veillonella***,
***Rothia***, ***Haemophilus***, ***Granulicatella***, ***Staphylococcus***.

*Rothia* and *Granulicatella* are near-exclusively oral organisms. Finding them
elevated in **stool** means oral bacteria are translocating to the gut and
establishing there.

This connects directly to the periodontal disease association in RA and the
citrullination hypothesis for disease onset. Partially visible in EDA (four oral
genera on raw mean differences); formal testing found five and put *Veillonella* at
the top of the entire ranking.

> **This is the strongest and most biologically coherent result in the project.**

### Secondary — *R. gnavus*

d = 0.251, and the **largest raw CLR difference in the dataset at +1.31** — roughly
50% above the next largest. Well documented as expanded in IBD; produces an
inflammatory polysaccharide driving TNF-α secretion.

Note the dissociation between raw difference and effect size: *R. gnavus* has by far
the biggest mean shift but ranks 4th on Cliff's delta, because its within-group
variance is high. Ranking by raw difference would have been misleading.

### Asymmetry — enrichment dominates

The volcano plot makes this unmissable: red points fill the right side, blue points
are three isolated dots on the left.

Only **three** genera reach small effect on the depleted side against 43 on the
enriched side.

> **RA in this cohort is defined by what is gained, not what is lost.**

### Counter-observations — keep, do not quietly drop

**The SCFA-depletion story is not strongly supported.** Butyrate producers appear in
the depleted list but almost all at negligible effect. The expected narrative is
weaker here than the literature would predict.

***Faecalibacterium* is absent from the depleted results entirely.** The single
most-reported RA-depleted genus in the published literature, and it does not appear.
Flat in EDA, not a top hit here. Could be cohort-specific, dietary, or the published
association may be inflated by publication bias. **Report as a non-replication of a
well-known finding, do not explain away.**

***Klebsiella* and *Bilophila* are both DEPLETED in RA** despite both being textbook
pro-inflammatory organisms. *Bilophila* is in fact the single strongest depletion
signal in the dataset. Cuts directly against a simple "dysbiosis = more harmful
bacteria" narrative.

### Replication — the pre-registered rule resolves

Same pipeline on 310 non-Beijing samples (215 RA / 95 HC).

| | Count | % |
| --- | --- | --- |
| Beijing hits (q<0.05, \|d\|≥0.147) | 46 | — |
| Same direction in non-Beijing | **45** | **98%** |
| Also formally significant | **38** | 83% |

Under the null, directional agreement would be ~50%. 98% is strong. And 38 formally
significant in a set of 310 samples at 69% RA imbalance is more than expected — the
prediction going in was that power would fail. It did not.

> **The pre-registered rule set in Notebook 03 resolves: findings replicate, report
> as robust, retain the stratified design.** Spec 2.10 becomes optional.

### Unexpected — effects are LARGER in the replication cohort

| Genus | d (Beijing) | d (non-Beijing) |
| --- | --- | --- |
| *Veillonella* | 0.274 | 0.336 |
| *Erysipelatoclostridium* | 0.257 | 0.370 |
| *Granulicatella* | 0.225 | 0.348 |
| *Rothia* | 0.246 | 0.312 |
| *Terrisporobacter* | 0.268 | 0.315 |

**This is the opposite of regression to the mean.** If the Beijing hits were partly
noise, effects in the replication set should shrink. They amplify instead, across
nearly the whole top 15.

Plausible explanation: the non-Beijing cohort has more severe or less-treated RA, and
microbiome signal scales with disease severity. **Cannot be confirmed without
clinical severity data** — exactly the question spec 2.5 would answer, and this
result substantially raises its priority.

### The two non-replications

***Lachnospiraceae NC2004 group*** — d 0.243 → 0.027, q = 0.775. Effect essentially
vanishes. Genuinely Beijing-specific.

***Clostridium sensu stricto 1*** — d 0.158 → −0.002, q = 0.985. Flagged by the code
as a "direction flip," but **this framing overstates it**: −0.0017 is zero with a
sign attached, not a reversal. The effect is absent, not contradictory. State it that
way rather than letting "flip" imply the cohorts disagree.

---

## Notebook 06 — Machine Learning

### Setup

Random Forest, 500 trees, `max_features="sqrt"` (~11 of 128 genera per split, to
decorrelate trees), `min_samples_leaf=5` as an overfitting guard. Stratified 5-fold
CV on the Beijing cohort (1,928 samples, 42.5% RA — balanced enough that no class
weighting was needed).

### Result

**AUC per fold:** 0.755, 0.779, 0.814, 0.757, 0.760
**Mean: 0.773 ± 0.022.** Cross-validated overall AUC 0.772, matching the fold mean
almost exactly and confirming stability rather than a lucky split.

### Reconciling 0.77 AUC with 0.91% PERMANOVA R²

These look contradictory and are not. PERMANOVA measures what fraction of *total*
community variation disease status accounts for across all dimensions, most of which
are individual variability unrelated to disease. The classifier only needs the
specific directions that separate groups, and can stack many small signals.

Concretely: *Veillonella* alone (d=0.274) would give roughly AUC 0.59. Combining
dozens of such genera reaches 0.77.

> **That gap is the multivariate gain, quantified — and it is the entire
> justification for running ML alongside univariate testing.**

Calibration: 0.5 = chance, 0.7–0.8 = acceptable discrimination, clinical diagnostics
generally want 0.85+. **Real signal, not a usable diagnostic.**

### Feature importance — only 7/20 overlap with differential abundance

Top Gini features: *Klebsiella* (0.0284), *Bilophila* (0.0215), *Veillonella*
(0.0211), *Blautia* (0.0201), *Lachnospiraceae NC2004* (0.0173), *R. gnavus*
(0.0169).

Thirteen genera rank highly for the model but not in univariate testing, including
*Klebsiella* (Gini rank 1, but d = −0.103, **negligible** univariately).

**Two competing explanations, inseparable using Gini:**

1. **Genuine multivariate signal** — a genus uninformative alone can be useful in
   combination. This is what ML is for, and *Klebsiella* fits the pattern exactly.
2. **Gini cardinality bias** — Gini importances favour high-variance features with
   many distinct split points, independent of predictive value. Several RF-only
   genera (*Klebsiella*, *Escherichia-Shigella*, *Enterobacter*) are Proteobacteria:
   highly variable and zero-inflated, exactly the profile that inflates Gini.

Supporting evidence for explanation 2: the oral taxa cluster largely disappeared from
the Gini top 20 despite being the strongest univariate finding. These genera are
low-abundance and low-variance, precisely what Gini penalises.

**→ SHAP promoted from Phase 2 (spec 2.3) to a Phase 1 completion requirement.**

### SHAP — resolving the ambiguity

> **Answer: both explanations are partly right, which neither hypothesis predicted.**

**Gini bias confirmed for the oral taxa.** All five climb substantially under SHAP:

| Genus | Gini rank | SHAP rank | Shift | Cliff's d |
| --- | --- | --- | --- | --- |
| *Haemophilus* | 85 | 39 | +46 | 0.235 |
| *Staphylococcus* | 89 | 45 | +44 | 0.228 |
| *Rothia* | 58 | 26 | +32 | 0.246 |
| *Granulicatella* | 39 | 30 | +9 | 0.225 |
| *Veillonella* | 3 | 2 | +1 | 0.274 |

Four genera with near-identical effect sizes (0.225–0.246) were scattered between
Gini ranks 39 and 89. That scatter was the metric, not the biology.

**But *Klebsiella* holds rank 1 under SHAP too**, despite being univariately
negligible. It survives the bias-free metric, so this is genuine multivariate signal.
Same for *Subdoligranulum* (SHAP rank 8, d = 0.041) and *Lachnoclostridium* (rank 10,
d = −0.036).

**The overlap barely moves: 7/20 → 8/20.** The least predicted result. If Gini bias
were the whole story, SHAP should have pulled the ranking much closer to differential
abundance. It did not.

> **Conclusion: the univariate/multivariate disagreement is mostly real.** These
> methods genuinely find different things in this data and neither is wrong. The
> univariate analysis found the oral taxa story; the model found *Klebsiella* and a
> set of genera that only matter in combination. Reporting either alone would have
> missed half the picture — a concrete justification for having run ML, not a
> procedural one.

Caveat: the oral taxa recover substantially but still do not crack the SHAP top 20
except *Veillonella*. The model relies more on combination effects than on the
strongest univariate signal.

### Beeswarm — what the rankings could not show

***Klebsiella* is asymmetric.** High abundance spreads far left (strongly toward HC);
low abundance clusters tight and slightly right. The model uses presence as fairly
strong evidence *against* RA, absence as weak evidence for it.

> **This explains the univariate/multivariate disagreement directly.** Mann-Whitney
> averages over both behaviours and sees little; the model uses only the informative
> tail.

*Enterobacter* and *Lactococcus* show the identical shape — red tails left, blue
clustered at zero. Three Proteobacteria-type genera behaving the same way. Treat as
one pattern rather than three findings.

***Veillonella* is the cleanest signal in the dataset.** Blue left, red right, clean
gradient, high abundance → RA. The only genus where univariate and multivariate views
fully agree.

**All directions match notebook 05.** *Bilophila* depleted ✓, *R. gnavus* enriched ✓,
*Blautia* enriched ✓, *Dorea* depleted ✓. No contradictions, which matters given the
ranking disagreement.

**Most distributions are tight around zero with long tails** — individual genera
contribute little for most patients. The 0.773 AUC comes from aggregating many weak
signals, consistent with the broad-shallow conclusion from every other analysis.

### Phase 1 complete

Six notebooks: acquisition → preprocessing → EDA → diversity → differential abundance
→ ML. All on real data (Li et al. 2025, n=2,238), with the stratified design, FDR
correction, effect-size reporting and replication check carried throughout.

---

## Open Questions

- Do the richer clinical variables (DAS28, CRP, ESR, medication) exist in the other
  Figshare files? Needed for spec 2.5 disease-activity regression, whose priority was
  raised by the larger-effects-in-replication finding.
- Pin the scikit-bio version in `environment.yml` to prevent Gotcha #3 recurring on
  environment rebuild.
- Does a confounders-only baseline model (spec 2.2) approach the 0.773 AUC? Until
  this is run, the claim that the microbiome adds signal beyond demographics is
  untested.
