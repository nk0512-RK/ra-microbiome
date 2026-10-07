# Project Journal

Running notes. Decisions, things that broke, things that didn't go the way I expected.
Newest at the bottom.

---

## Setup

conda env `ra-microbiome`, Python 3.10, macOS. pandas / numpy / scipy / sklearn /
scikit-bio / matplotlib / seaborn. Added `pyreadr` later for the .rds files, `shap`
much later.

Built the full folder tree up front including dirs for stages I haven't started, and
stubbed all the `src/` modules. Reasoning: when PICRUSt2 output needs somewhere to go
in a month, the folder's already there.

All paths and params in `config.py`. Notebooks import from it, never hardcode.

---

## 01 — Data acquisition

Started drafting a synthetic dataset to move fast. Scrapped it — a portfolio project
on made-up data proves nothing.

Found Li et al. 2025 (*Scientific Data*) instead. 2,238 samples, already
QIIME2-processed, on Figshare. Bigger and better documented than Scher 2013, which
was the first candidate.

Data comes pre-rarefied — every sample exactly 10,000 reads. The authors already
subsampled to even depth, so that's one preprocessing step I don't need.

ASV table loads ASVs-in-rows. Transposed.

Taxonomy arrives as Greengenes strings (`k__Bacteria;p__Firmicutes;...`). Split on
`;`, stripped the prefixes, one column per rank.

Alignment was clean — 2,238 samples in both tables, no orphans. Taxonomy had 729
extra ASVs not in the feature table, probably filtered during the authors' QC but
left in the reference. Dropped to the intersection.

### The dtype bug

`disease_status` came back all NaN. `group_numeric` held 1 and 2,
`.map({1: "HC", 2: "RA"})` matched nothing.

Thought they were floats (1.0 ≠ 1 in a dict lookup). Wrong. `.dtype` said `object` —
they were strings. `.astype(float).astype(int)` fixed it but not for the reason I
assumed.

Lesson: `object` is pandas' bin for anything that isn't a clean number. Breaks
numeric ops silently instead of erroring. Check `.dtype` first when something numeric
misbehaves.

After the fix: 1,204 HC / 1,034 RA. Matches the paper.

### Confounder flagged

Beijing is 1,928 of 2,238 samples. 86%. Noting it now, dealing with it later.

---

## 02 — Preprocessing

### Rollup before filtering, not after

Originally specced it the other way. Changed my mind: a genus can be dozens of
individually rare ASVs that are collectively common, and an ASV-level filter would
delete all of them and lose the genus entirely.

General rule — filter at the level you're going to analyse at, not upstream.

### Dropping unassigned taxa

Removed ASVs labelled `Unknown` or `uncultured`. These are real bacteria and
discarding them is a real cost, but I can't interpret a feature I can't name, can't
compare it to literature, can't do anything with it.

Goes in limitations rather than happening silently.

39,868 ASVs → 445 genera. 5,702 dropped (~14%) but mean reads only fell 10,000 →
9,612. So the unassigned ones were almost all rare. Good trade.

### Quotes

`NameError: name 'Genus' is not defined` — wrote `tax[Genus]` instead of
`tax["Genus"]`. Quotes = literal text, no quotes = variable name.

Also nearly shipped `"Unkown"` misspelled. Wouldn't have errored, just matched
nothing. Worse than a crash.

### Dtype check

`.describe()` said float64 which looked suspicious after the earlier bug. Checked
`genus_counts.dtypes.value_counts()` — all 445 columns int64. `.describe()` just
returns floats because means and SDs aren't integers. Fine.

Keeping the habit of checking though.

### Sanity check

Most prevalent genera after rollup: *Bacteroides* 99.4%, *Blautia* 97.9%,
*Faecalibacterium* 95.8%, *Streptococcus* 95.2%, *Escherichia-Shigella* 94.4%,
*Bifidobacterium* 93.2%.

That's what a healthy human gut should look like. Pipeline's working.

### Threshold: 10% → 5%

At 10%: 445 → 109 genera, 336 dropped, worth only ~1.5% of reads. So the filter was
cutting genuinely near-empty stuff.

But 10% is aggressive for a dataset this size. Thresholds like that partly exist to
protect underpowered studies and n=2,238 isn't underpowered. At 5% a genus still
needs ~112 samples to survive. And 336 dropped genera is a lot of chances to have
binned a rare-but-real RA taxon before testing it.

Changed to 0.05. 10% and 20% comparisons go to Phase 2 as a sensitivity analysis.

Effect of the change: 109 → 128 genera. Only 19 recovered, worth 0.6% of reads. So
there's a gap in the prevalence distribution — genera are mostly either widespread or
near-absent. Suggests the sensitivity analysis will find results are stable, which
would be a good thing to report.

### Outlier — not dropping it

One sample kept only 4,572 of 10,000 reads (46%) vs a cohort mean of 9,467. Its gut
is dominated by bacteria that are either unnameable or rare here. That's a real
biological property, not a data error.

Not excluding it. Dropping samples for looking weird has no stopping point, and
outliers in disease cohorts are often the interesting ones. Exclusion criteria need
to be outcome-blind and set before you look at results.

Will check whether it shows up as a beta diversity outlier later.

### skbio rename

`from skbio.stats.composition import clr, multiplicative_replacement` → ImportError.
Function was renamed. It's `multi_replace` in this version.

Found it with `[f for f in dir(skbio.stats.composition) if not f.startswith("_")]`.

An ImportError for a function that definitely exists usually means a rename, not a
broken install. Should pin the skbio version in environment.yml so this doesn't come
back.

### Why CLR and not something else

Looked at the alternatives:

- **ALR** — depends on which feature you pick as denominator. Arbitrary.
- **ILR** — mathematically the cleanest, genuinely independent coordinates. But the
  coordinates are combinations of taxa rather than taxa, so you can't say
  "*Prevotella* was enriched." Kills it for me.
- **PhILR** — ILR along a phylogenetic tree, more interpretable. Needs a tree I don't
  have.
- **Rarefaction + relative abundance** — still compositional, same problems. Fine for
  plots.
- **TSS/CSS/TMM/DESeq2** — from RNA-seq, and there's ongoing argument about whether
  those assumptions transfer to sparser microbiome data.

Went with CLR. Genus names stay interpretable, it's compositionally valid, and
ALDEx2 uses it internally so Phase 1 and Phase 2 stay consistent.

Weakness for the limitations section: CLR can't tell a structural zero ("not there")
from a sampling zero ("there, missed it"). Zero replacement treats both the same.
There's no clean fix for this.

### Done

`genus_counts.csv` and `genus_clr.csv`, both 2,238 × 128.

Row sums all exactly 0.0, which is a mathematical property of CLR so that confirms
it applied correctly. No infs, no NaNs. Range −2.74 to 9.50 — asymmetric because CLR
compresses the low end and leaves the top open.

---

## 03 — EDA

Look before testing. If I run the stats first I'll read the plots to fit whatever I
already found.

### Recruitment site is confounded with disease

| Location | HC | RA | RA % |
| --- | --- | --- | --- |
| Beijing | 1,109 | 819 | 42.5% |
| Shanxi | 46 | 83 | 64.3% |
| Henan | 40 | 65 | 61.9% |
| Guangdong | 9 | 20 | 69.0% |
| Hubei | 0 | 47 | 100% |

Hubei has zero controls. For those 47 samples "from Hubei" and "has RA" are the same
statement. That's not a power problem, it's a logical one — no covariate fixes it.

And the gradient runs the same way everywhere. Beijing 42.5% RA, everywhere else
62–100%. Gut composition varies a lot with geography and diet, so a classifier could
learn "not Beijing → probably RA" and get a decent AUC without learning anything
about bacteria.

Probably just how they recruited: controls at their Beijing base, patients from
collaborating hospitals elsewhere. Normal practice, this exact side effect.

### Decision — stratify

Primary analysis on Beijing only. 1,928 samples, 1,109 HC / 819 RA. Decent balance,
internally consistent for site.

Other 310 become a replication set.

Rule, set now before I've seen any results: if findings replicate in the non-Beijing
samples, report as robust. If not, rerun on the full cohort with location as a
covariate and decide which is better justified then.

Alternatives I considered. Location-as-covariate on everything keeps all the samples
but can't rescue Hubei, since a covariate adjusts for something that overlaps between
groups and Hubei has zero overlap. Dropping only Hubei is defensible but doesn't
touch the broader gradient.

Costs 310 samples. Worth it.

### Phylum composition

| Phylum | HC | RA |
| --- | --- | --- |
| Firmicutes | 65.11 | 66.56 |
| Bacteroidetes | 23.77 | 21.86 |
| Proteobacteria | 7.51 | 7.03 |
| Actinobacteria | 2.96 | 4.05 |
| Fusobacteria | 0.29 | 0.26 |
| Tenericutes | 0.27 | 0.18 |

F+B ≈ 89%, which is textbook. Main value of this plot is confirming the pipeline
works end to end.

F:B ratio 2.74 → 3.05, right direction per the literature but small. Actinobacteria
is the biggest proportional shift (+37% relative), interesting because it contains
*Collinsella* which keeps showing up in RA papers.

Not a finding though. These are means of proportions, no dispersion, no test. And
phylum aggregation hides opposing shifts — Firmicutes has both *Faecalibacterium*
(anti-inflammatory) and *Streptococcus* (up in inflammation). One rises, one falls,
phylum total barely moves.

Going to genus level.

### Genus heatmap

Top 30 by mean CLR, all 1,928 Beijing samples, ordered HC then RA so the groups form
blocks. Used CLR not counts — raw counts span orders of magnitude and *Bacteroides*
would saturate the colour scale.

Can't see any difference between the blocks.

That's honest rather than a failed plot. 1,928 samples at a few pixels each, a genus
would have to differ dramatically *and* consistently to show up. Real microbiome
effects are distributional shifts, not on/off. The eye can't do this, which is the
entire reason differential abundance testing exists.

Keeping the figure as a "no gross structure" result.

### Mean CLR differences

Not a test. No p-values, no correction, no effect sizes. Just want to know where to
look in notebook 05.

Enriched in RA: *R. gnavus* +1.31, *Erysipelatoclostridium* +0.86, *Veillonella*
+0.84, *Tyzzerella 4* +0.81, *Bifidobacterium* +0.62, *Lactobacillus* +0.53,
*Blautia* +0.43.

Depleted: *Klebsiella* −0.76, *Bilophila* −0.59, *Ruminococcaceae UCG-002* −0.58,
*Alistipes* −0.57, *Parasutterella* −0.51, *Dorea* −0.44.

*R. gnavus* is way out in front, 50% above the next one. It's well documented as
expanded in IBD and it makes an inflammatory polysaccharide that drives TNF-α. Came
out without me looking for it.

Four oral genera enriched independently — *Veillonella*, *Streptococcus*,
*Granulicatella*, *Lactobacillus*. Oral taxa in the gut is a documented thing in RA
and connects to the periodontal/citrullination hypothesis. Four separate genera
pointing the same way isn't noise.

Depleted list is mostly Ruminococcaceae and Lachnospiraceae, which is where the
butyrate producers live. Consistent with SCFA depletion but diffuse.

Two things that don't fit, keeping them:

*Faecalibacterium* is flat. 4.33 HC vs 4.29 RA. It's the most-reported RA-depleted
genus in the literature and it's doing nothing here.

*Klebsiella* and *Bilophila* are both depleted, and both usually get framed as
pro-inflammatory.

### Abundance plot

Top 20 with SE bars. Added this because the difference table gives point estimates
with no sense of spread, and a big mean difference means nothing if within-group
variance is bigger.

Eyeballing the error bars: *Bacteroides* overlaps almost completely, *Blautia* and
*Bifidobacterium* separate cleanly.

### Read depth check

Retention varies 4,572–10,000 after filtering. If RA systematically kept fewer reads,
depth becomes a confounder.

| Group | n | mean | std | min | median |
| --- | --- | --- | --- | --- | --- |
| HC | 1,109 | 9,502 | 423 | 4,572 | 9,609 |
| RA | 819 | 9,574 | 336 | 5,580 | 9,645 |

Mann-Whitney p = 1.52e-05.

Which is a trap. Highly significant, and the actual difference in medians is 36 reads
out of ~9,600. 0.4%.

If I'd only looked at the p-value I'd have concluded depth was a confounder and
started adjusting for it. Looking at magnitude: note it, move on. This is exactly why
the project rules say effect size alongside significance.

Outlier is HC0745, 4,572 reads. A healthy control, not RA. Not excluding it per the
earlier decision. Will check it in the PCoA.

Side note: 8 of the 10 lowest-depth samples are HC, and HC has higher variance
overall (423 vs 336). Probably because RA patients came through clinics with tighter
criteria and the controls are a looser group.

### EDA summary

Phylum (too coarse) → heatmap (right level, wrong tool) → numbers (found structure) →
abundance plot with SE (does the difference beat the variance).

Each step because the previous one hit a wall, not because it was next on a list.

---

## 04 — Diversity

### Alpha

Ran three metrics rather than one. Observed counts all genera equally, Simpson weights
the dominant ones, Shannon is in between. If they agree it's robust; if they disagree
it tells you where the signal is.

On counts, not CLR — these metrics need abundances and CLR has negatives. Only place
in the project that goes back to the count matrix.

| Metric | HC | RA | p | Cliff's d | Effect |
| --- | --- | --- | --- | --- | --- |
| Observed | 53.9 | 47.8 | <0.0001 | −0.252 | small |
| Shannon | 2.695 | 2.597 | <0.0001 | −0.172 | small |
| Simpson | 0.859 | 0.858 | 0.0001 | −0.105 | negligible |

All three significant, one interpretable. Simpson at p=0.0001 and d=−0.105 is the
clearest example yet of significance being worthless at this n.

The ordering is the actual finding. Simpson weights dominant taxa and shows nothing →
the abundant core is unchanged. Observed counts everything and shows the biggest
effect → genera are missing. Since Simpson's flat, the missing ones must be
low-abundance. Shannon in between, exactly where it should be.

So: modestly reduced richness, concentrated in rare taxa, dominant community intact.
d = −0.252 means a random control has more genera than a random RA patient about 63%
of the time.

Not "RA patients have depleted gut diversity." That overstates it.

Three metrics instead of one is what got me the precision. Shannon alone would have
given "significantly reduced, small effect" and missed where the loss sits.

### Beta

Aitchison (Euclidean on CLR) rather than Bray-Curtis. Bray-Curtis is what most papers
use but it works on relative abundances, which brings back every problem CLR was
meant to solve.

1.86M pairwise distances. Mean 27.28, range 6.08–42.68.

### PCoA

PC1 16.1%, PC2 6.3%, first five 33.7% cumulative. Normal for microbiome data.

Worth remembering when reading the plot — I'm looking at 22.4% of the variation. Two
thirds of the structure is in dimensions I can't see, so "the plot looks mixed" isn't
strong evidence either way. PERMANOVA uses the whole distance matrix.

There is a gradient along PC1. Left edge below PC2=0 is mostly HC, right edge past
PC1=10 is mostly RA, middle is mixed. Distributional shift, not clusters.

HC0745 sits at the 77th percentile on PC1, 11th on PC2. Unusual but inside the main
cloud. So it's odd in read retention and alpha diversity but normal in overall
community structure. Supports keeping it.

### PERMANOVA

pseudo-F 17.74, p = 0.001 (floor for 999 permutations), R² = 0.0091.

Need both numbers. RA and HC are reliably distinguishable as groups, and disease
status explains 0.91% of the variation. The other 99% is diet, age, genetics,
medication, stool consistency, noise. Not contradictory — at n=1,928 you can detect
a consistent small shift with high confidence.

Typical microbiome disease R² is 1–5%. Studies reporting much higher usually have a
confounder doing the work. 0.91% is low end, consistent with everything else here.

### PERMDISP

F = 2.61, p = 0.101. Not significant.

Ran this because PERMANOVA can't tell a centroid shift from a dispersion difference,
and alpha diversity showed RA with tighter distributions so dispersion was a real
alternative. It's ruled out. The PERMANOVA result is a genuine compositional shift.

A lot of published microbiome papers skip this.

### Location comparison — revises what I said earlier

PERMANOVA with location as grouping, full cohort: pseudo-F 4.63, p = 0.001,
R² = 0.0082.

This partly undercuts how I framed the stratification decision. I implied site might
be a dominant source of structure swamping the disease signal. It isn't. Location
0.82%, disease 0.91%. Comparable, disease slightly ahead.

Caveats on comparing them directly: different sample sets, different group counts,
and location's effect is diluted by Beijing being 86%.

Does the decision still hold? Yes but on narrower grounds. The real argument was
always structural — Hubei is 47/0 and no covariate fixes that. That's untouched. The
magnitude concern is downgraded.

---

## 05 — Differential abundance

Mann-Whitney per genus across all 128, BH correction, Cliff's delta. Beijing primary,
non-Beijing replication.

FDR rather than Bonferroni because 128 tests at p<0.05 gives ~6 false positives by
chance, and Bonferroni controls the probability of *any* false positive which is too
harsh for 128 correlated features. BH controls the expected proportion among the
significant ones.

Effect size governs interpretation. Already seen a 0.4% difference hit p=1.5e-05 in
this project. q-value says "is it real", Cliff's delta says "does it matter".

### Numbers

128 tested. 96 significant at q<0.05. 46 with non-negligible effect. **Zero** with
medium or above.

75% of genera clearing FDR looks dramatic if you report it bare. It isn't. Half have
effects too small to interpret and nothing in the whole dataset reaches medium.
Biggest is *Veillonella* at 0.274.

So: broad shallow restructuring, not a few dramatic shifts. Third independent line of
evidence saying the same thing after the PERMANOVA R² and the alpha diversity.

### Enriched

| Genus | Δ CLR | d |
| --- | --- | --- |
| *Veillonella* | +0.83 | 0.274 |
| *Terrisporobacter* | +0.49 | 0.268 |
| *Erysipelatoclostridium* | +0.86 | 0.257 |
| *R. gnavus* | +1.31 | 0.251 |
| *Rothia* | +0.27 | 0.246 |
| *Lachnospiraceae NC2004* | +0.34 | 0.243 |
| *Haemophilus* | +0.18 | 0.235 |
| *Tyzzerella 4* | +0.81 | 0.234 |
| *Staphylococcus* | +0.21 | 0.228 |
| *Granulicatella* | +0.39 | 0.225 |
| *Blautia* | +0.43 | 0.222 |

### Depleted

| Genus | Δ CLR | d | Effect |
| --- | --- | --- | --- |
| *Bilophila* | −0.59 | −0.221 | small |
| *Alistipes* | −0.57 | −0.159 | small |
| *Lachnospiraceae UCG-010* | −0.46 | −0.150 | small |
| *Ruminococcaceae UCG-002* | −0.58 | −0.137 | negligible |
| *Parasutterella* | −0.50 | −0.127 | negligible |
| *Dorea* | −0.44 | −0.106 | negligible |
| *Klebsiella* | −0.76 | −0.103 | negligible |

### Oral taxa

Five of the top 15 enriched are oral-associated: *Veillonella*, *Rothia*,
*Haemophilus*, *Granulicatella*, *Staphylococcus*.

*Rothia* and *Granulicatella* are basically only found in the mouth. Finding them
elevated in stool means translocation to the gut.

Connects to the periodontal disease association in RA and the citrullination
hypothesis. I'd partially spotted this in EDA on raw means — formal testing found
five and put *Veillonella* at the top of the whole ranking.

Strongest and most coherent result in the project.

### R. gnavus

d = 0.251 but the biggest raw shift in the dataset at +1.31, about 50% above the next.

Worth noting the dissociation: biggest mean difference, 4th on effect size, because
its within-group variance is high. Ranking by raw difference would have been
misleading.

### Asymmetry

Three genera reach small effect on the depleted side. 43 on the enriched side. The
volcano plot makes it obvious — red fills the right, blue is three dots on the left.

RA here is defined by what's gained, not what's lost.

### Things that don't fit

SCFA depletion isn't well supported. Butyrate producers are in the depleted list but
almost all at negligible effect. Weaker than the literature would predict.

*Faecalibacterium* doesn't appear at all. Most-reported RA-depleted genus in the
field, absent from my results. Could be cohort-specific, could be diet, could be that
the published association is inflated by publication bias. Reporting it as a
non-replication rather than explaining it away.

*Klebsiella* and *Bilophila* both depleted despite both being textbook
pro-inflammatory. *Bilophila* is the strongest depletion signal in the dataset. Cuts
against "dysbiosis = more harmful bacteria."

### Replication

Same pipeline on the 310 non-Beijing samples (215 RA / 95 HC).

46 Beijing hits. 45 same direction. 38 also formally significant.

98% directional agreement where the null would give ~50%. And 38 reaching
significance in 310 badly imbalanced samples is better than I expected — I'd assumed
power would fail.

Pre-registered rule resolves: findings replicate, keep the stratified design.

### Effects are bigger in the replication cohort

| Genus | d Beijing | d non-Beijing |
| --- | --- | --- |
| *Veillonella* | 0.274 | 0.336 |
| *Erysipelatoclostridium* | 0.257 | 0.370 |
| *Granulicatella* | 0.225 | 0.348 |
| *Rothia* | 0.246 | 0.312 |
| *Terrisporobacter* | 0.268 | 0.315 |

Opposite of regression to the mean. If the Beijing hits were partly noise these
should shrink. They're amplifying across nearly the whole top 15.

Best guess is the non-Beijing cohort has more severe or less-treated RA and the
signal scales with severity. Can't confirm without clinical data. This is exactly
what spec 2.5 would test and it bumps that item up the list.

### Two non-replications

*Lachnospiraceae NC2004* — 0.243 → 0.027, q = 0.775. Effect vanishes. Beijing-specific.

*Clostridium sensu stricto 1* — 0.158 → −0.002, q = 0.985. My code flagged this as a
"direction flip" but that's overstating it. −0.0017 is zero with a sign on it. The
effect is absent, not reversed. Should say it that way.

---

## 06 — Machine learning

RF, 500 trees, `max_features="sqrt"` (~11 genera per split), `min_samples_leaf=5`.
Stratified 5-fold on Beijing. 42.5% RA so no class weighting needed.

### Result

Folds: 0.755, 0.779, 0.814, 0.757, 0.760. Mean 0.773 ± 0.022. Cross-validated overall
0.772, basically identical to the fold mean so it's stable.

### 0.77 AUC vs 0.91% R²

Looks contradictory, isn't. PERMANOVA measures what fraction of *total* variation
disease explains across all dimensions, most of which is individual variability. The
classifier only needs the directions that separate groups and can stack many small
signals.

*Veillonella* alone would give maybe 0.59. Dozens of genera together gets 0.77. That
gap is the multivariate gain, and it's the whole reason to run ML alongside
univariate tests.

For scale: 0.5 chance, 0.7–0.8 acceptable discrimination, clinical diagnostics want
0.85+. Real signal, not a diagnostic.

### Only 7/20 overlap with differential abundance

Top Gini: *Klebsiella* 0.0284, *Bilophila* 0.0215, *Veillonella* 0.0211, *Blautia*
0.0201, *Lachnospiraceae NC2004* 0.0173, *R. gnavus* 0.0169.

13 genera rank high for the model but not univariately. *Klebsiella* is rank 1 with
d = −0.103, which is negligible.

Two explanations and Gini can't separate them.

One, it's real multivariate signal — a genus uninformative alone being useful in
combination. That's what ML is for and *Klebsiella* fits.

Two, Gini bias. Gini favours high-variance features with lots of split points
regardless of predictive value. Several of the RF-only genera (*Klebsiella*,
*Escherichia-Shigella*, *Enterobacter*) are Proteobacteria — variable and
zero-inflated, exactly what inflates Gini.

Evidence for the second: the oral taxa almost disappeared from the Gini top 20
despite being my strongest univariate finding. They're low-abundance and low-variance,
which is what Gini penalises.

So SHAP moves from Phase 2 to a Phase 1 requirement. Can't honestly report this
ranking without it.

### SHAP

Both explanations turn out to be partly right, which I wasn't expecting.

Gini bias confirmed for the oral taxa:

| Genus | Gini | SHAP | Shift | d |
| --- | --- | --- | --- | --- |
| *Haemophilus* | 85 | 39 | +46 | 0.235 |
| *Staphylococcus* | 89 | 45 | +44 | 0.228 |
| *Rothia* | 58 | 26 | +32 | 0.246 |
| *Granulicatella* | 39 | 30 | +9 | 0.225 |
| *Veillonella* | 3 | 2 | +1 | 0.274 |

Four genera with nearly identical effect sizes were scattered between Gini ranks 39
and 89. That spread was the metric, not biology.

But *Klebsiella* holds rank 1 under SHAP too. Survives the unbiased metric despite
being univariately negligible, so that's genuine multivariate signal. Same for
*Subdoligranulum* (rank 8, d = 0.041) and *Lachnoclostridium* (rank 10, d = −0.036).

And the overlap barely moves — 7/20 to 8/20. Least expected result of the lot. If
Gini bias were the whole story SHAP should have pulled the rankings much closer
together.

So the univariate/multivariate disagreement is mostly real. The methods find
different things and neither is wrong. Univariate found the oral taxa, the model
found *Klebsiella* and a set of genera that only matter in combination. Reporting
either alone misses half of it.

Caveat: the oral taxa recover but still don't crack the SHAP top 20 except
*Veillonella*. The model leans more on combination effects than on my strongest
univariate signal.

### Beeswarm

*Klebsiella* is asymmetric. High abundance spreads far left (toward HC), low
abundance clusters tight and slightly right. So the model treats presence as fairly
strong evidence against RA and absence as weak evidence for it.

That explains the disagreement directly. Mann-Whitney averages over both behaviours
and sees nothing; the model only uses the informative tail.

*Enterobacter* and *Lactococcus* have the same shape. Three Proteobacteria-type
genera behaving identically — one pattern, not three findings.

*Veillonella* is the cleanest thing in the plot. Blue left, red right, clean gradient.
Only genus where univariate and multivariate fully agree.

All directions match notebook 05 — *Bilophila* depleted, *R. gnavus* enriched,
*Blautia* enriched, *Dorea* depleted. No contradictions, which matters given how much
the rankings disagreed.

Most distributions are tight around zero with long tails. Individual genera do very
little for most patients. The 0.773 comes from aggregating weak signals, which is
consistent with everything else in this project.

---

## Open

- Are DAS28/CRP/ESR in the other Figshare files? Needed for the disease activity
  regression, which the larger-effects-in-replication result made more interesting.
- Pin skbio version in environment.yml.
- Haven't run the confounders-only baseline yet. Until I do, "the microbiome adds
  signal beyond demographics" is an assumption, not a result.
