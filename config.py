"""
Central configuration for ra-microbiome project.
Edit paths and parameters here — nowhere else.
"""
import os

# ── Cohorts ────────────────────────────────────────────────────────────────
# To add a new cohort later: add an entry here and set active=True.
COHORTS = {
    "scher_2013": {
        "feature_table": "data/processed/scher_feature_table.csv",
        "metadata":      "data/processed/scher_metadata.csv",
        "active":        True
    },
    # "zhang_2015": {
    #     "feature_table": "data/processed/zhang_feature_table.csv",
    #     "metadata":      "data/processed/zhang_metadata.csv",
    #     "active":        False
    # }
}

# ── Analysis parameters ────────────────────────────────────────────────────
PARAMS = {
    "min_prevalence":  0.05,   # drop ASVs present in <5% of samples
    "taxonomy_level":  "Genus",
    "random_state":    42,
    "n_cv_folds":      5,
    "fdr_method":      "fdr_bh"
}

# ── Output paths ───────────────────────────────────────────────────────────
PATHS = {
    "figures": "results/figures/",
    "tables":  "results/tables/"
}

# ── Ensure output directories exist ───────────────────────────────────────
for path in PATHS.values():
    os.makedirs(path, exist_ok=True)
