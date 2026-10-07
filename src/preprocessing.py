"""
Preprocessing module: filtering, CLR transform, taxonomy rollup.
"""

def filter_features(feature_table, min_prevalence=0.10):
    """Remove ASVs present in fewer than min_prevalence fraction of samples."""
    pass

def clr_transform(feature_table):
    """Apply centred log-ratio transform to compositional data."""
    pass

def rollup_taxonomy(feature_table, taxonomy, level="Genus"):
    """Aggregate ASV-level counts to a given taxonomic level."""
    pass
