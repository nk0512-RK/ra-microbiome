"""
Visualisation module: all reusable plotting functions.
"""

def plot_phylum_barplot(feature_table, metadata, group_col):
    """Stacked bar chart of phylum-level composition by group."""
    pass

def plot_genus_heatmap(feature_table, metadata, top_n=20):
    """Heatmap of top N genera across samples."""
    pass

def plot_pcoa(distance_matrix, metadata, colour_col):
    """2D PCoA scatter plot coloured by metadata variable."""
    pass

def plot_volcano(results_df, lfc_col, pval_col):
    """Volcano plot: log fold change vs significance."""
    pass
