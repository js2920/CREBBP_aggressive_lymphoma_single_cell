#!/usr/bin/env python3
"""
Wilcoxon Rank Sum Test for Discriminating Genes per Leiden Cluster
===================================================================
Equivalent to the R script wilcoxon_rank_A.R but for AnnData objects.

Extracts top discriminating genes for each Leiden cluster (resolution 1.0)
using Wilcoxon rank sum test, generates dotplot and saves results.

Author: Generated script
Date: 2024
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

# Input AnnData object
ADATA_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad"

# Output directory
OUTPUT_DIR = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/wilcoxon_rank_sum"

# Leiden cluster column (resolution 1.0)
LEIDEN_KEY = "leiden_1.0"  # Adjust if the column name differs

# Analysis parameters
TOP_N_GENES_PER_CLUSTER = 25  # Top genes per cluster for gene panel
MAX_GENES_PANEL = 100         # Maximum genes in dotplot panel
MIN_CELLS_PER_GROUP = 2       # Minimum cells required per cluster
LOGFC_THRESHOLD = 0.0         # Log fold change threshold (0 = no filter)

# ══════════════════════════════════════════════════════════════════════════════
# 0. Setup output directory
# ══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output directory: {OUTPUT_DIR}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load AnnData object
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nLoading AnnData from: {ADATA_PATH}")
adata = sc.read_h5ad(ADATA_PATH)
print(f"Loaded object with {adata.n_obs} cells and {adata.n_vars} genes")

# Check available columns
print(f"\nAvailable obs columns: {list(adata.obs.columns)}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Identify Leiden cluster column
# ══════════════════════════════════════════════════════════════════════════════

# Try to find the leiden 1.0 resolution column
leiden_candidates = [col for col in adata.obs.columns if 'leiden' in col.lower()]
print(f"\nLeiden-related columns found: {leiden_candidates}")

# Select the appropriate column
if LEIDEN_KEY in adata.obs.columns:
    cluster_key = LEIDEN_KEY
elif 'leiden_res1.0' in adata.obs.columns:
    cluster_key = 'leiden_res1.0'
elif 'leiden' in adata.obs.columns:
    cluster_key = 'leiden'
elif len(leiden_candidates) > 0:
    # Find one with "1.0" or "1" in name, or take first
    res10_cols = [c for c in leiden_candidates if '1.0' in c or '_1' in c]
    cluster_key = res10_cols[0] if res10_cols else leiden_candidates[0]
else:
    raise ValueError("No Leiden cluster column found in adata.obs")

print(f"\nUsing cluster column: '{cluster_key}'")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Prepare data for DE analysis
# ══════════════════════════════════════════════════════════════════════════════

# Convert cluster labels to string for consistency
adata.obs['cluster_for_de'] = adata.obs[cluster_key].astype(str)

# Count cells per cluster
cluster_counts = adata.obs['cluster_for_de'].value_counts().sort_index()
print(f"\nCells per cluster:")
print(cluster_counts)

# Filter clusters with sufficient cells
valid_clusters = cluster_counts[cluster_counts >= MIN_CELLS_PER_GROUP].index.tolist()
print(f"\nClusters with >= {MIN_CELLS_PER_GROUP} cells: {len(valid_clusters)}")

if len(valid_clusters) < 2:
    raise ValueError("Fewer than two clusters have sufficient cells - DE not meaningful")

# Subset to valid clusters
adata_sub = adata[adata.obs['cluster_for_de'].isin(valid_clusters)].copy()
print(f"Subset to {adata_sub.n_obs} cells in {len(valid_clusters)} clusters")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Run Wilcoxon Rank Sum Test (each cluster vs rest)
# ══════════════════════════════════════════════════════════════════════════════

print("\nRunning Wilcoxon rank sum test...")

# Ensure we have normalized data
# Check if data looks normalized (values typically between 0-10 for log-normalized)
data_max = adata_sub.X.max() if hasattr(adata_sub.X, 'max') else np.max(adata_sub.X.toarray())
print(f"Data max value: {data_max:.2f}")

if data_max > 100:
    print("Data appears to be counts, normalizing...")
    # Store raw if not already
    if adata_sub.raw is None:
        adata_sub.raw = adata_sub.copy()
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)

# Run rank_genes_groups with Wilcoxon test
sc.tl.rank_genes_groups(
    adata_sub,
    groupby='cluster_for_de',
    method='wilcoxon',
    pts=True,  # Calculate percentage of cells expressing
    key_added='wilcoxon_de'
)

print("Differential expression analysis complete!")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Extract DE results to DataFrame
# ══════════════════════════════════════════════════════════════════════════════

print("\nExtracting DE results...")

# Get all results
de_results = []
groups = adata_sub.uns['wilcoxon_de']['names'].dtype.names

for group in groups:
    n_genes = len(adata_sub.uns['wilcoxon_de']['names'][group])
    
    group_df = pd.DataFrame({
        'gene': adata_sub.uns['wilcoxon_de']['names'][group],
        'scores': adata_sub.uns['wilcoxon_de']['scores'][group],
        'logfoldchanges': adata_sub.uns['wilcoxon_de']['logfoldchanges'][group],
        'pvals': adata_sub.uns['wilcoxon_de']['pvals'][group],
        'pvals_adj': adata_sub.uns['wilcoxon_de']['pvals_adj'][group],
        'cluster': group
    })
    
    # Add percentage expressed if available
    if 'pts' in adata_sub.uns['wilcoxon_de']:
        group_df['pct_expressed'] = adata_sub.uns['wilcoxon_de']['pts'][group]
    if 'pts_rest' in adata_sub.uns['wilcoxon_de']:
        group_df['pct_expressed_rest'] = adata_sub.uns['wilcoxon_de']['pts_rest'][group]
    
    de_results.append(group_df)

markers_df = pd.concat(de_results, ignore_index=True)

# Apply logFC threshold if specified
if LOGFC_THRESHOLD > 0:
    markers_df = markers_df[np.abs(markers_df['logfoldchanges']) >= LOGFC_THRESHOLD]

print(f"Total DE results: {len(markers_df)} gene-cluster pairs")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Build gene panel (top N per cluster)
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nSelecting top {TOP_N_GENES_PER_CLUSTER} genes per cluster...")

# Get top genes per cluster by adjusted p-value, then by log fold change
top_genes_per_cluster = (
    markers_df
    .sort_values(['cluster', 'pvals_adj', 'logfoldchanges'], 
                 ascending=[True, True, False])
    .groupby('cluster')
    .head(TOP_N_GENES_PER_CLUSTER)
)

# Get unique genes maintaining order
gene_panel = top_genes_per_cluster['gene'].drop_duplicates().tolist()

if len(gene_panel) > MAX_GENES_PANEL:
    print(f"Truncating gene panel from {len(gene_panel)} to {MAX_GENES_PANEL} genes")
    gene_panel = gene_panel[:MAX_GENES_PANEL]

print(f"Final gene panel: {len(gene_panel)} unique genes")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Save DE results
# ══════════════════════════════════════════════════════════════════════════════

# Save full results
full_results_path = os.path.join(OUTPUT_DIR, "wilcoxon_DEgenes_all_clusters.csv")
markers_df.to_csv(full_results_path, index=False)
print(f"\nFull DE results saved to: {full_results_path}")

# Save top genes per cluster
top_genes_path = os.path.join(OUTPUT_DIR, "wilcoxon_top_genes_per_cluster.csv")
top_genes_per_cluster.to_csv(top_genes_path, index=False)
print(f"Top genes per cluster saved to: {top_genes_path}")

# Save gene panel
gene_panel_path = os.path.join(OUTPUT_DIR, "gene_panel.txt")
with open(gene_panel_path, 'w') as f:
    f.write('\n'.join(gene_panel))
print(f"Gene panel saved to: {gene_panel_path}")

# ══════════════════════════════════════════════════════════════════════════════
# 8. Generate Dot Plot
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating dot plot...")

if len(gene_panel) > 0:
    # Set up figure
    n_genes = len(gene_panel)
    n_clusters = len(valid_clusters)
    
    # Calculate figure size
    fig_width = max(12, n_genes * 0.25)
    fig_height = max(6, n_clusters * 0.4)
    
    # Sort clusters numerically if possible
    try:
        sorted_clusters = sorted(valid_clusters, key=lambda x: float(x))
    except ValueError:
        sorted_clusters = sorted(valid_clusters)
    
    # Create dot plot
    sc.pl.dotplot(
        adata_sub,
        var_names=gene_panel,
        groupby='cluster_for_de',
        categories_order=sorted_clusters,
        standard_scale='var',  # Scale gene expression across clusters
        dendrogram=False,
        show=False,
        save=False
    )
    
    # Save figure
    dotplot_path = os.path.join(OUTPUT_DIR, "dotplot_top_DE_genes.pdf")
    plt.savefig(dotplot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Dot plot saved to: {dotplot_path}")
    
    # Also save as PNG
    sc.pl.dotplot(
        adata_sub,
        var_names=gene_panel,
        groupby='cluster_for_de',
        categories_order=sorted_clusters,
        standard_scale='var',
        dendrogram=False,
        show=False,
        save=False
    )
    dotplot_png_path = os.path.join(OUTPUT_DIR, "dotplot_top_DE_genes.png")
    plt.savefig(dotplot_png_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Dot plot (PNG) saved to: {dotplot_png_path}")
    
    # ══════════════════════════════════════════════════════════════════════════
    # 9. Additional visualizations
    # ══════════════════════════════════════════════════════════════════════════
    
    # Heatmap of top genes
    print("\nGenerating heatmap...")
    # Reorder categories in the obs column for proper ordering in heatmap
    adata_sub.obs['cluster_for_de'] = pd.Categorical(
        adata_sub.obs['cluster_for_de'],
        categories=sorted_clusters,
        ordered=True
    )
    sc.pl.heatmap(
        adata_sub,
        var_names=gene_panel[:50] if len(gene_panel) > 50 else gene_panel,  # Limit for readability
        groupby='cluster_for_de',
        standard_scale='var',
        show=False,
        save=False
    )
    heatmap_path = os.path.join(OUTPUT_DIR, "heatmap_top_DE_genes.pdf")
    plt.savefig(heatmap_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Heatmap saved to: {heatmap_path}")
    
    # Rank genes groups plot (scanpy style)
    print("\nGenerating rank genes groups plot...")
    sc.pl.rank_genes_groups(
        adata_sub,
        key='wilcoxon_de',
        n_genes=10,
        sharey=False,
        show=False,
        save=False
    )
    rank_plot_path = os.path.join(OUTPUT_DIR, "rank_genes_groups.pdf")
    plt.savefig(rank_plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Rank genes groups plot saved to: {rank_plot_path}")
    
    # Rank genes groups dotplot (more compact visualization)
    print("\nGenerating rank genes groups dotplot...")
    sc.pl.rank_genes_groups_dotplot(
        adata_sub,
        key='wilcoxon_de',
        n_genes=5,
        standard_scale='var',
        show=False,
        save=False
    )
    rank_dotplot_path = os.path.join(OUTPUT_DIR, "rank_genes_groups_dotplot.pdf")
    plt.savefig(rank_dotplot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Rank genes groups dotplot saved to: {rank_dotplot_path}")

else:
    print("WARNING: No genes met selection criteria - plots skipped")

# ══════════════════════════════════════════════════════════════════════════════
# 10. Summary statistics
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"Total cells analyzed: {adata_sub.n_obs}")
print(f"Total clusters: {len(valid_clusters)}")
print(f"Clusters: {', '.join(sorted_clusters)}")
print(f"Total DE gene-cluster pairs: {len(markers_df)}")
print(f"Genes in panel: {len(gene_panel)}")
print(f"\nOutput files saved to: {OUTPUT_DIR}")

# Count significant genes per cluster
sig_genes = markers_df[markers_df['pvals_adj'] < 0.05]
sig_per_cluster = sig_genes.groupby('cluster').size()
print(f"\nSignificant genes (adj. p < 0.05) per cluster:")
for cluster in sorted_clusters:
    if cluster in sig_per_cluster.index:
        print(f"  Cluster {cluster}: {sig_per_cluster[cluster]}")
    else:
        print(f"  Cluster {cluster}: 0")

print("\nAnalysis complete!")



