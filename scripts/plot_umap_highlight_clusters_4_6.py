#!/usr/bin/env python3
"""
Plot UMAP highlighting leiden clusters 4 and 6
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"

# Clusters to highlight (as strings to match leiden_1.0 format)
HIGHLIGHT_CLUSTERS = ['4', '6']

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Check if leiden_1.0 exists
if 'leiden_1.0' not in adata.obs.columns:
    print("ERROR: leiden_1.0 not found in obs columns!")
    print("Available columns:", sorted(adata.obs.columns))
    exit(1)

# Check if UMAP exists
if 'X_umap' not in adata.obsm:
    print("ERROR: X_umap not found in obsm!")
    print("Available obsm keys:", list(adata.obsm.keys()))
    exit(1)

print(f"\nLeiden clusters present: {sorted(adata.obs['leiden_1.0'].unique())}")
print(f"Highlighting clusters: {HIGHLIGHT_CLUSTERS}")

# Check if highlight clusters exist
for cluster in HIGHLIGHT_CLUSTERS:
    # leiden_1.0 might be stored as string or int, so convert to match
    cluster_val = str(cluster) if isinstance(adata.obs['leiden_1.0'].iloc[0], str) else int(cluster)
    n_cells = (adata.obs['leiden_1.0'].astype(str) == str(cluster)).sum()
    print(f"  Cluster {cluster}: {n_cells:,} cells")

# Create figure
fig, ax = plt.subplots(figsize=(14, 12))

# Get all leiden clusters
all_clusters = sorted(adata.obs['leiden_1.0'].unique())
n_clusters = len(all_clusters)

# Color palette for all clusters (use muted colors for non-highlighted)
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')

# Assign colors: highlighted clusters get bright colors, others get muted
cluster_colors = {}
for cluster in all_clusters:
    cluster_str = str(cluster)
    if cluster_str in HIGHLIGHT_CLUSTERS:
        # Bright, saturated colors for highlighted clusters
        if cluster_str == '4':
            cluster_colors[cluster_str] = (1.0, 0.0, 0.0)  # Bright red
        elif cluster_str == '6':
            cluster_colors[cluster_str] = (0.0, 0.0, 1.0)  # Bright blue
    else:
        # Muted gray for non-highlighted clusters
        cluster_colors[cluster_str] = (0.7, 0.7, 0.7)  # Light gray

# Plot all clusters
print("\nPlotting UMAP...")
for cluster in all_clusters:
    cluster_str = str(cluster)
    mask = adata.obs['leiden_1.0'].astype(str) == cluster_str
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    color = cluster_colors.get(cluster_str, (0.7, 0.7, 0.7))
    
    # Highlighted clusters: larger size, full opacity, with edge
    # Non-highlighted: smaller size, lower opacity, no edge
    if cluster_str in HIGHLIGHT_CLUSTERS:
        ax.scatter(coords[:, 0], coords[:, 1], 
                  c=[color], s=50, alpha=1.0, 
                  marker='o', edgecolors='black', linewidths=1.5,
                  rasterized=True, label=f'Cluster {cluster_str} (highlighted)', zorder=3)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], 
                  c=[color], s=10, alpha=0.3, 
                  marker='o', edgecolors='none',
                  rasterized=True, zorder=1)

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title(f'UMAP: Highlighted Leiden Clusters {HIGHLIGHT_CLUSTERS}\n' +
             f'(All other clusters shown in gray)', 
             fontsize=14, fontweight='bold')

# Add legend
ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

plt.tight_layout()

# Save
output_path = FIGDIR / f"umap_leiden_1.0_highlight_clusters_{'_'.join(map(str, HIGHLIGHT_CLUSTERS))}.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

# Also save to the mouse_human_integration directory if it exists
integration_figdir = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/figures_individual_samples")
if integration_figdir.exists():
    integration_output = integration_figdir / f"umap_all_disease_state_highlight_clusters_{'_'.join(map(str, HIGHLIGHT_CLUSTERS))}.png"
    fig.savefig(integration_output, dpi=300, bbox_inches='tight')
    print(f"Also saved → {integration_output}")

plt.close()

print("\nDone!")


