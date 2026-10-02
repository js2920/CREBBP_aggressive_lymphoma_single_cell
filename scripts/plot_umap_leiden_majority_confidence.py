#!/usr/bin/env python3
"""
Plot UMAP showing:
1. Leiden clusters (1.0 resolution) as background
2. Majority-voted celltype annotation for each cluster
3. Confidence encoded by size and shape
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Paths
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"
COMPOSITION_CSV = FIGDIR / "leiden_1.0_celltype_composition.csv"

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)

# Load composition data to get majority celltypes and percentages
print("Loading cluster composition data...")
comp_df = pd.read_csv(COMPOSITION_CSV)
comp_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_celltype']))
conf_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_pct'] / 100.0))

# Map to cells
print("Mapping majority celltypes and confidence to cells...")
adata.obs['majority_celltype'] = adata.obs['leiden_1.0'].astype(str).map(comp_dict)
adata.obs['annotation_confidence'] = adata.obs['leiden_1.0'].astype(str).map(conf_dict)

# Fill NaN values
adata.obs['majority_celltype'] = adata.obs['majority_celltype'].fillna('Unknown')
adata.obs['annotation_confidence'] = adata.obs['annotation_confidence'].fillna(0.0)

conf_min = adata.obs['annotation_confidence'].min()
conf_max = adata.obs['annotation_confidence'].max()
print(f"Confidence range: {conf_min:.3f} - {conf_max:.3f}")
# Calculate size range with the new formula
size_min = 3 + 22 * (conf_min ** 0.7)
size_max = 3 + 22 * (conf_max ** 0.7)
print(f"Size range: {size_min:.1f} - {size_max:.1f} pixels (emphasizing high confidence)")
print(f"Number of celltypes: {adata.obs['majority_celltype'].nunique()}")
print(f"Number of leiden clusters: {adata.obs['leiden_1.0'].nunique()}")

# Get unique celltypes and leiden clusters
celltypes = sorted(adata.obs['majority_celltype'].unique())
leiden_clusters = sorted(adata.obs['leiden_1.0'].unique())

# Color palette for celltypes
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')

celltype_colors = {}
for i, ct in enumerate(celltypes):
    celltype_colors[ct] = base_cmap(i % 20)[:3]

# Color palette for leiden clusters - use distinct colors for each
# Use a larger colormap and cycle through to ensure all clusters get different colors
try:
    leiden_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    leiden_cmap = plt.cm.get_cmap('tab20')
leiden_colors = {}
for i, cluster in enumerate(leiden_clusters):
    # Use tab20 and cycle if needed, or use Set3 for more pastel
    if len(leiden_clusters) <= 20:
        leiden_colors[str(cluster)] = leiden_cmap(i)[:3]
    else:
        # For more than 20 clusters, use Set3 and cycle
        try:
            set3_cmap = plt.colormaps['Set3']
        except (AttributeError, KeyError):
            set3_cmap = plt.cm.get_cmap('Set3')
        leiden_colors[str(cluster)] = set3_cmap(i % 12)[:3]

# Create figure with space for legend on the right
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.25], hspace=0.3)
ax = fig.add_subplot(gs[0, 0])
ax_legend = fig.add_subplot(gs[0, 1])
ax_legend.axis('off')

# Plot cells colored by leiden cluster, size/shape by confidence
print("Plotting cells with leiden cluster colors and confidence encoding...")
for cluster in leiden_clusters:
    cluster_str = str(cluster)
    mask = adata.obs['leiden_1.0'].astype(str) == cluster_str
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Get leiden cluster color
    leiden_color = leiden_colors[cluster_str]
    
    # Encode confidence by size and shape
    sizes = []
    markers_list = []
    
    for conf in confidences:
        # Size: use wider range and emphasize high confidence cells
        # Use a power function to make high confidence cells more prominent
        # conf^0.7 makes the scaling more dramatic for high values
        conf_adj = conf ** 0.7  # Emphasize high confidence
        # Wider range: 3 to 25 pixels for much better visibility
        size = 3 + 22 * conf_adj  # Range from 3 (low) to 25 (high)
        sizes.append(size)
        
        # Shape based on confidence
        if conf > 0.7:
            marker = 'o'  # Circle for high confidence
        elif conf > 0.4:
            marker = 's'  # Square for medium confidence
        else:
            marker = '^'  # Triangle for low confidence
        markers_list.append(marker)
    
    sizes = np.array(sizes)
    
    # Plot by marker type
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_list) == marker_type
        if mask_marker.sum() == 0:
            continue
        ax.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                  color=leiden_color, s=sizes[mask_marker], 
                  alpha=0.8, marker=marker_type,
                  edgecolors='white', linewidths=0.5, rasterized=True, zorder=2)

# Note: Cluster labels removed - information shown in legend instead

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title('UMAP: Leiden Clusters (1.0) with Majority-Voted Celltype Annotations\n' +
             'Colors = Leiden Clusters | Size/Shape = Confidence', 
             fontsize=13, fontweight='bold')

# Add shape legend for confidence (on the plot)
# Use larger marker sizes to reflect the new size range
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=15, label='High confidence (>70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=10, label='Medium confidence (40-70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=6, label='Low confidence (<40%)', markeredgecolor='black', linewidth=0.5),
]
shape_legend = ax.legend(handles=shape_legend_elements, loc='lower right', 
                        fontsize=9, title='Confidence Level', title_fontsize=10,
                        framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax.add_artist(shape_legend)

# Add leiden cluster legend on the right side
print("Creating leiden cluster legend...")
leiden_legend_elements = []
for cluster in leiden_clusters:
    cluster_str = str(cluster)
    color = leiden_colors[cluster_str]
    
    # Get majority celltype and confidence for this cluster
    majority_ct = comp_dict.get(cluster_str, 'Unknown')
    confidence_pct = conf_dict.get(cluster_str, 0.0) * 100
    
    # Truncate long celltype names for cleaner legend
    if len(majority_ct) > 20:
        majority_ct_display = majority_ct[:17] + "..."
    else:
        majority_ct_display = majority_ct
    
    # Create label with cluster info
    label = f"L{cluster}: {majority_ct_display} ({confidence_pct:.0f}%)"
    
    leiden_legend_elements.append(
        mpatches.Patch(facecolor=color, edgecolor='black', linewidth=1.5, label=label)
    )

ax_legend.legend(handles=leiden_legend_elements, loc='center left', 
                fontsize=7.5, title='Leiden Clusters (1.0)\n[Majority Celltype (Confidence %)]', 
                title_fontsize=9, framealpha=0.95, 
                bbox_to_anchor=(0, 0.5), handlelength=1.5, handletextpad=0.5)

# Use constrained_layout instead of tight_layout for better handling of subplots
fig.set_constrained_layout(True)

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_annotation_confidence.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

print("\nDone!")


