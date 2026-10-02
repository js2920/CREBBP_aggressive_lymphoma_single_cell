#!/usr/bin/env python3
"""
Plot UMAP in separate panels for each condition:
1. WT_B_cells
2. Crebbp_B_cells
3. Pre_malignant
4. Malignant (combined Malignant + Matched_malignant)

Each panel shows:
- Leiden clusters (1.0 resolution) colored by cluster
- Confidence encoded by size and shape
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

# Define conditions to plot - using exact condition names from data
conditions_to_plot = {
    'WT_B_cells': ['WT_B_cells'],
    'Crebbp_B_cells': ['Crebbp_B_cells'],
    'Pre_malignant': ['Pre_malignant'],
    'Malignant (combined)': ['Malignant', 'Matched_malignant']  # Combine both malignant types
}

# Get unique leiden clusters (same across all conditions)
leiden_clusters = sorted(adata.obs['leiden_1.0'].unique())

# Color palette for leiden clusters - use distinct colors for each
try:
    leiden_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    leiden_cmap = plt.cm.get_cmap('tab20')
leiden_colors = {}
for i, cluster in enumerate(leiden_clusters):
    if len(leiden_clusters) <= 20:
        leiden_colors[str(cluster)] = leiden_cmap(i)[:3]
    else:
        try:
            set3_cmap = plt.colormaps['Set3']
        except (AttributeError, KeyError):
            set3_cmap = plt.cm.get_cmap('Set3')
        leiden_colors[str(cluster)] = set3_cmap(i % 12)[:3]

# Create figure with 2x2 grid for conditions + space for legend
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.3], hspace=0.3, wspace=0.3)

axes = []
for i, (cond_name, cond_values) in enumerate(conditions_to_plot.items()):
    row = i // 2
    col = i % 2
    ax = fig.add_subplot(gs[row, col])
    axes.append((ax, cond_name, cond_values))

# Plot each condition in its own panel
print("\nPlotting conditions in separate panels...")
for ax, cond_name, cond_values in axes:
    print(f"  Processing {cond_name}...")
    
    # Filter cells for this condition - use exact condition column values
    if len(cond_values) == 1:
        condition_mask = adata.obs['condition'] == cond_values[0]
    else:
        # Combine multiple conditions (e.g., Malignant + Matched_malignant)
        condition_mask = adata.obs['condition'].isin(cond_values)
    
    adata_cond = adata[condition_mask].copy()
    
    if adata_cond.n_obs == 0:
        print(f"    Warning: No cells found for {cond_name}")
        ax.text(0.5, 0.5, f'No cells\nfor {cond_name}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(cond_name, fontsize=14, fontweight='bold')
        continue
    
    # Get sample information
    if 'sample_id' in adata_cond.obs.columns:
        unique_samples = adata_cond.obs['sample_id'].unique()
        sample_info = f" ({len(unique_samples)} samples)"
    else:
        sample_info = ""
    
    print(f"    {adata_cond.n_obs:,} cells from condition(s): {', '.join(cond_values)}")
    if 'sample_id' in adata_cond.obs.columns:
        print(f"    Samples: {', '.join(sorted(unique_samples))}")
    
    # Create title with condition name and cell count
    title = f"{cond_name}\n{adata_cond.n_obs:,} cells{sample_info}"
    
    # Plot cells colored by leiden cluster, size/shape by confidence
    for cluster in leiden_clusters:
        cluster_str = str(cluster)
        mask = adata_cond.obs['leiden_1.0'].astype(str) == cluster_str
        if mask.sum() == 0:
            continue
        
        coords = adata_cond.obsm['X_umap'][mask]
        confidences = adata_cond.obs.loc[mask, 'annotation_confidence'].values
        
        # Get leiden cluster color
        leiden_color = leiden_colors[cluster_str]
        
        # Encode confidence by size and shape
        sizes = []
        markers_list = []
        
        for conf in confidences:
            # Size: use wider range and emphasize high confidence cells
            conf_adj = conf ** 0.7  # Emphasize high confidence
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
    
    ax.set_xlabel('UMAP 1', fontsize=11)
    ax.set_ylabel('UMAP 2', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_aspect('equal')

# Add shared shape legend for confidence (on the last plot)
ax_last = axes[-1][0]
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=15, label='High confidence (>70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=10, label='Medium confidence (40-70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=6, label='Low confidence (<40%)', markeredgecolor='black', linewidth=0.5),
]
shape_legend = ax_last.legend(handles=shape_legend_elements, loc='lower right', 
                            fontsize=9, title='Confidence Level', title_fontsize=10,
                            framealpha=0.9, bbox_to_anchor=(0.98, 0.02))

# Add leiden cluster legend on the right side
print("\nCreating leiden cluster legend...")
ax_legend = fig.add_subplot(gs[:, 2])
ax_legend.axis('off')

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

# Use constrained_layout for better handling of subplots
fig.set_constrained_layout(True)

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_annotation_confidence_by_condition.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

print("\nDone!")


