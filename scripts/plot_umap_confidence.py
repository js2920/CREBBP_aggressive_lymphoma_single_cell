#!/usr/bin/env python3
"""
Generate UMAP with color intensity proportional to annotation confidence
Based on majority_celltype percentage within each leiden cluster
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from colorsys import rgb_to_hsv, hsv_to_rgb

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
conf_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_pct'] / 100.0))  # Convert to 0-1

# Map to cells
print("Mapping majority celltypes and confidence to cells...")
adata.obs['majority_celltype'] = adata.obs['leiden_1.0'].astype(str).map(comp_dict)
adata.obs['annotation_confidence'] = adata.obs['leiden_1.0'].astype(str).map(conf_dict)

# Fill NaN values (if any clusters not in composition file)
adata.obs['majority_celltype'] = adata.obs['majority_celltype'].fillna('Unknown')
adata.obs['annotation_confidence'] = adata.obs['annotation_confidence'].fillna(0.0)

conf_min = adata.obs['annotation_confidence'].min()
conf_max = adata.obs['annotation_confidence'].max()
print(f"Confidence range: {conf_min:.3f} - {conf_max:.3f}")
print(f"Size range: {1 + 14 * conf_min:.1f} - {1 + 14 * conf_max:.1f} pixels (directly proportional)")
print(f"Number of celltypes: {adata.obs['majority_celltype'].nunique()}")

# Get unique celltypes and create a color palette
celltypes = sorted(adata.obs['majority_celltype'].unique())
n_types = len(celltypes)

# Use a qualitative colormap (tab20 or similar, cycled if needed)
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')  # Fallback for older matplotlib

colors = {}
for i, ct in enumerate(celltypes):
    colors[ct] = base_cmap(i % 20)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Plot each celltype with varying intensity based on confidence
# Use color intensity (saturation) where darker = higher confidence
for celltype in celltypes:
    mask = adata.obs['majority_celltype'] == celltype
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Use full base color for all cells (no intensity variation)
    base_color = colors[celltype][:3]  # RGB tuple
    
    # Encode confidence only by size and shape, not color
    sizes = []
    markers_list = []
    
    # Define markers for different confidence levels
    # High conf (>0.7): circle 'o'
    # Medium conf (0.4-0.7): square 's'  
    # Low conf (<0.4): triangle '^'
    
    for conf in confidences:
        # Size: directly proportional to confidence score
        # Use wider range for better visibility: 1 to 15
        # conf=0.258 (min) -> size ~4, conf=0.931 (max) -> size ~14
        size = 1 + 14 * conf  # Directly proportional: size = 1 + 14 * confidence
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
    
    # Plot by marker type for better control
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_list) == marker_type
        if mask_marker.sum() == 0:
            continue
        # Use same full color for all points of this celltype
        ax.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                  color=base_color, s=sizes[mask_marker], 
                  alpha=0.7, marker=marker_type,
                  label=celltype if len(celltypes) <= 20 and marker_type == 'o' else None,
                  edgecolors='black', linewidths=0.3, rasterized=True)

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title('UMAP colored by majority celltype\n(Larger + Circle = High Conf | Smaller + Triangle = Low Conf)', 
             fontsize=13, fontweight='bold')

# Add legend for cell types (limit to avoid overcrowding)
if len(celltypes) <= 20:
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, 
              title='Cell Type', title_fontsize=9, framealpha=0.9)

# Add shape legend for confidence
from matplotlib.lines import Line2D
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=8, label='High confidence (>70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=6, label='Medium confidence (40-70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=5, label='Low confidence (<40%)', markeredgecolor='black'),
]
shape_legend = ax.legend(handles=shape_legend_elements, loc='lower right', 
                         fontsize=9, title='Confidence Level', title_fontsize=10,
                         framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax.add_artist(shape_legend)  # Keep both legends

# Add colorbar for confidence (reference scale)
# Note: confidence is encoded by size and shape, not color
sm = plt.cm.ScalarMappable(cmap=plt.cm.Greys, 
                           norm=plt.Normalize(vmin=0, vmax=100))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Annotation Confidence\n(% majority in cluster)\n(Encoded by size & shape)', 
               rotation=270, labelpad=25, fontsize=10)

plt.tight_layout()

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_celltype_confidence.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nSaved → {output_path}")

# Also save PDF and SVG
fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

# ===== ALTERNATIVE: HSV-based visualization (even more dramatic) =====
print("\nCreating alternative HSV-based visualization...")
fig2, ax2 = plt.subplots(figsize=(12, 10))

for celltype in celltypes:
    mask = adata.obs['majority_celltype'] == celltype
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Use full base color for all cells (no intensity variation)
    base_color = colors[celltype][:3]  # RGB tuple
    
    # Encode confidence only by size and shape, not color
    sizes_hsv = []
    markers_hsv = []
    
    for conf in confidences:
        # Size: directly proportional to confidence score
        # Use wider range for better visibility: 1 to 15
        # conf=0.258 (min) -> size ~4, conf=0.931 (max) -> size ~14
        size = 1 + 14 * conf  # Directly proportional: size = 1 + 14 * confidence
        sizes_hsv.append(size)
        
        # Shape based on confidence
        if conf > 0.7:
            marker = 'o'  # Circle for high confidence
        elif conf > 0.4:
            marker = 's'  # Square for medium confidence
        else:
            marker = '^'  # Triangle for low confidence
        markers_hsv.append(marker)
    
    sizes_hsv = np.array(sizes_hsv)
    
    # Plot by marker type
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_hsv) == marker_type
        if mask_marker.sum() == 0:
            continue
        ax2.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                   color=base_color, s=sizes_hsv[mask_marker], 
                   alpha=0.7, marker=marker_type,
                   label=celltype if len(celltypes) <= 20 and marker_type == 'o' else None,
                   edgecolors='black', linewidths=0.3, rasterized=True)

ax2.set_xlabel('UMAP 1', fontsize=12)
ax2.set_ylabel('UMAP 2', fontsize=12)
ax2.set_title('UMAP colored by majority celltype\n(Larger + Circle = High Conf | Smaller + Triangle = Low Conf)', 
             fontsize=13, fontweight='bold')

if len(celltypes) <= 20:
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, 
              title='Cell Type', title_fontsize=9, framealpha=0.9)

# Add shape legend for HSV version too
shape_legend_elements2 = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=8, label='High confidence (>70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=6, label='Medium confidence (40-70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=5, label='Low confidence (<40%)', markeredgecolor='black'),
]
shape_legend2 = ax2.legend(handles=shape_legend_elements2, loc='lower right', 
                           fontsize=9, title='Confidence Level', title_fontsize=10,
                           framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax2.add_artist(shape_legend2)

# Add colorbar for confidence (reference scale)
sm2 = plt.cm.ScalarMappable(cmap=plt.cm.Greys, 
                            norm=plt.Normalize(vmin=0, vmax=100))
sm2.set_array([])
cbar2 = plt.colorbar(sm2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.set_label('Annotation Confidence\n(% majority in cluster)\n(Encoded by size & shape)', 
               rotation=270, labelpad=25, fontsize=10)

plt.tight_layout()

# Save alternative version
output_path_hsv = FIGDIR / "umap_leiden_1.0_majority_celltype_confidence_hsv.png"
fig2.savefig(output_path_hsv, dpi=300, bbox_inches='tight')
print(f"Saved → {output_path_hsv}")
fig2.savefig(output_path_hsv.with_suffix('.pdf'), bbox_inches='tight')
fig2.savefig(output_path_hsv.with_suffix('.svg'), bbox_inches='tight')

plt.close()

print("\nDone!")


