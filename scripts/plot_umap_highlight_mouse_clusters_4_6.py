#!/usr/bin/env python3
"""
Highlight cells from mouse-only leiden clusters 4 and 6 on the human-mouse integration UMAP
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
MOUSE_ONLY_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
HUMAN_MOUSE_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")
FIGDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/figures_individual_samples")
FIGDIR.mkdir(parents=True, exist_ok=True)

# Clusters to highlight
HIGHLIGHT_CLUSTERS = ['4', '6']

print("="*80)
print("Step 1: Loading mouse-only integration")
print("="*80)
adata_mouse = sc.read_h5ad(MOUSE_ONLY_H5AD)
print(f"Mouse-only: {adata_mouse.shape[0]:,} cells × {adata_mouse.shape[1]:,} genes")

# Check leiden clusters
if 'leiden_1.0' not in adata_mouse.obs.columns:
    print("ERROR: leiden_1.0 not found in mouse-only data!")
    exit(1)

print(f"\nLeiden clusters in mouse-only: {sorted(adata_mouse.obs['leiden_1.0'].unique())}")

# Get cell names for clusters 4 and 6
highlight_cells = []
for cluster in HIGHLIGHT_CLUSTERS:
    mask = adata_mouse.obs['leiden_1.0'].astype(str) == cluster
    n_cells = mask.sum()
    print(f"  Cluster {cluster}: {n_cells:,} cells")
    cluster_cells = adata_mouse.obs_names[mask].tolist()
    highlight_cells.extend(cluster_cells)

highlight_cells = set(highlight_cells)
print(f"\nTotal unique cells to highlight: {len(highlight_cells):,}")

print("\n" + "="*80)
print("Step 2: Loading human-mouse integration")
print("="*80)
adata_integrated = sc.read_h5ad(HUMAN_MOUSE_H5AD)
print(f"Human-mouse integration: {adata_integrated.shape[0]:,} cells × {adata_integrated.shape[1]:,} genes")

# Check if UMAP exists
if 'X_umap' not in adata_integrated.obsm:
    print("ERROR: X_umap not found in human-mouse integration!")
    print("Available obsm keys:", list(adata_integrated.obsm.keys()))
    exit(1)

# Find matching cells
print("\n" + "="*80)
print("Step 3: Matching cells between datasets")
print("="*80)

# Check cell name format - might need to match by barcode or full name
mouse_cell_names = set(adata_mouse.obs_names)
integrated_cell_names = set(adata_integrated.obs_names)

# Try direct matching first
matching_cells = highlight_cells & integrated_cell_names
print(f"Direct matches: {len(matching_cells):,} cells")

# If not many matches, try matching by barcode (part before first underscore or dash)
if len(matching_cells) < len(highlight_cells) * 0.5:
    print("\nTrying barcode matching...")
    # Extract barcodes (part before separator)
    mouse_barcodes = {}
    for cell in highlight_cells:
        # Try different separators
        for sep in ['-', '_']:
            if sep in cell:
                barcode = cell.split(sep)[0]
                mouse_barcodes[barcode] = cell
                break
    
    integrated_barcodes = {}
    for cell in adata_integrated.obs_names:
        for sep in ['-', '_']:
            if sep in cell:
                barcode = cell.split(sep)[0]
                integrated_barcodes[barcode] = cell
                break
    
    # Match by barcode
    matching_barcodes = set(mouse_barcodes.keys()) & set(integrated_barcodes.keys())
    matching_cells = {integrated_barcodes[b] for b in matching_barcodes if b in integrated_barcodes}
    print(f"Barcode matches: {len(matching_cells):,} cells")

if len(matching_cells) == 0:
    print("ERROR: No matching cells found!")
    print(f"Sample mouse cell names: {list(highlight_cells)[:5]}")
    print(f"Sample integrated cell names: {list(integrated_cell_names)[:5]}")
    exit(1)

print(f"\nFinal matching cells to highlight: {len(matching_cells):,}")

# Create mask for highlighted cells
highlight_mask = adata_integrated.obs_names.isin(matching_cells)

print("\n" + "="*80)
print("Step 4: Creating separate density plot visualizations")
print("="*80)

# Get all coordinates for extent calculation
coords_all = adata_integrated.obsm['X_umap']
x_min, x_max = coords_all[:, 0].min(), coords_all[:, 0].max()
y_min, y_max = coords_all[:, 1].min(), coords_all[:, 1].max()

# Identify which cells belong to which cluster
cluster_4_cells = set()
cluster_6_cells = set()

for cluster in HIGHLIGHT_CLUSTERS:
    mask = adata_mouse.obs['leiden_1.0'].astype(str) == cluster
    cluster_cells = set(adata_mouse.obs_names[mask])
    
    if cluster == '4':
        cluster_4_cells = cluster_cells & matching_cells
    elif cluster == '6':
        cluster_6_cells = cluster_cells & matching_cells

print(f"Cluster 4 cells: {len(cluster_4_cells):,}")
print(f"Cluster 6 cells: {len(cluster_6_cells):,}")

# Identify tonsil cells for contour overlay
tonsil_mask = (adata_integrated.obs['disease_state'] == 'Tonsil_Normal')
coords_tonsil = adata_integrated.obsm['X_umap'][tonsil_mask]
print(f"Tonsil cells for contour: {tonsil_mask.sum():,}")

# Create separate plots for each cluster
for cluster_num, cluster_cells_set in [('4', cluster_4_cells), ('6', cluster_6_cells)]:
    if len(cluster_cells_set) == 0:
        print(f"Skipping cluster {cluster_num} - no matching cells")
        continue
    
    # Get coordinates for this cluster
    mask_cluster = adata_integrated.obs_names.isin(cluster_cells_set)
    coords_cluster = adata_integrated.obsm['X_umap'][mask_cluster]
    
    from scipy.stats import gaussian_kde
    
    # Create figure with white background
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_facecolor('white')
    
    # Plot all OTHER cells in gray first (background)
    other_mask = ~mask_cluster
    ax.scatter(coords_all[other_mask, 0], coords_all[other_mask, 1],
              c='lightgray', s=5, alpha=0.3,
              marker='o', edgecolors='none', rasterized=True, zorder=1)
    
    # Compute KDE density for the cluster cells
    xy = np.vstack([coords_cluster[:, 0], coords_cluster[:, 1]])
    try:
        kde = gaussian_kde(xy, bw_method=0.15)
        density = kde(xy)
    except Exception:
        density = np.ones(len(coords_cluster))
    
    # Sort by density so densest points are plotted on top
    idx = density.argsort()
    x_sorted = coords_cluster[idx, 0]
    y_sorted = coords_cluster[idx, 1]
    density_sorted = density[idx]
    
    # Plot cluster cells colored by density (viridis)
    sc_plot = ax.scatter(x_sorted, y_sorted,
                        c=density_sorted, cmap='viridis', s=30, alpha=0.9,
                        edgecolors='none', rasterized=True, zorder=3)
    
    # Add colorbar
    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Cell Density (KDE)', rotation=270, labelpad=20, fontsize=12)
    
    # Overlay tonsil cell density contours in black
    if len(coords_tonsil) > 0:
        xy_tonsil = np.vstack([coords_tonsil[:, 0], coords_tonsil[:, 1]])
        kde_tonsil = gaussian_kde(xy_tonsil, bw_method=0.15)
        
        # Create grid for contour
        xx, yy = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        z_tonsil = kde_tonsil(positions).reshape(xx.shape)
        
        # Draw contours
        contour = ax.contour(xx, yy, z_tonsil, levels=6, colors='black', 
                            linewidths=1.5, alpha=0.8, zorder=4)
        ax.clabel(contour, inline=False, fontsize=0)  # no labels on contour lines
        
        # Add legend entry for contour
        from matplotlib.lines import Line2D
        contour_legend = Line2D([0], [0], color='black', linewidth=1.5, 
                               label='Tonsil cell density')
        ax.legend(handles=[contour_legend], loc='lower right', fontsize=10, framealpha=0.9)
    
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'Human-Mouse Integration UMAP:\nDensity of Mouse Leiden Cluster {cluster_num} (viridis) + Tonsil Contours (black)\n({len(cluster_cells_set):,} cells)', 
                 fontsize=14, fontweight='bold')
    
    # Set axis limits to match full UMAP
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    
    # Save
    output_path = FIGDIR / f"umap_all_disease_state_highlight_mouse_cluster_{cluster_num}_density.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {output_path}")
    
    fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    print(f"Saved → {output_path.with_suffix('.pdf')}")
    print(f"Saved → {output_path.with_suffix('.svg')}")
    
    plt.close()

print("\nDone!")


