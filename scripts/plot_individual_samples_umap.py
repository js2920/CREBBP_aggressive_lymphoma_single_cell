#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Individual Samples on UMAP
================================

Loads the integrated AnnData and creates individual UMAP plots
highlighting each sample while showing others in gray.

"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================== CONFIGURATION =================================
# Input file (integrated scVI output)
INTEGRATED_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Output directory for figures
FIGDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/figures_individual_samples")
FIGDIR.mkdir(parents=True, exist_ok=True)

# Group keys to plot by (sample_batch, disease_state, etc.)
GROUP_KEYS = ["sample_batch", "disease_state", "species"]

# Color for highlighted sample
HIGHLIGHT_COLOR = "#d62728"  # Red
BACKGROUND_COLOR = "#d3d3d3"  # Light gray

# Figure settings
FIGSIZE = (10, 9)
POINT_SIZE = 25
DPI = 300

# ============================== LOAD DATA =====================================
print("=" * 80)
print("PLOT INDIVIDUAL SAMPLES ON UMAP")
print("=" * 80)

print(f"\nLoading: {INTEGRATED_H5AD}")
adata = sc.read_h5ad(INTEGRATED_H5AD)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if "X_umap" not in adata.obsm:
    print("  Computing UMAP...")
    if "X_scvi" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=30)
    else:
        sc.pp.neighbors(adata, n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.2, spread=1.5)

# Get UMAP coordinates for consistent axis limits
umap_coords = adata.obsm["X_umap"]
x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
x_margin = (x_max - x_min) * 0.05
y_margin = (y_max - y_min) * 0.05

# ============================== PLOTTING FUNCTIONS ============================
def plot_sample_highlighted(adata, group_key, group_value, outdir, 
                            highlight_color=HIGHLIGHT_COLOR, 
                            background_color=BACKGROUND_COLOR):
    """Plot UMAP with one sample highlighted, others in gray."""
    labels = adata.obs[group_key].astype(str)
    mask = labels == group_value
    n_cells = mask.sum()
    
    fig, ax = plt.subplots(figsize=FIGSIZE)
    
    # Plot background cells first (gray)
    ax.scatter(
        umap_coords[~mask, 0],
        umap_coords[~mask, 1],
        c=background_color,
        s=POINT_SIZE * 0.6,
        alpha=0.3,
        rasterized=True,
        label="Other"
    )
    
    # Plot highlighted cells on top
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=highlight_color,
        s=POINT_SIZE,
        alpha=0.8,
        rasterized=True,
        label=group_value
    )
    
    # Set consistent axis limits
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Styling
    ax.set_title(f"{group_value}\n({n_cells:,} cells)", fontsize=14, fontweight="bold")
    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    # Save
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", group_value)
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_{group_key}_{safe_name}.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    
    return n_cells


def plot_sample_only(adata, group_key, group_value, outdir, color=HIGHLIGHT_COLOR):
    """Plot UMAP showing ONLY cells from one sample (no background)."""
    labels = adata.obs[group_key].astype(str)
    mask = labels == group_value
    n_cells = mask.sum()
    
    fig, ax = plt.subplots(figsize=FIGSIZE)
    
    # Plot only the selected cells
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.7,
        rasterized=True
    )
    
    # Set consistent axis limits (same as full UMAP)
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Styling
    ax.set_title(f"{group_value} only\n({n_cells:,} cells)", fontsize=14, fontweight="bold")
    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    # Save
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", group_value)
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_{group_key}_{safe_name}_only.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    
    return n_cells


def create_grid_plot(adata, group_key, outdir, ncols=4):
    """Create a grid of all samples in one figure."""
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    n_groups = len(unique_groups)
    nrows = (n_groups + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    axes = axes.flatten() if n_groups > 1 else [axes]
    
    for idx, grp in enumerate(unique_groups):
        ax = axes[idx]
        mask = labels == grp
        n_cells = mask.sum()
        
        # Background
        ax.scatter(
            umap_coords[~mask, 0],
            umap_coords[~mask, 1],
            c=BACKGROUND_COLOR,
            s=5,
            alpha=0.2,
            rasterized=True
        )
        
        # Highlighted
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=HIGHLIGHT_COLOR,
            s=8,
            alpha=0.7,
            rasterized=True
        )
        
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_title(f"{grp}\n(n={n_cells:,})", fontsize=9)
        ax.axis("off")
    
    # Hide empty subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].axis("off")
    
    plt.tight_layout()
    
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_grid_{group_key}.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ Grid plot saved: umap_grid_{group_key}")


# ============================== GENERATE PLOTS ================================
for group_key in GROUP_KEYS:
    if group_key not in adata.obs.columns:
        print(f"\n  (skip) '{group_key}' not found in adata.obs")
        continue
    
    print(f"\n{'='*80}")
    print(f"Plotting by: {group_key}")
    print("="*80)
    
    # Create subdirectory for this grouping
    outdir = FIGDIR / group_key
    outdir.mkdir(parents=True, exist_ok=True)
    
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    print(f"  Found {len(unique_groups)} unique groups")
    
    # Plot each group
    for grp in unique_groups:
        n1 = plot_sample_highlighted(adata, group_key, grp, outdir)
        n2 = plot_sample_only(adata, group_key, grp, outdir)
        print(f"    ✓ {grp}: {n1:,} cells")
    
    # Create grid plot
    create_grid_plot(adata, group_key, outdir)

# ============================== SUMMARY PLOT ==================================
print(f"\n{'='*80}")
print("Creating summary plots")
print("="*80)

# Disease state with custom colors
disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000",
    "Tonsil_GC_B": "#4169E1",
    "Tonsil_Normal": "#4169E1"
}

# Species colors
species_colors = {
    "mouse": "#98FB98",
    "human": "#6495ED"
}

# Full UMAP colored by disease state
fig, ax = plt.subplots(figsize=FIGSIZE)
for ds in adata.obs["disease_state"].unique():
    mask = adata.obs["disease_state"] == ds
    color = disease_colors.get(ds, "#808080")
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.6,
        label=f"{ds} ({mask.sum():,})",
        rasterized=True
    )
ax.set_xlim(x_min - x_margin, x_max + x_margin)
ax.set_ylim(y_min - y_margin, y_max + y_margin)
ax.set_title("All Samples by Disease State", fontsize=14, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for fmt in ["png", "pdf", "svg"]:
    fig.savefig(FIGDIR / f"umap_all_disease_state.{fmt}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  ✓ umap_all_disease_state")

# Full UMAP colored by species
fig, ax = plt.subplots(figsize=FIGSIZE)
for sp in adata.obs["species"].unique():
    mask = adata.obs["species"] == sp
    color = species_colors.get(sp, "#808080")
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.6,
        label=f"{sp} ({mask.sum():,})",
        rasterized=True
    )
ax.set_xlim(x_min - x_margin, x_max + x_margin)
ax.set_ylim(y_min - y_margin, y_max + y_margin)
ax.set_title("All Samples by Species", fontsize=14, fontweight="bold")
ax.legend(loc="upper right", fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for fmt in ["png", "pdf", "svg"]:
    fig.savefig(FIGDIR / f"umap_all_species.{fmt}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  ✓ umap_all_species")

# ============================== DONE ==========================================
print(f"\n{'='*80}")
print("COMPLETE!")
print("="*80)
print(f"  Output directory: {FIGDIR}")
print(f"  Total figures generated: {len(list(FIGDIR.rglob('*.png')))}")
print("\nDONE.\n")



