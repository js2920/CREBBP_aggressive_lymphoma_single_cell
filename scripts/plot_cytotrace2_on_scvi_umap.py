#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot CytoTRACE2 Scores on scVI UMAP
====================================

This script loads the integrated data with CytoTRACE2 scores and generates
publication-quality visualizations of differentiation potential.

Author: J
Date: 2025-12-11
"""

import os
import re
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# ============================== CONFIGURATION ================================
# Input: Integrated data with CytoTRACE2 scores
# Try CT2 output first, fall back to main integrated file
CT2_OUTPUT = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/CytoTRACE2/integrated_with_cytotrace2.h5ad")
MAIN_OUTPUT = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Output directory
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/figures_cytotrace2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================== HELPER FUNCTIONS =============================
def save_figure(fig, basename, dpi=300):
    """Save figure in PNG, PDF, and SVG formats."""
    for fmt in ['png', 'pdf', 'svg']:
        filepath = OUTPUT_DIR / f'{basename}.{fmt}'
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
    plt.close(fig)
    print(f"  ✓ {basename} (.png, .pdf, .svg)")


# ============================== MAIN =========================================
print("=" * 84)
print("PLOTTING CYTOTRACE2 SCORES ON scVI UMAP")
print("=" * 84)

# Load data
if CT2_OUTPUT.exists():
    print(f"Loading: {CT2_OUTPUT}")
    adata = sc.read_h5ad(CT2_OUTPUT)
elif MAIN_OUTPUT.exists():
    print(f"Loading: {MAIN_OUTPUT}")
    adata = sc.read_h5ad(MAIN_OUTPUT)
else:
    raise FileNotFoundError("No integrated data found!")

print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Check for CT2 scores
if "cytotrace2_score" not in adata.obs:
    raise ValueError("CytoTRACE2 scores not found in adata.obs! Run CytoTRACE2 first.")

print(f"  CytoTRACE2 scores: found")
print(f"  Score range: {adata.obs['cytotrace2_score'].min():.3f} - {adata.obs['cytotrace2_score'].max():.3f}")

# Color palettes
disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000",
    "Tonsil_Normal": "#4169E1"
}
species_colors = {"mouse": "#98FB98", "human": "#6495ED"}

# Custom colormaps for CT2
ct2_cmap = "viridis"
ct2_cmap_r = "viridis_r"

print(f"\nOutput directory: {OUTPUT_DIR}\n")
print("Generating figures...")

# ==================== 1. Basic CT2 UMAP ====================
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=15, frameon=False,
           cmap=ct2_cmap, title="CytoTRACE2 Score\n(Higher = Less Differentiated)")
save_figure(fig, "umap_cytotrace2_score")

# ==================== 2. CT2 UMAP (reversed) ====================
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=15, frameon=False,
           cmap=ct2_cmap_r, title="CytoTRACE2 Score\n(Darker = Less Differentiated)")
save_figure(fig, "umap_cytotrace2_score_reversed")

# ==================== 3. CT2 split by disease state ====================
if "disease_state" in adata.obs:
    disease_states = sorted(adata.obs["disease_state"].unique())
    n_states = len(disease_states)
    
    fig, axes = plt.subplots(1, n_states, figsize=(5*n_states, 5))
    if n_states == 1:
        axes = [axes]
    
    umap_coords = adata.obsm["X_umap"]
    vmin, vmax = adata.obs["cytotrace2_score"].quantile([0.01, 0.99])
    
    for idx, ds in enumerate(disease_states):
        ax = axes[idx]
        mask = adata.obs["disease_state"] == ds
        
        # Plot background in gray
        ax.scatter(umap_coords[~mask, 0], umap_coords[~mask, 1],
                   c="#e0e0e0", s=5, alpha=0.3, rasterized=True)
        
        # Plot disease state with CT2 colors
        scatter = ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                            c=adata.obs.loc[mask, "cytotrace2_score"],
                            cmap=ct2_cmap, s=15, alpha=0.8, vmin=vmin, vmax=vmax,
                            rasterized=True)
        
        ax.set_title(f"{ds}\n(n={mask.sum():,})", fontsize=12)
        ax.axis("off")
        
        if idx == n_states - 1:
            plt.colorbar(scatter, ax=ax, label="CT2 Score", shrink=0.8)
    
    plt.tight_layout()
    save_figure(fig, "umap_cytotrace2_by_disease_state")

# ==================== 4. CT2 split by species ====================
if "species" in adata.obs:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    umap_coords = adata.obsm["X_umap"]
    vmin, vmax = adata.obs["cytotrace2_score"].quantile([0.01, 0.99])
    
    for idx, sp in enumerate(["human", "mouse"]):
        ax = axes[idx]
        if sp not in adata.obs["species"].values:
            ax.axis("off")
            continue
            
        mask = adata.obs["species"] == sp
        
        # Plot background
        ax.scatter(umap_coords[~mask, 0], umap_coords[~mask, 1],
                   c="#e0e0e0", s=5, alpha=0.3, rasterized=True)
        
        # Plot species with CT2
        scatter = ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                            c=adata.obs.loc[mask, "cytotrace2_score"],
                            cmap=ct2_cmap, s=15, alpha=0.8, vmin=vmin, vmax=vmax,
                            rasterized=True)
        
        ax.set_title(f"{sp.upper()}\n(n={mask.sum():,})", fontsize=14, fontweight="bold")
        ax.axis("off")
        plt.colorbar(scatter, ax=ax, label="CT2 Score", shrink=0.8)
    
    plt.tight_layout()
    save_figure(fig, "umap_cytotrace2_by_species")

# ==================== 5. Violin plot by disease state ====================
if "disease_state" in adata.obs:
    fig, ax = plt.subplots(figsize=(12, 6))
    order = sorted(adata.obs["disease_state"].unique())
    colors = [disease_colors.get(d, "#808080") for d in order]
    
    sns.violinplot(data=adata.obs, x="disease_state", y="cytotrace2_score",
                   order=order, palette=colors, ax=ax, inner="box")
    
    ax.set_xlabel("Disease State", fontsize=12)
    ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
    ax.set_title("Differentiation Potential by Disease State\n(Higher = Less Differentiated)", fontsize=14)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_figure(fig, "violin_cytotrace2_by_disease")

# ==================== 6. Violin plot by species ====================
if "species" in adata.obs:
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(data=adata.obs, x="species", y="cytotrace2_score",
                   palette=species_colors, ax=ax, inner="box")
    ax.set_xlabel("Species", fontsize=12)
    ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
    ax.set_title("Differentiation Potential by Species", fontsize=14)
    save_figure(fig, "violin_cytotrace2_by_species")

# ==================== 7. Boxplot by sample ====================
if "sample_batch" in adata.obs:
    fig, ax = plt.subplots(figsize=(14, 6))
    order = sorted(adata.obs["sample_batch"].unique())
    
    # Color by disease state
    sample_to_disease = adata.obs.groupby("sample_batch")["disease_state"].first().to_dict()
    colors = [disease_colors.get(sample_to_disease.get(s, ""), "#808080") for s in order]
    
    sns.boxplot(data=adata.obs, x="sample_batch", y="cytotrace2_score",
                order=order, palette=colors, ax=ax)
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
    ax.set_title("Differentiation Potential by Sample", fontsize=14)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_figure(fig, "boxplot_cytotrace2_by_sample")

# ==================== 8. Histogram ====================
fig, ax = plt.subplots(figsize=(10, 6))
scores = adata.obs["cytotrace2_score"].dropna()
ax.hist(scores, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
ax.axvline(scores.median(), color="red", linestyle="--", linewidth=2, label=f"Median: {scores.median():.3f}")
ax.axvline(scores.mean(), color="green", linestyle="--", linewidth=2, label=f"Mean: {scores.mean():.3f}")
ax.set_xlabel("CytoTRACE2 Score", fontsize=12)
ax.set_ylabel("Cell Count", fontsize=12)
ax.set_title("CytoTRACE2 Score Distribution", fontsize=14)
ax.legend(fontsize=11)
save_figure(fig, "histogram_cytotrace2_score")

# ==================== 9. Histogram split by disease state ====================
if "disease_state" in adata.obs:
    fig, ax = plt.subplots(figsize=(12, 6))
    for ds in sorted(adata.obs["disease_state"].unique()):
        mask = adata.obs["disease_state"] == ds
        scores = adata.obs.loc[mask, "cytotrace2_score"].dropna()
        ax.hist(scores, bins=40, alpha=0.5, label=f"{ds} (n={len(scores):,})",
                color=disease_colors.get(ds, "#808080"), density=True)
    
    ax.set_xlabel("CytoTRACE2 Score", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("CytoTRACE2 Score Distribution by Disease State", fontsize=14)
    ax.legend(fontsize=10)
    save_figure(fig, "histogram_cytotrace2_by_disease")

# ==================== 10. Potency categories ====================
if "cytotrace2_potency" in adata.obs:
    fig, ax = plt.subplots(figsize=(10, 6))
    potency_counts = adata.obs["cytotrace2_potency"].value_counts()
    potency_counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="black")
    ax.set_xlabel("Potency Category", fontsize=12)
    ax.set_ylabel("Cell Count", fontsize=12)
    ax.set_title("CytoTRACE2 Potency Distribution", fontsize=14)
    plt.xticks(rotation=45, ha="right")
    save_figure(fig, "barplot_potency_categories")
    
    # Potency by disease state (stacked bar)
    if "disease_state" in adata.obs:
        crosstab = pd.crosstab(adata.obs["disease_state"], adata.obs["cytotrace2_potency"])
        crosstab_pct = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
        
        fig, ax = plt.subplots(figsize=(12, 7))
        crosstab_pct.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
        ax.set_xlabel("Disease State", fontsize=12)
        ax.set_ylabel("Percentage", fontsize=12)
        ax.set_title("Potency Distribution by Disease State", fontsize=14)
        ax.legend(title="Potency", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        save_figure(fig, "stacked_bar_potency_by_disease")

# ==================== 11. CT2 vs other scores ====================
for score_col in ["S_score", "G2M_score", "oxphos_score", "bcr_score"]:
    if score_col in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Sample for speed
        n_sample = min(10000, adata.n_obs)
        idx = np.random.choice(adata.n_obs, n_sample, replace=False)
        
        x = adata.obs["cytotrace2_score"].iloc[idx]
        y = adata.obs[score_col].iloc[idx]
        
        ax.scatter(x, y, c="#404040", s=5, alpha=0.3, rasterized=True)
        
        # Add correlation
        valid = ~(x.isna() | y.isna())
        if valid.sum() > 10:
            corr = np.corrcoef(x[valid], y[valid])[0, 1]
            ax.text(0.05, 0.95, f"r = {corr:.3f}", transform=ax.transAxes,
                    fontsize=12, verticalalignment='top')
        
        ax.set_xlabel("CytoTRACE2 Score", fontsize=12)
        ax.set_ylabel(score_col.replace("_", " ").title(), fontsize=12)
        ax.set_title(f"CT2 vs {score_col.replace('_', ' ').title()}", fontsize=14)
        save_figure(fig, f"scatter_ct2_vs_{score_col}")

# ==================== Statistics ====================
print("\n" + "=" * 84)
print("STATISTICS")
print("=" * 84)

scores = adata.obs["cytotrace2_score"].dropna()
print(f"\nOverall CT2 statistics:")
print(f"  Mean: {scores.mean():.4f}")
print(f"  Median: {scores.median():.4f}")
print(f"  Std: {scores.std():.4f}")
print(f"  Min: {scores.min():.4f}")
print(f"  Max: {scores.max():.4f}")

if "disease_state" in adata.obs:
    print("\nCT2 by disease state:")
    stats = adata.obs.groupby("disease_state")["cytotrace2_score"].agg(
        ["mean", "median", "std", "count"]
    ).round(4)
    print(stats.to_string())
    stats.to_csv(OUTPUT_DIR / "statistics_ct2_by_disease.csv")

if "species" in adata.obs:
    print("\nCT2 by species:")
    stats = adata.obs.groupby("species")["cytotrace2_score"].agg(
        ["mean", "median", "std", "count"]
    ).round(4)
    print(stats.to_string())
    stats.to_csv(OUTPUT_DIR / "statistics_ct2_by_species.csv")

# ==================== DONE ====================
print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTPUT_DIR}")
print(f"  Total figures generated: multiple")
print("\nDONE.\n")



