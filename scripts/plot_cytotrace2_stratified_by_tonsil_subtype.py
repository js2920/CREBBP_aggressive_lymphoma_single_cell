#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stratified CytoTRACE2 Analysis by Tonsil Cell Subtype
=====================================================

This script analyzes CytoTRACE2 scores stratified by:
- Tonsil subtypes (DZ proliferative, DZ non-proliferative, LZ, Memory B, etc.)
- DLBCL
- Mouse malignant

Helps determine if specific tonsil populations drive the high CT2 scores.

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
from scipy import stats

# ============================== CONFIGURATION ================================
# Input files
CT2_OUTPUT = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/CytoTRACE2/integrated_with_cytotrace2.h5ad")
MAIN_OUTPUT = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Tonsil source files (to get original annotations)
TONSIL_GC_PATH = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_GCBC_RNA.h5ad")
TONSIL_MBC_PATH = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad")

# Output
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/figures_cytotrace2_stratified")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================== HELPER FUNCTIONS =============================
def save_figure(fig, basename, dpi=300):
    """Save figure in PNG, PDF, and SVG formats."""
    for fmt in ['png', 'pdf', 'svg']:
        filepath = OUTPUT_DIR / f'{basename}.{fmt}'
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
    plt.close(fig)
    print(f"  ✓ {basename} (.png, .pdf, .svg)")


def find_label_column(adata):
    """Find cell type annotation column."""
    candidates = [
        "annotation_20230508", "annotation_20220414", "cell_type", "CellType",
        "celltype", "label", "celltype_l2", "celltype_l1"
    ]
    for col in candidates:
        if col in adata.obs.columns:
            return col
    for col in adata.obs.columns:
        if re.search(r"(cell.?type|annotation|label)", col, re.I):
            return col
    return None


# ============================== MAIN =========================================
print("=" * 84)
print("STRATIFIED CYTOTRACE2 ANALYSIS BY TONSIL CELL SUBTYPE")
print("=" * 84)

# Load integrated data
if CT2_OUTPUT.exists():
    print(f"Loading: {CT2_OUTPUT}")
    adata = sc.read_h5ad(CT2_OUTPUT)
elif MAIN_OUTPUT.exists():
    print(f"Loading: {MAIN_OUTPUT}")
    adata = sc.read_h5ad(MAIN_OUTPUT)
else:
    raise FileNotFoundError("No integrated data found!")

print(f"  Loaded: {adata.n_obs:,} cells")

if "cytotrace2_score" not in adata.obs:
    raise ValueError("CytoTRACE2 scores not found!")

# ==================== Load tonsil annotations ====================
print("\nLoading tonsil annotations from source files...")

# Create mapping from cell barcode to original annotation
tonsil_annotations = {}

# Load GC annotations
if TONSIL_GC_PATH.exists():
    print(f"  Loading GC annotations from: {TONSIL_GC_PATH.name}")
    tonsil_gc = sc.read_h5ad(TONSIL_GC_PATH, backed="r")
    label_col = find_label_column(tonsil_gc)
    if label_col:
        for bc, label in zip(tonsil_gc.obs_names, tonsil_gc.obs[label_col]):
            tonsil_annotations[str(bc)] = str(label)
        print(f"    Found {len(tonsil_annotations):,} annotations")
    if hasattr(tonsil_gc, 'file') and tonsil_gc.file is not None:
        try: tonsil_gc.file.close()
        except: pass
    del tonsil_gc

# Load MBC annotations
if TONSIL_MBC_PATH.exists():
    print(f"  Loading MBC annotations from: {TONSIL_MBC_PATH.name}")
    tonsil_mbc = sc.read_h5ad(TONSIL_MBC_PATH, backed="r")
    label_col = find_label_column(tonsil_mbc)
    if label_col:
        n_before = len(tonsil_annotations)
        for bc, label in zip(tonsil_mbc.obs_names, tonsil_mbc.obs[label_col]):
            tonsil_annotations[str(bc)] = str(label)
        print(f"    Added {len(tonsil_annotations) - n_before:,} MBC annotations")
    if hasattr(tonsil_mbc, 'file') and tonsil_mbc.file is not None:
        try: tonsil_mbc.file.close()
        except: pass
    del tonsil_mbc

print(f"  Total tonsil annotations: {len(tonsil_annotations):,}")

# ==================== Add annotations to integrated data ====================
print("\nMapping annotations to integrated data...")

# Map tonsil cell types
adata.obs["tonsil_subtype"] = adata.obs_names.map(
    lambda x: tonsil_annotations.get(str(x), None)
)

# Create unified cell type column
def get_cell_category(row):
    if row["disease_state"] == "Tonsil_Normal" and pd.notna(row.get("tonsil_subtype")):
        return row["tonsil_subtype"]
    elif row["disease_state"] == "DLBCL":
        return "DLBCL"
    elif "Mouse" in str(row["disease_state"]):
        return row["disease_state"]
    else:
        return row["disease_state"]

adata.obs["cell_category"] = adata.obs.apply(get_cell_category, axis=1)

# Print distribution
print("\nCell category distribution:")
cat_counts = adata.obs["cell_category"].value_counts()
for cat, count in cat_counts.items():
    pct = 100 * count / adata.n_obs
    print(f"  {cat}: {count:,} ({pct:.1f}%)")

# ==================== Color palettes ====================
# Custom colors for tonsil subtypes
tonsil_colors = {
    # Dark Zone (proliferating)
    "DZ late Sphase": "#e74c3c",
    "DZ early Sphase": "#c0392b",
    "DZ late G2Mphase": "#e67e22",
    "DZ early G2Mphase": "#d35400",
    # Dark Zone (non-proliferating)
    "DZ non proliferative": "#9b59b6",
    "DZ cell cycle exit": "#8e44ad",
    "GC DZ Noproli": "#7d3c98",
    # Light Zone
    "LZ": "#3498db",
    "LZ proliferative": "#2980b9",
    "PC committed Light Zone GCBC": "#1abc9c",
    # Memory B cells
    "MBC FCRL5+": "#27ae60",
    "Early MBC": "#2ecc71",
    # Malignant
    "DLBCL": "#8B0000",
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
}

# Get categories present in data
categories_present = [c for c in adata.obs["cell_category"].unique() if pd.notna(c)]
palette = {c: tonsil_colors.get(c, "#808080") for c in categories_present}

# ==================== FIGURE 1: Violin plot by cell category ====================
print("\nGenerating figures...")

# Sort categories by median CT2 score
cat_medians = adata.obs.groupby("cell_category")["cytotrace2_score"].median().sort_values(ascending=False)
order = list(cat_medians.index)

fig, ax = plt.subplots(figsize=(16, 8))
sns.violinplot(data=adata.obs, x="cell_category", y="cytotrace2_score",
               order=order, palette=palette, ax=ax, inner="box", cut=0)

ax.set_xlabel("Cell Type", fontsize=12)
ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
ax.set_title("CytoTRACE2 Score by Cell Type\n(Sorted by Median, Higher = Less Differentiated)", fontsize=14)
plt.xticks(rotation=60, ha="right", fontsize=9)

# Add horizontal line at DLBCL median for reference
if "DLBCL" in cat_medians.index:
    dlbcl_median = cat_medians["DLBCL"]
    ax.axhline(dlbcl_median, color="#8B0000", linestyle="--", alpha=0.7, label=f"DLBCL median: {dlbcl_median:.3f}")
    ax.legend(loc="upper right")

plt.tight_layout()
save_figure(fig, "violin_ct2_by_cell_category_all")

# ==================== FIGURE 2: Boxplot with statistics ====================
fig, ax = plt.subplots(figsize=(16, 8))
sns.boxplot(data=adata.obs, x="cell_category", y="cytotrace2_score",
            order=order, palette=palette, ax=ax)

ax.set_xlabel("Cell Type", fontsize=12)
ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
ax.set_title("CytoTRACE2 Score by Cell Type (Boxplot)", fontsize=14)
plt.xticks(rotation=60, ha="right", fontsize=9)

# Add sample sizes
for i, cat in enumerate(order):
    n = (adata.obs["cell_category"] == cat).sum()
    ax.text(i, ax.get_ylim()[0] - 0.02, f"n={n:,}", ha="center", fontsize=7, rotation=90)

plt.tight_layout()
save_figure(fig, "boxplot_ct2_by_cell_category_all")

# ==================== FIGURE 3: Grouped by category type ====================
# Group into: DZ proliferative, DZ non-proliferative, LZ, Memory, Malignant
def get_broad_category(cat):
    cat_lower = str(cat).lower()
    if "dlbcl" in cat_lower:
        return "DLBCL"
    elif "mouse" in cat_lower:
        return "Mouse Malignant"
    elif any(x in cat_lower for x in ["sphase", "g2mphase", "proliferative"]):
        if "dz" in cat_lower or "dark" in cat_lower:
            return "DZ (Proliferating)"
        elif "lz" in cat_lower or "light" in cat_lower:
            return "LZ (Proliferating)"
        else:
            return "Proliferating"
    elif "dz" in cat_lower or "dark" in cat_lower:
        return "DZ (Non-proliferating)"
    elif "lz" in cat_lower or "light" in cat_lower or "pc committed" in cat_lower:
        return "LZ"
    elif "mbc" in cat_lower or "memory" in cat_lower:
        return "Memory B"
    else:
        return "Other"

adata.obs["broad_category"] = adata.obs["cell_category"].apply(get_broad_category)

broad_order = ["DZ (Proliferating)", "DZ (Non-proliferating)", "LZ (Proliferating)", "LZ", 
               "Memory B", "DLBCL", "Mouse Malignant", "Other"]
broad_order = [c for c in broad_order if c in adata.obs["broad_category"].values]

broad_colors = {
    "DZ (Proliferating)": "#e74c3c",
    "DZ (Non-proliferating)": "#9b59b6",
    "LZ (Proliferating)": "#2980b9",
    "LZ": "#3498db",
    "Memory B": "#27ae60",
    "DLBCL": "#8B0000",
    "Mouse Malignant": "#FF6B6B",
    "Other": "#808080"
}

fig, ax = plt.subplots(figsize=(12, 7))
sns.violinplot(data=adata.obs, x="broad_category", y="cytotrace2_score",
               order=broad_order, palette=broad_colors, ax=ax, inner="box")

ax.set_xlabel("Cell Category", fontsize=12)
ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
ax.set_title("CytoTRACE2 Score by Broad Category\n(Higher = Less Differentiated)", fontsize=14)
plt.xticks(rotation=45, ha="right")

# Add sample sizes
for i, cat in enumerate(broad_order):
    n = (adata.obs["broad_category"] == cat).sum()
    med = adata.obs.loc[adata.obs["broad_category"] == cat, "cytotrace2_score"].median()
    ax.text(i, ax.get_ylim()[1] + 0.01, f"n={n:,}\nmed={med:.2f}", ha="center", fontsize=8)

plt.tight_layout()
save_figure(fig, "violin_ct2_by_broad_category")

# ==================== FIGURE 4: UMAP colored by cell category ====================
if "X_umap" in adata.obsm:
    fig, ax = plt.subplots(figsize=(14, 10))
    sc.pl.umap(adata, color="cell_category", ax=ax, show=False, frameon=False,
               legend_loc="right margin", legend_fontsize=7, s=10, palette=palette,
               title="Cell Categories on scVI UMAP")
    save_figure(fig, "umap_cell_categories")

# ==================== FIGURE 5: UMAP colored by broad category ====================
if "X_umap" in adata.obsm:
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color="broad_category", ax=ax, show=False, frameon=False,
               legend_loc="right margin", legend_fontsize=10, s=10, palette=broad_colors,
               title="Broad Categories on scVI UMAP")
    save_figure(fig, "umap_broad_categories")

# ==================== FIGURE 6: Histogram overlays ====================
fig, ax = plt.subplots(figsize=(12, 6))
for cat in broad_order:
    mask = adata.obs["broad_category"] == cat
    scores = adata.obs.loc[mask, "cytotrace2_score"].dropna()
    if len(scores) > 0:
        ax.hist(scores, bins=40, alpha=0.4, label=f"{cat} (n={len(scores):,})",
                color=broad_colors.get(cat, "#808080"), density=True)

ax.set_xlabel("CytoTRACE2 Score", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("CT2 Score Distribution by Broad Category", fontsize=14)
ax.legend(fontsize=9, loc="upper left")
save_figure(fig, "histogram_ct2_by_broad_category")

# ==================== FIGURE 7: Tonsil only - by subtype ====================
tonsil_mask = adata.obs["disease_state"] == "Tonsil_Normal"
if tonsil_mask.sum() > 0:
    tonsil_data = adata.obs[tonsil_mask].copy()
    
    # Get tonsil subtypes
    tonsil_subtypes = tonsil_data["tonsil_subtype"].dropna().unique()
    
    if len(tonsil_subtypes) > 0:
        # Sort by median
        subtype_medians = tonsil_data.groupby("tonsil_subtype")["cytotrace2_score"].median().sort_values(ascending=False)
        tonsil_order = list(subtype_medians.index)
        
        fig, ax = plt.subplots(figsize=(14, 7))
        sns.violinplot(data=tonsil_data, x="tonsil_subtype", y="cytotrace2_score",
                       order=tonsil_order, palette=tonsil_colors, ax=ax, inner="box")
        
        ax.set_xlabel("Tonsil Subtype", fontsize=12)
        ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
        ax.set_title("CytoTRACE2 Score by Tonsil Subtype\n(Sorted by Median)", fontsize=14)
        plt.xticks(rotation=60, ha="right", fontsize=9)
        
        # Add DLBCL median reference line
        if "DLBCL" in cat_medians.index:
            dlbcl_median = cat_medians["DLBCL"]
            ax.axhline(dlbcl_median, color="#8B0000", linestyle="--", alpha=0.7, 
                       label=f"DLBCL median: {dlbcl_median:.3f}")
            ax.legend(loc="upper right")
        
        plt.tight_layout()
        save_figure(fig, "violin_ct2_tonsil_subtypes_only")

# ==================== FIGURE 8: Compare specific populations ====================
# Compare LZ (non-proliferating) vs DLBCL vs Mouse
compare_cats = ["LZ", "DLBCL", "Mouse_Malignant", "Mouse_Matched_malignant"]
compare_cats = [c for c in compare_cats if c in adata.obs["cell_category"].values]

if len(compare_cats) > 1:
    compare_mask = adata.obs["cell_category"].isin(compare_cats)
    compare_data = adata.obs[compare_mask].copy()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=compare_data, x="cell_category", y="cytotrace2_score",
                   order=compare_cats, palette=palette, ax=ax, inner="box")
    
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
    ax.set_title("CT2: LZ (Non-proliferating) vs Malignant\n(Excluding Proliferating Cells)", fontsize=14)
    
    # Add statistics
    for i, cat in enumerate(compare_cats):
        n = (compare_data["cell_category"] == cat).sum()
        med = compare_data.loc[compare_data["cell_category"] == cat, "cytotrace2_score"].median()
        ax.text(i, ax.get_ylim()[1] + 0.01, f"n={n:,}\nmed={med:.2f}", ha="center", fontsize=9)
    
    plt.tight_layout()
    save_figure(fig, "violin_ct2_lz_vs_malignant")

# ==================== STATISTICS ====================
print("\n" + "=" * 84)
print("STATISTICS")
print("=" * 84)

# By cell category
print("\nCT2 by Cell Category:")
stats_df = adata.obs.groupby("cell_category")["cytotrace2_score"].agg(
    ["count", "mean", "median", "std", "min", "max"]
).round(4).sort_values("median", ascending=False)
print(stats_df.to_string())
stats_df.to_csv(OUTPUT_DIR / "statistics_ct2_by_cell_category.csv")

# By broad category
print("\nCT2 by Broad Category:")
broad_stats = adata.obs.groupby("broad_category")["cytotrace2_score"].agg(
    ["count", "mean", "median", "std"]
).round(4).sort_values("median", ascending=False)
print(broad_stats.to_string())
broad_stats.to_csv(OUTPUT_DIR / "statistics_ct2_by_broad_category.csv")

# Statistical tests
print("\n" + "-" * 40)
print("Statistical comparisons (Mann-Whitney U):")
print("-" * 40)

# Compare DLBCL vs each tonsil subtype
if "DLBCL" in adata.obs["cell_category"].values:
    dlbcl_scores = adata.obs.loc[adata.obs["cell_category"] == "DLBCL", "cytotrace2_score"].dropna()
    
    comparisons = []
    for cat in adata.obs["cell_category"].unique():
        if cat == "DLBCL" or "Mouse" in str(cat):
            continue
        cat_scores = adata.obs.loc[adata.obs["cell_category"] == cat, "cytotrace2_score"].dropna()
        if len(cat_scores) > 10:
            stat, pval = stats.mannwhitneyu(dlbcl_scores, cat_scores, alternative='two-sided')
            comparisons.append({
                "comparison": f"DLBCL vs {cat}",
                "n_dlbcl": len(dlbcl_scores),
                "n_other": len(cat_scores),
                "median_dlbcl": dlbcl_scores.median(),
                "median_other": cat_scores.median(),
                "U_statistic": stat,
                "p_value": pval,
                "significant": pval < 0.05
            })
    
    if comparisons:
        comp_df = pd.DataFrame(comparisons)
        comp_df = comp_df.sort_values("p_value")
        print(comp_df.to_string(index=False))
        comp_df.to_csv(OUTPUT_DIR / "statistical_tests_dlbcl_vs_tonsil.csv", index=False)

# ==================== DONE ====================
print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTPUT_DIR}")
print(f"\n  Key findings to check:")
print(f"    - Which tonsil subtypes have highest CT2?")
print(f"    - Is LZ (non-proliferating) still higher than DLBCL?")
print(f"    - Are DZ proliferating cells driving the high tonsil scores?")
print("\nDONE.\n")




