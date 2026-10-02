#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Downstream Plotting Script for Human DLBCL + Mouse Malignant Integration
=========================================================================

Generates:
1. Violin plots for CytoTRACE2 score distributions by disease_state
2. Violin plots for BCR score distributions by disease_state
3. Violin plots for OXPHOS score distributions by disease_state
4. Majority cell type (from Geneformer) per Leiden cluster (resolution 1.0)
5. Comparisons by species (mouse vs human)

Outputs saved as PNG, SVG, and PDF.

Author: J
Date: 2025-12-02
"""

import os
os.environ["OMP_NUM_THREADS"] = "4"
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# ============================== PATHS ========================================
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration")
FIGDIR = OUTDIR / "Geneformer" / "figures"
H5AD_PATH = OUTDIR / "integrated_human_dlbcl_mouse_malignant.h5ad"
GENEFORMER_H5AD = OUTDIR / "Geneformer" / "integrated_with_geneformer_predictions.h5ad"

FIGDIR.mkdir(parents=True, exist_ok=True)

# ============================== LOAD DATA ====================================
print("=" * 84)
print("Loading integrated AnnData...")
print("=" * 84)

if not H5AD_PATH.exists():
    raise FileNotFoundError(f"AnnData file not found: {H5AD_PATH}")

adata = sc.read_h5ad(H5AD_PATH)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  Disease states: {adata.obs['disease_state'].unique().tolist()}")
print(f"  Species: {adata.obs['species'].unique().tolist()}")

# Load Geneformer predictions
print("\n  Loading Geneformer predictions...")
if GENEFORMER_H5AD.exists():
    adata_gf = sc.read_h5ad(GENEFORMER_H5AD)
    print(f"  Geneformer file: {adata_gf.n_obs:,} cells")
    
    # Transfer geneformer_predicted_celltype to main adata
    if 'geneformer_predicted_celltype' in adata_gf.obs.columns:
        # Match by cell barcode index
        common_cells = adata.obs_names.intersection(adata_gf.obs_names)
        print(f"  Matching cells: {len(common_cells):,}")
        
        adata.obs['geneformer_predicted_celltype'] = pd.NA
        adata.obs.loc[common_cells, 'geneformer_predicted_celltype'] = \
            adata_gf.obs.loc[common_cells, 'geneformer_predicted_celltype'].values
        
        if 'geneformer_confidence' in adata_gf.obs.columns:
            adata.obs['geneformer_confidence'] = pd.NA
            adata.obs.loc[common_cells, 'geneformer_confidence'] = \
                adata_gf.obs.loc[common_cells, 'geneformer_confidence'].values
        
        print(f"  ✓ Transferred geneformer_predicted_celltype")
        n_celltypes = adata.obs['geneformer_predicted_celltype'].dropna().nunique()
        print(f"  Unique cell types: {n_celltypes}")
    else:
        print(f"  (warn) geneformer_predicted_celltype not found in Geneformer file")
else:
    print(f"  (warn) Geneformer file not found: {GENEFORMER_H5AD}")
    print(f"  Run geneformer_predict_human_mouse_integration.py first!")

# ============================== COLOR PALETTES ===============================
disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000"
}

species_colors = {
    "mouse": "#98FB98",
    "human": "#6495ED"
}

# Define order for consistent plotting
disease_order = ["Mouse_Matched_malignant", "Mouse_Malignant", "DLBCL"]
disease_order = [d for d in disease_order if d in adata.obs['disease_state'].unique()]

species_order = ["mouse", "human"]
species_order = [s for s in species_order if s in adata.obs['species'].unique()]

# ============================== HELPER FUNCTION ==============================
def save_figure(fig, figdir, filename_base):
    """Save figure in PNG, SVG, and PDF formats."""
    for ext in ['png', 'svg', 'pdf']:
        filepath = figdir / f"{filename_base}.{ext}"
        fig.savefig(filepath, dpi=300, bbox_inches="tight", format=ext)
    print(f"  ✓ {filename_base}.{{png,svg,pdf}}")


def save_violin_plot(adata, score_key, groupby, order, palette, title, filename_base, figdir):
    """Save a violin plot for a given score by group."""
    if score_key not in adata.obs.columns:
        print(f"  (skip) '{score_key}' not found in adata.obs")
        return
    
    # Prepare data
    df = adata.obs[[score_key, groupby]].copy()
    df = df.dropna(subset=[score_key])
    df[score_key] = df[score_key].astype(float)
    
    # Filter to order categories present
    present_order = [c for c in order if c in df[groupby].unique()]
    df = df[df[groupby].isin(present_order)]
    
    if len(df) == 0:
        print(f"  (skip) No data for '{score_key}'")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create violin plot with seaborn
    sns.violinplot(
        data=df,
        x=groupby,
        y=score_key,
        order=present_order,
        palette=[palette.get(c, "#999999") for c in present_order],
        inner="box",
        ax=ax,
        cut=0,
        scale="width"
    )
    
    # Styling
    ax.set_xlabel(groupby.replace("_", " ").title(), fontsize=14, fontweight='bold')
    ax.set_ylabel(score_key.replace("_", " ").title(), fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=12)
    ax.tick_params(axis='y', labelsize=11)
    
    # Add sample sizes
    for i, cat in enumerate(present_order):
        n = (df[groupby] == cat).sum()
        ax.text(i, ax.get_ylim()[0] - 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                f"n={n:,}", ha='center', va='top', fontsize=10, color='gray')
    
    plt.tight_layout()
    save_figure(fig, figdir, filename_base)
    plt.close(fig)


# ============================== VIOLIN PLOTS BY DISEASE STATE ================
print("\n" + "=" * 84)
print("Generating violin plots by disease state...")
print("=" * 84)

# --- CytoTRACE2 Score Violin ---
save_violin_plot(
    adata, 
    score_key="cytotrace2_score",
    groupby="disease_state",
    order=disease_order,
    palette=disease_colors,
    title="CytoTRACE2 Score Distribution by Disease State",
    filename_base="violin_cytotrace2_score_by_disease_state",
    figdir=FIGDIR
)

# --- BCR Score Violin ---
save_violin_plot(
    adata,
    score_key="bcr_score",
    groupby="disease_state",
    order=disease_order,
    palette=disease_colors,
    title="BCR Signaling Score Distribution by Disease State",
    filename_base="violin_bcr_score_by_disease_state",
    figdir=FIGDIR
)

# --- OXPHOS Score Violin ---
save_violin_plot(
    adata,
    score_key="oxphos_score",
    groupby="disease_state",
    order=disease_order,
    palette=disease_colors,
    title="OXPHOS Score Distribution by Disease State",
    filename_base="violin_oxphos_score_by_disease_state",
    figdir=FIGDIR
)

# ============================== VIOLIN PLOTS BY SPECIES ======================
print("\n" + "=" * 84)
print("Generating violin plots by species...")
print("=" * 84)

# --- CytoTRACE2 Score by Species ---
save_violin_plot(
    adata, 
    score_key="cytotrace2_score",
    groupby="species",
    order=species_order,
    palette=species_colors,
    title="CytoTRACE2 Score Distribution by Species",
    filename_base="violin_cytotrace2_score_by_species",
    figdir=FIGDIR
)

# --- BCR Score by Species ---
save_violin_plot(
    adata,
    score_key="bcr_score",
    groupby="species",
    order=species_order,
    palette=species_colors,
    title="BCR Signaling Score Distribution by Species",
    filename_base="violin_bcr_score_by_species",
    figdir=FIGDIR
)

# --- OXPHOS Score by Species ---
save_violin_plot(
    adata,
    score_key="oxphos_score",
    groupby="species",
    order=species_order,
    palette=species_colors,
    title="OXPHOS Score Distribution by Species",
    filename_base="violin_oxphos_score_by_species",
    figdir=FIGDIR
)

# ============================== COMBINED MULTI-PANEL FIGURES =================
print("\n" + "=" * 84)
print("Generating combined multi-panel figures...")
print("=" * 84)

scores_to_plot = []
if "cytotrace2_score" in adata.obs.columns:
    scores_to_plot.append(("cytotrace2_score", "CytoTRACE2 Score"))
if "bcr_score" in adata.obs.columns:
    scores_to_plot.append(("bcr_score", "BCR Score"))
if "oxphos_score" in adata.obs.columns:
    scores_to_plot.append(("oxphos_score", "OXPHOS Score"))

# Combined by disease_state
if scores_to_plot:
    fig, axes = plt.subplots(1, len(scores_to_plot), figsize=(6*len(scores_to_plot), 7))
    if len(scores_to_plot) == 1:
        axes = [axes]
    
    for ax, (score_key, label) in zip(axes, scores_to_plot):
        df = adata.obs[[score_key, 'disease_state']].dropna()
        df[score_key] = df[score_key].astype(float)
        present_order = [c for c in disease_order if c in df['disease_state'].unique()]
        
        sns.violinplot(
            data=df,
            x='disease_state',
            y=score_key,
            order=present_order,
            palette=[disease_colors.get(c, "#999999") for c in present_order],
            inner="box",
            ax=ax,
            cut=0,
            scale="width"
        )
        ax.set_xlabel("")
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.tick_params(axis='x', rotation=45, labelsize=10)
    
    plt.suptitle("Score Distributions by Disease State", fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_figure(fig, FIGDIR, "violin_all_scores_by_disease_state_combined")
    plt.close(fig)

# Combined by species
if scores_to_plot:
    fig, axes = plt.subplots(1, len(scores_to_plot), figsize=(5*len(scores_to_plot), 7))
    if len(scores_to_plot) == 1:
        axes = [axes]
    
    for ax, (score_key, label) in zip(axes, scores_to_plot):
        df = adata.obs[[score_key, 'species']].dropna()
        df[score_key] = df[score_key].astype(float)
        present_order = [c for c in species_order if c in df['species'].unique()]
        
        sns.violinplot(
            data=df,
            x='species',
            y=score_key,
            order=present_order,
            palette=[species_colors.get(c, "#999999") for c in present_order],
            inner="box",
            ax=ax,
            cut=0,
            scale="width"
        )
        ax.set_xlabel("")
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.tick_params(axis='x', rotation=0, labelsize=11)
    
    plt.suptitle("Score Distributions by Species", fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_figure(fig, FIGDIR, "violin_all_scores_by_species_combined")
    plt.close(fig)

# ============================== LEIDEN CLUSTER CELL TYPE COMPOSITION =========
print("\n" + "=" * 84)
print("Analyzing Leiden cluster cell type composition (resolution 1.0)...")
print("=" * 84)

leiden_key = "leiden_1.0"
if leiden_key not in adata.obs.columns:
    print(f"  (warn) '{leiden_key}' not found, trying 'leiden'...")
    leiden_key = "leiden" if "leiden" in adata.obs.columns else None

celltype_key = "geneformer_predicted_celltype"
confidence_key = "geneformer_confidence"
CONFIDENCE_THRESHOLD = 0.8

if leiden_key and celltype_key in adata.obs.columns:
    # Filter to cells with celltype annotation AND high confidence (>0.8)
    has_annotation = adata.obs[celltype_key].notna()
    
    if confidence_key in adata.obs.columns:
        # Convert confidence to numeric and filter
        adata.obs[confidence_key] = pd.to_numeric(adata.obs[confidence_key], errors='coerce')
        high_confidence = adata.obs[confidence_key] > CONFIDENCE_THRESHOLD
        mask = has_annotation & high_confidence
        adata_annotated = adata[mask].copy()
        print(f"  Cells with Geneformer annotation: {has_annotation.sum():,}")
        print(f"  Cells with confidence > {CONFIDENCE_THRESHOLD}: {adata_annotated.n_obs:,}")
    else:
        adata_annotated = adata[has_annotation].copy()
        print(f"  Cells with Geneformer annotation: {adata_annotated.n_obs:,}")
        print(f"  (warn) No confidence scores found, using all annotated cells")
    
    # Get cluster composition by cell type
    cluster_celltype = adata_annotated.obs.groupby([leiden_key, celltype_key]).size().unstack(fill_value=0)
    cluster_celltype_pct = cluster_celltype.div(cluster_celltype.sum(axis=1), axis=0) * 100
    
    # Find majority cell type per cluster
    majority_celltype = cluster_celltype.idxmax(axis=1)
    majority_pct = cluster_celltype.max(axis=1) / cluster_celltype.sum(axis=1) * 100
    
    # Create summary DataFrame
    cluster_summary = pd.DataFrame({
        'cluster': majority_celltype.index,
        'majority_celltype': majority_celltype.values,
        'majority_pct': majority_pct.values,
        'n_cells': cluster_celltype.sum(axis=1).values
    })
    cluster_summary = cluster_summary.sort_values('cluster', key=lambda x: x.astype(int))
    
    # Save summary
    cluster_summary.to_csv(FIGDIR / "leiden_1.0_celltype_composition.csv", index=False)
    print(f"  ✓ leiden_1.0_celltype_composition.csv")
    print("\n  Cluster Cell Type Summary:")
    print(cluster_summary.to_string(index=False))
    
    # Generate color palette for ALL cell types (not just high-confidence subset)
    all_celltypes = adata.obs[celltype_key].dropna().unique().tolist()
    n_types = len(all_celltypes)
    cmap = plt.cm.get_cmap('tab20', max(n_types, 20))
    celltype_colors = {ct: cmap(i % 20) for i, ct in enumerate(sorted(all_celltypes))}
    
    # --- Plot: Stacked bar chart of cluster composition by cell type ---
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Sort clusters numerically
    cluster_order = sorted(cluster_celltype_pct.index, key=lambda x: int(x))
    cluster_celltype_pct = cluster_celltype_pct.loc[cluster_order]
    
    cluster_celltype_pct.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[celltype_colors.get(c, "#999999") for c in cluster_celltype_pct.columns],
        edgecolor='white',
        linewidth=0.5
    )
    
    ax.set_xlabel("Leiden Cluster (res=1.0)", fontsize=14, fontweight='bold')
    ax.set_ylabel("Percentage of Cells", fontsize=14, fontweight='bold')
    ax.set_title("Cell Type Composition per Leiden Cluster (Geneformer)", fontsize=16, fontweight='bold')
    ax.legend(title="Cell Type", bbox_to_anchor=(1.02, 1), loc='upper left', framealpha=0.9, fontsize=8)
    ax.tick_params(axis='x', rotation=0, labelsize=10)
    ax.set_ylim(0, 100)
    
    # Add cell count labels on top
    for i, clust in enumerate(cluster_order):
        n = cluster_celltype.loc[clust].sum()
        ax.text(i, 102, f"{int(n)}", ha='center', va='bottom', fontsize=8, rotation=90)
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "leiden_1.0_celltype_composition_stacked")
    plt.close(fig)
    
    # --- Plot: UMAP colored by majority cell type per cluster ---
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Create a new column with cluster annotated by majority cell type
    adata.obs['cluster_majority_celltype'] = adata.obs[leiden_key].map(
        lambda x: f"C{x}: {majority_celltype[x]}" if x in majority_celltype.index else f"C{x}: Unknown"
    )
    
    # Generate colors based on majority cell type
    unique_clusters = sorted(adata.obs['cluster_majority_celltype'].unique(), 
                              key=lambda x: int(x.split(':')[0].replace('C', '')))
    cluster_palette = {}
    for clust_label in unique_clusters:
        cluster_num = clust_label.split(':')[0].replace('C', '')
        if cluster_num in majority_celltype.index:
            ct = majority_celltype[cluster_num]
            cluster_palette[clust_label] = celltype_colors.get(ct, "#999999")
        else:
            cluster_palette[clust_label] = "#999999"
    
    sc.pl.umap(
        adata,
        color='cluster_majority_celltype',
        palette=cluster_palette,
        ax=ax,
        show=False,
        frameon=False,
        title="Leiden Clusters (res=1.0) colored by Majority Cell Type (Geneformer)",
        legend_loc='right margin',
        legend_fontsize=8,
        s=15
    )
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "umap_leiden_1.0_majority_celltype")
    plt.close(fig)
    
    # Clean up temp column
    adata.obs.drop(columns=['cluster_majority_celltype'], inplace=True)
    
    # --- Plot: UMAP colored by cell type directly ---
    fig, ax = plt.subplots(figsize=(14, 10))
    
    sc.pl.umap(
        adata,
        color=celltype_key,
        palette=celltype_colors,
        ax=ax,
        show=False,
        frameon=False,
        title="Geneformer Predicted Cell Type",
        legend_loc='right margin',
        legend_fontsize=8,
        s=15
    )
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "umap_celltype_geneformer")
    plt.close(fig)
    
    # --- Plot: Cluster composition by species ---
    cluster_species = adata_annotated.obs.groupby([leiden_key, 'species']).size().unstack(fill_value=0)
    cluster_species_pct = cluster_species.div(cluster_species.sum(axis=1), axis=0) * 100
    cluster_species_pct = cluster_species_pct.loc[cluster_order]
    
    fig, ax = plt.subplots(figsize=(14, 6))
    cluster_species_pct.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[species_colors.get(c, "#999999") for c in cluster_species_pct.columns],
        edgecolor='white',
        linewidth=0.5
    )
    
    ax.set_xlabel("Leiden Cluster (res=1.0)", fontsize=14, fontweight='bold')
    ax.set_ylabel("Percentage of Cells", fontsize=14, fontweight='bold')
    ax.set_title("Species Composition per Leiden Cluster", fontsize=16, fontweight='bold')
    ax.legend(title="Species", loc='upper right', framealpha=0.9)
    ax.tick_params(axis='x', rotation=0, labelsize=10)
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "leiden_1.0_species_composition_stacked")
    plt.close(fig)
    
    # --- Plot: Cluster composition by disease_state ---
    cluster_disease = adata_annotated.obs.groupby([leiden_key, 'disease_state']).size().unstack(fill_value=0)
    cluster_disease_pct = cluster_disease.div(cluster_disease.sum(axis=1), axis=0) * 100
    cluster_disease_pct = cluster_disease_pct.loc[cluster_order]
    
    fig, ax = plt.subplots(figsize=(14, 6))
    cluster_disease_pct.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[disease_colors.get(c, "#999999") for c in cluster_disease_pct.columns],
        edgecolor='white',
        linewidth=0.5
    )
    
    ax.set_xlabel("Leiden Cluster (res=1.0)", fontsize=14, fontweight='bold')
    ax.set_ylabel("Percentage of Cells", fontsize=14, fontweight='bold')
    ax.set_title("Disease State Composition per Leiden Cluster", fontsize=16, fontweight='bold')
    ax.legend(title="Disease State", loc='upper right', framealpha=0.9)
    ax.tick_params(axis='x', rotation=0, labelsize=10)
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "leiden_1.0_disease_state_composition_stacked")
    plt.close(fig)

else:
    if not leiden_key:
        print("  (warn) No Leiden clustering found in adata.obs")
    if celltype_key not in adata.obs.columns:
        print(f"  (warn) '{celltype_key}' not found in adata.obs")
        print(f"  Run geneformer_predict_human_mouse_integration.py first!")
        print(f"  Available columns: {list(adata.obs.columns)}")

# ============================== STATISTICAL SUMMARY ==========================
print("\n" + "=" * 84)
print("Statistical Summary")
print("=" * 84)

# By disease_state
summary_stats = []
for score_key, score_label in scores_to_plot:
    if score_key in adata.obs.columns:
        for ds in disease_order:
            if ds in adata.obs['disease_state'].unique():
                vals = adata.obs.loc[adata.obs['disease_state'] == ds, score_key].dropna().astype(float)
                summary_stats.append({
                    'score': score_label,
                    'disease_state': ds,
                    'n_cells': len(vals),
                    'mean': vals.mean(),
                    'median': vals.median(),
                    'std': vals.std(),
                    'min': vals.min(),
                    'max': vals.max()
                })

if summary_stats:
    stats_df = pd.DataFrame(summary_stats)
    stats_df.to_csv(FIGDIR / "score_statistics_by_disease_state.csv", index=False)
    print(f"  ✓ score_statistics_by_disease_state.csv")
    print("\n  By Disease State:")
    print(stats_df.to_string(index=False))

# By species
summary_stats_species = []
for score_key, score_label in scores_to_plot:
    if score_key in adata.obs.columns:
        for sp in species_order:
            if sp in adata.obs['species'].unique():
                vals = adata.obs.loc[adata.obs['species'] == sp, score_key].dropna().astype(float)
                summary_stats_species.append({
                    'score': score_label,
                    'species': sp,
                    'n_cells': len(vals),
                    'mean': vals.mean(),
                    'median': vals.median(),
                    'std': vals.std(),
                    'min': vals.min(),
                    'max': vals.max()
                })

if summary_stats_species:
    stats_df_sp = pd.DataFrame(summary_stats_species)
    stats_df_sp.to_csv(FIGDIR / "score_statistics_by_species.csv", index=False)
    print(f"\n  ✓ score_statistics_by_species.csv")
    print("\n  By Species:")
    print(stats_df_sp.to_string(index=False))

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  All figures saved to: {FIGDIR}")
print("\nDONE.\n")



