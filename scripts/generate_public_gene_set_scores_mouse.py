#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate UMAP and Violin Plots with Public/Database BCR and OXPHOS Gene Sets
============================================================================

This script:
1. Loads the integrated mouse data
2. Computes BCR and OXPHOS scores using public gene sets (MSigDB, Reactome, KEGG)
3. Generates UMAP plots and violin plots for comparison with manually curated sets

Usage:
    conda activate <your_scanpy_env>  # or appropriate environment
    python generate_public_gene_set_scores_mouse.py

Author: J
Date: 2025-01-XX
"""

import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
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
INPUT_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/figures")
FIGDIR_PUBLIC = OUTDIR / "public_gene_sets"
FIGDIR_PUBLIC.mkdir(parents=True, exist_ok=True)

# ============================== PUBLIC GENE SETS =============================
# MSigDB Hallmark: BCR Signaling (mouse gene symbols, lowercase)
# Based on HALLMARK_BCR_SIGNALING_PATHWAY
MSIGDB_BCR_MOUSE = [
    'Cd79a', 'Cd79b', 'Cd19', 'Cd22', 'Cd72', 'Cr2', 'Ms4a1', 'Cd40', 'Cd40lg',
    'Btk', 'Lyn', 'Syk', 'Blk', 'Blnk', 'Pik3cd', 'Pik3ap1', 'Pik3r1', 'Pik3r2',
    'Plcg2', 'Prkcb', 'Prkca', 'Plcg1', 'Nfkb1', 'Nfkb2', 'Rel', 'Rela', 'Relb',
    'Nfatc1', 'Nfatc2', 'Nfatc3', 'Bcl10', 'Card11', 'Malt1', 'Map3k7', 'Ikbkb',
    'Ikbkg', 'Chuk', 'Ptpn6', 'Ptprc', 'Vav1', 'Vav2', 'Vav3', 'Grb2', 'Sos1',
    'Sos2', 'Hras', 'Kras', 'Nras', 'Raf1', 'Map2k1', 'Map2k2', 'Mapk1', 'Mapk3',
    'Mapk8', 'Mapk9', 'Mapk14', 'Fos', 'Jun', 'Junb', 'Jund', 'Egr1', 'Egr2',
    'Egr3', 'Ighm', 'Ighd', 'Ighg1', 'Ighg2a', 'Ighg2b', 'Ighg2c', 'Ighg3', 'Ighe',
    'Igha', 'Cd81', 'Cd82', 'Cd86', 'Cd80', 'Il4', 'Il4ra', 'Il13', 'Il13ra1'
]

# Reactome: B Cell Receptor Signaling (mouse)
REACTOME_BCR_MOUSE = [
    'Cd79a', 'Cd79b', 'Cd19', 'Cd22', 'Cd72', 'Cr2', 'Ms4a1', 'Btk', 'Lyn', 'Syk',
    'Blk', 'Blnk', 'Pik3cd', 'Pik3ap1', 'Pik3r1', 'Pik3r2', 'Plcg2', 'Prkcb',
    'Nfkb1', 'Nfkb2', 'Rel', 'Rela', 'Nfatc1', 'Nfatc2', 'Bcl10', 'Card11', 'Malt1',
    'Map3k7', 'Ikbkb', 'Ikbkg', 'Chuk', 'Ptpn6', 'Ptprc', 'Vav1', 'Vav2', 'Vav3',
    'Grb2', 'Sos1', 'Sos2', 'Hras', 'Kras', 'Nras', 'Raf1', 'Map2k1', 'Map2k2',
    'Mapk1', 'Mapk3', 'Ighm', 'Ighd', 'Ighg1', 'Ighg2a', 'Ighg2b', 'Ighg2c', 'Ighg3',
    'Ighe', 'Igha', 'Cd40', 'Cd40lg', 'Tnfrsf13b', 'Tnfrsf13c', 'Tnfrsf17'
]

# MSigDB Hallmark: OXPHOS (mouse gene symbols, lowercase)
# Based on HALLMARK_OXIDATIVE_PHOSPHORYLATION
MSIGDB_OXPHOS_MOUSE = [
    'Cox4i1', 'Cox5a', 'Cox5b', 'Cox6a1', 'Cox6b1', 'Cox6c', 'Cox7a1', 'Cox7a2',
    'Cox7b', 'Cox7c', 'Cox8a', 'Cox8b', 'Cox8c', 'Cyc1', 'Cycs', 'Ndufa1', 'Ndufa2',
    'Ndufa3', 'Ndufa4', 'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9', 'Ndufa10',
    'Ndufa11', 'Ndufa12', 'Ndufa13', 'Ndufab1', 'Ndufb1', 'Ndufb2', 'Ndufb3', 'Ndufb4',
    'Ndufb5', 'Ndufb6', 'Ndufb7', 'Ndufb8', 'Ndufb9', 'Ndufb10', 'Ndufb11', 'Ndufc1',
    'Ndufc2', 'Ndufs1', 'Ndufs2', 'Ndufs3', 'Ndufs4', 'Ndufs5', 'Ndufs6', 'Ndufs7',
    'Ndufs8', 'Ndufv1', 'Ndufv2', 'Ndufv3', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Sdhaf1',
    'Sdhaf2', 'Uqcr10', 'Uqcr11', 'Uqcrb', 'Uqcrc1', 'Uqcrc2', 'Uqcrfs1', 'Uqcrh',
    'Uqcrq', 'Atp5f1a', 'Atp5f1b', 'Atp5f1c', 'Atp5f1d', 'Atp5f1e', 'Atp5mc1', 'Atp5mc2',
    'Atp5mc3', 'Atp5me', 'Atp5mf', 'Atp5mg', 'Atp5pb', 'Atp5pd', 'Atp5pf', 'Atp5po',
    'Atp5if1', 'Atp5j', 'Atp5j2', 'Atp5l', 'Atp5o', 'Atp5s', 'Atp6v0a1', 'Atp6v0a2',
    'Atp6v0a4', 'Atp6v0b', 'Atp6v0c', 'Atp6v0d1', 'Atp6v0d2', 'Atp6v0e1', 'Atp6v0e2',
    'Atp6v1a', 'Atp6v1b1', 'Atp6v1b2', 'Atp6v1c1', 'Atp6v1c2', 'Atp6v1d', 'Atp6v1e1',
    'Atp6v1e2', 'Atp6v1f', 'Atp6v1g1', 'Atp6v1g2', 'Atp6v1g3', 'Atp6v1h'
]

# Reactome: Respiratory Electron Transport (mouse)
REACTOME_OXPHOS_MOUSE = [
    'Cox4i1', 'Cox5a', 'Cox5b', 'Cox6a1', 'Cox6b1', 'Cox6c', 'Cox7a2', 'Cox7b',
    'Cox7c', 'Cox8a', 'Cyc1', 'Cycs', 'Ndufa1', 'Ndufa2', 'Ndufa3', 'Ndufa4',
    'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9', 'Ndufa10', 'Ndufa11', 'Ndufa12',
    'Ndufa13', 'Ndufab1', 'Ndufb1', 'Ndufb2', 'Ndufb3', 'Ndufb4', 'Ndufb5', 'Ndufb6',
    'Ndufb7', 'Ndufb8', 'Ndufb9', 'Ndufb10', 'Ndufb11', 'Ndufc1', 'Ndufc2', 'Ndufs1',
    'Ndufs2', 'Ndufs3', 'Ndufs4', 'Ndufs5', 'Ndufs6', 'Ndufs7', 'Ndufs8', 'Ndufv1',
    'Ndufv2', 'Ndufv3', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Uqcr10', 'Uqcr11', 'Uqcrb',
    'Uqcrc1', 'Uqcrc2', 'Uqcrfs1', 'Uqcrh', 'Uqcrq', 'Atp5f1a', 'Atp5f1b', 'Atp5f1c',
    'Atp5f1d', 'Atp5f1e', 'Atp5mc1', 'Atp5mc2', 'Atp5mc3', 'Atp5me', 'Atp5mf',
    'Atp5mg', 'Atp5pb', 'Atp5pd', 'Atp5pf', 'Atp5po'
]

# Fetch actual KEGG gene sets using gseapy
def fetch_kegg_gene_sets():
    """Fetch actual KEGG pathway gene sets for mouse."""
    try:
        import gseapy as gp
        
        # Get KEGG gene set library (correct API: name first, then organism)
        gs = gp.get_library(name='KEGG_2019_Mouse', organism='Mouse')
        
        # Find BCR and OXPHOS pathways by searching keys
        bcr_pathway = None
        oxphos_pathway = None
        
        for pathway_name in gs.keys():
            pathway_lower = pathway_name.lower()
            # Match exact KEGG pathway names
            if pathway_lower == 'b cell receptor signaling pathway':
                bcr_pathway = pathway_name
            elif pathway_lower == 'oxidative phosphorylation':
                oxphos_pathway = pathway_name
        
        print(f"  Found KEGG pathways:")
        if bcr_pathway:
            print(f"    BCR: {bcr_pathway}")
        if oxphos_pathway:
            print(f"    OXPHOS: {oxphos_pathway}")
        
        # Extract gene sets
        kegg_genesets = {}
        
        if bcr_pathway and bcr_pathway in gs:
            # Convert to lowercase and handle gene name format
            kegg_genesets['bcr'] = [str(g).lower() for g in gs[bcr_pathway] if g]
            print(f"    BCR genes: {len(kegg_genesets['bcr'])} (sample: {kegg_genesets['bcr'][:5]})")
        
        if oxphos_pathway and oxphos_pathway in gs:
            # Convert to lowercase and handle gene name format
            kegg_genesets['oxphos'] = [str(g).lower() for g in gs[oxphos_pathway] if g]
            print(f"    OXPHOS genes: {len(kegg_genesets['oxphos'])} (sample: {kegg_genesets['oxphos'][:5]})")
        
        return kegg_genesets if kegg_genesets else None
        
    except ImportError:
        print("  WARNING: gseapy not available, using manually curated KEGG gene sets")
        return None
    except Exception as e:
        print(f"  WARNING: Error fetching KEGG gene sets: {e}")
        import traceback
        traceback.print_exc()
        print("  Using manually curated KEGG gene sets")
        return None

# Fallback: Manually curated KEGG gene sets (used if fetch fails)
# KEGG: B Cell Receptor Signaling Pathway (mouse) - mmu04662
KEGG_BCR_MOUSE_FALLBACK = [
    'Cd79a', 'Cd79b', 'Cd19', 'Cd22', 'Cd72', 'Cr2', 'Ms4a1', 'Btk', 'Lyn', 'Syk',
    'Blk', 'Blnk', 'Pik3cd', 'Pik3ap1', 'Pik3r1', 'Pik3r2', 'Pik3r3', 'Plcg2', 'Prkcb',
    'Nfkb1', 'Nfkb2', 'Rel', 'Rela', 'Nfatc1', 'Nfatc2', 'Bcl10', 'Card11', 'Malt1',
    'Map3k7', 'Ikbkb', 'Ikbkg', 'Chuk', 'Ptpn6', 'Ptprc', 'Vav1', 'Vav2', 'Vav3',
    'Grb2', 'Sos1', 'Sos2', 'Hras', 'Kras', 'Nras', 'Raf1', 'Map2k1', 'Map2k2',
    'Mapk1', 'Mapk3', 'Ighm', 'Ighd', 'Ighg1', 'Ighg2a', 'Ighg2b', 'Ighg2c', 'Ighg3',
    'Cd40', 'Cd40lg', 'Tnfrsf13b', 'Tnfrsf13c', 'Tnfrsf17', 'Tnf', 'Tnfrsf1a', 'Tnfrsf1b'
]

# KEGG: Oxidative Phosphorylation (mouse) - mmu00190
KEGG_OXPHOS_MOUSE_FALLBACK = [
    'Cox4i1', 'Cox5a', 'Cox5b', 'Cox6a1', 'Cox6b1', 'Cox6c', 'Cox7a2', 'Cox7b',
    'Cox7c', 'Cox8a', 'Cox8b', 'Cox8c', 'Cyc1', 'Cycs', 'Ndufa1', 'Ndufa2',
    'Ndufa3', 'Ndufa4', 'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9', 'Ndufa10',
    'Ndufa11', 'Ndufa12', 'Ndufa13', 'Ndufab1', 'Ndufb1', 'Ndufb2', 'Ndufb3', 'Ndufb4',
    'Ndufb5', 'Ndufb6', 'Ndufb7', 'Ndufb8', 'Ndufb9', 'Ndufb10', 'Ndufb11', 'Ndufc1',
    'Ndufc2', 'Ndufs1', 'Ndufs2', 'Ndufs3', 'Ndufs4', 'Ndufs5', 'Ndufs6', 'Ndufs7',
    'Ndufs8', 'Ndufv1', 'Ndufv2', 'Ndufv3', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Sdhaf1',
    'Sdhaf2', 'Uqcr10', 'Uqcr11', 'Uqcrb', 'Uqcrc1', 'Uqcrc2', 'Uqcrfs1', 'Uqcrh',
    'Uqcrq', 'Atp5f1a', 'Atp5f1b', 'Atp5f1c', 'Atp5f1d', 'Atp5f1e', 'Atp5mc1', 'Atp5mc2',
    'Atp5mc3', 'Atp5me', 'Atp5mf', 'Atp5mg', 'Atp5pb', 'Atp5pd', 'Atp5pf', 'Atp5po',
    'Atp5if1', 'Atp5j', 'Atp5j2', 'Atp5l', 'Atp5o', 'Atp5s'
]

# ============================== HELPER FUNCTIONS =============================
def save_figure(fig, basename, dpi=300):
    """Save figure in PNG, PDF, and SVG formats."""
    for fmt in ['png', 'pdf', 'svg']:
        filepath = FIGDIR_PUBLIC / f'{basename}.{fmt}'
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
    print(f"  ✓ {basename} (.png, .pdf, .svg)")

def compute_and_plot_scores(adata, gene_set_name, gene_list, score_name, title_suffix=""):
    """Compute gene set score and generate UMAP and violin plots."""
    print(f"\n  Processing {gene_set_name}...")
    
    # Filter to available genes - handle case-insensitive matching
    # Create lowercase mapping for case-insensitive lookup
    var_names_lower = {str(g).lower(): str(g) for g in adata.var_names}
    avail = set(adata.var_names)
    
    # Try exact match first, then case-insensitive
    genes_present = []
    for g in gene_list:
        g_str = str(g)
        if g_str in avail:
            genes_present.append(g_str)
        elif g_str.lower() in var_names_lower:
            genes_present.append(var_names_lower[g_str.lower()])
    
    if len(genes_present) < 5:
        print(f"    WARNING: Only {len(genes_present)}/{len(gene_list)} genes found, skipping...")
        return None
    
    print(f"    Found {len(genes_present)}/{len(gene_list)} genes")
    
    # Compute score
    sc.tl.score_genes(adata, gene_list=genes_present, score_name=score_name)
    
    # UMAP plot
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color=score_name, ax=ax, show=False, frameon=False,
               cmap="RdYlBu_r", s=20, title=f"{gene_set_name} {title_suffix}")
    save_figure(fig, f"umap_{score_name}")
    plt.close()
    
    # Violin plot by condition
    if 'condition' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        condition_order = sorted(adata.obs['condition'].unique())
        sns.violinplot(data=adata.obs, x='condition', y=score_name,
                       order=condition_order, ax=ax, inner='box', cut=0)
        ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{gene_set_name} Score', fontsize=12, fontweight='bold')
        ax.set_title(f'{gene_set_name} Score by Condition {title_suffix}', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_figure(fig, f"violin_{score_name}_by_condition")
        plt.close()
    
    # Violin plot by sample
    if 'sample_id' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(14, 6))
        sample_order = sorted(adata.obs['sample_id'].unique())
        sns.violinplot(data=adata.obs, x='sample_id', y=score_name,
                       order=sample_order, ax=ax, inner='box', cut=0)
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{gene_set_name} Score', fontsize=12, fontweight='bold')
        ax.set_title(f'{gene_set_name} Score by Sample {title_suffix}', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_figure(fig, f"violin_{score_name}_by_sample")
        plt.close()
    
    return score_name

# ============================== MAIN PIPELINE ================================
print("=" * 84)
print("PUBLIC GENE SET SCORES: BCR and OXPHOS")
print("=" * 84)
print(f"Input: {INPUT_H5AD}")
print(f"Output: {FIGDIR_PUBLIC}\n")

# Fetch KEGG gene sets first
print("Fetching KEGG gene sets from database...")
KEGG_GENESETS = fetch_kegg_gene_sets()

# Use fetched KEGG sets if available, otherwise use fallback
if KEGG_GENESETS and 'bcr' in KEGG_GENESETS:
    KEGG_BCR_MOUSE = KEGG_GENESETS['bcr']
    print(f"  Using KEGG BCR: {len(KEGG_BCR_MOUSE)} genes from database")
else:
    KEGG_BCR_MOUSE = KEGG_BCR_MOUSE_FALLBACK
    print(f"  Using manually curated KEGG BCR: {len(KEGG_BCR_MOUSE)} genes")

if KEGG_GENESETS and 'oxphos' in KEGG_GENESETS:
    KEGG_OXPHOS_MOUSE = KEGG_GENESETS['oxphos']
    print(f"  Using KEGG OXPHOS: {len(KEGG_OXPHOS_MOUSE)} genes from database")
else:
    KEGG_OXPHOS_MOUSE = KEGG_OXPHOS_MOUSE_FALLBACK
    print(f"  Using manually curated KEGG OXPHOS: {len(KEGG_OXPHOS_MOUSE)} genes")
print()

# Load data
print("Loading integrated mouse data...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Check UMAP exists
if "X_umap" not in adata.obsm:
    print("  WARNING: X_umap not found. Computing from X_scvi...")
    if "X_scvi" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=30)
        sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    else:
        raise ValueError("Neither X_umap nor X_scvi found")

# BCR scores from public gene sets
print("\n" + "=" * 84)
print("BCR SIGNALING SCORES (Public Gene Sets)")
print("=" * 84)

bcr_scores = {}
bcr_scores['msigdb'] = compute_and_plot_scores(
    adata, "MSigDB Hallmark BCR", MSIGDB_BCR_MOUSE, 
    "bcr_score_msigdb", "(MSigDB Hallmark)"
)
bcr_scores['reactome'] = compute_and_plot_scores(
    adata, "Reactome BCR", REACTOME_BCR_MOUSE,
    "bcr_score_reactome", "(Reactome)"
)
bcr_scores['kegg'] = compute_and_plot_scores(
    adata, "KEGG BCR", KEGG_BCR_MOUSE,
    "bcr_score_kegg", "(KEGG)"
)

# OXPHOS scores from public gene sets
print("\n" + "=" * 84)
print("OXPHOS SCORES (Public Gene Sets)")
print("=" * 84)

oxphos_scores = {}
oxphos_scores['msigdb'] = compute_and_plot_scores(
    adata, "MSigDB Hallmark OXPHOS", MSIGDB_OXPHOS_MOUSE,
    "oxphos_score_msigdb", "(MSigDB Hallmark)"
)
oxphos_scores['reactome'] = compute_and_plot_scores(
    adata, "Reactome OXPHOS", REACTOME_OXPHOS_MOUSE,
    "oxphos_score_reactome", "(Reactome)"
)
oxphos_scores['kegg'] = compute_and_plot_scores(
    adata, "KEGG OXPHOS", KEGG_OXPHOS_MOUSE,
    "oxphos_score_kegg", "(KEGG)"
)

# Comparison plots (if both custom and public scores exist)
print("\n" + "=" * 84)
print("COMPARISON PLOTS")
print("=" * 84)

if 'bcr_score' in adata.obs.columns and any(bcr_scores.values()):
    print("\n  Creating BCR comparison plots...")
    # Side-by-side UMAP comparison
    n_public = sum(1 for v in bcr_scores.values() if v is not None)
    if n_public > 0:
        fig, axes = plt.subplots(1, n_public + 1, figsize=(6*(n_public+1), 10))
        if n_public == 0:
            axes = [axes]
        
        # Custom score
        sc.pl.umap(adata, color='bcr_score', ax=axes[0], show=False, frameon=False,
                   cmap="RdYlBu_r", s=20, title="BCR Score (Custom)")
        
        # Public scores
        idx = 1
        for name, score_key in [('msigdb', 'bcr_score_msigdb'), 
                                ('reactome', 'bcr_score_reactome'),
                                ('kegg', 'bcr_score_kegg')]:
            if score_key in adata.obs.columns:
                sc.pl.umap(adata, color=score_key, ax=axes[idx], show=False, frameon=False,
                           cmap="RdYlBu_r", s=20, title=f"BCR Score ({name.upper()})")
                idx += 1
        
        plt.tight_layout()
        save_figure(fig, "umap_bcr_comparison")
        plt.close()

if 'oxphos_score' in adata.obs.columns and any(oxphos_scores.values()):
    print("\n  Creating OXPHOS comparison plots...")
    # Side-by-side UMAP comparison
    n_public = sum(1 for v in oxphos_scores.values() if v is not None)
    if n_public > 0:
        fig, axes = plt.subplots(1, n_public + 1, figsize=(6*(n_public+1), 10))
        if n_public == 0:
            axes = [axes]
        
        # Custom score
        sc.pl.umap(adata, color='oxphos_score', ax=axes[0], show=False, frameon=False,
                   cmap="RdYlBu_r", s=20, title="OXPHOS Score (Custom)")
        
        # Public scores
        idx = 1
        for name, score_key in [('msigdb', 'oxphos_score_msigdb'),
                                ('reactome', 'oxphos_score_reactome'),
                                ('kegg', 'oxphos_score_kegg')]:
            if score_key in adata.obs.columns:
                sc.pl.umap(adata, color=score_key, ax=axes[idx], show=False, frameon=False,
                           cmap="RdYlBu_r", s=20, title=f"OXPHOS Score ({name.upper()})")
                idx += 1
        
        plt.tight_layout()
        save_figure(fig, "umap_oxphos_comparison")
        plt.close()

# Save updated adata
output_h5ad = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated_with_public_scores.h5ad")
adata.write_h5ad(output_h5ad)
print(f"\n✓ Saved updated data: {output_h5ad}")

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  All figures saved to: {FIGDIR_PUBLIC}")
print("\nDONE.\n")


