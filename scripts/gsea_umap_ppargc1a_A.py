#!/usr/bin/env python3
"""
GSEA-style Gene Set Enrichment on scVI UMAP
============================================
Scores cells for PPARGC1A target genes and MOOTHA PGC gene sets,
then visualizes enrichment scores on UMAP.

Gene sets are human - we convert to mouse orthologs (capitalize first letter).

Author: Generated for CytoTRACE2 analysis
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# ===== CONFIGURATION =====
INPUT_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
GMT_PPARGC1A = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/PPARGC1A_TARGET_GENES.v2023.1.Hs.gmt")
GMT_MOOTHA = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/Mootha PGC.gmt")
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/figures/gsea_ppargc1a")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output formats
SAVE_FORMATS = ['png', 'svg', 'pdf']


def parse_gmt(gmt_path: Path) -> dict:
    """Parse GMT file and return gene set dictionary."""
    gene_sets = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                name = parts[0]
                # parts[1] is description/URL
                genes = parts[2:]
                gene_sets[name] = genes
    return gene_sets


def human_to_mouse_genes(human_genes: list) -> list:
    """
    Convert human gene symbols to mouse orthologs.
    Mouse genes are typically: First letter uppercase, rest lowercase.
    e.g., HSPA1A -> Hspa1a, ATP5MC3 -> Atp5mc3
    """
    mouse_genes = []
    for gene in human_genes:
        if gene.startswith('ENSG') or gene.startswith('LINC') or gene.startswith('MIR'):
            # Skip Ensembl IDs, lncRNAs, and miRNAs (may not have simple orthologs)
            continue
        # Convert: GENE -> Gene (first letter cap, rest lowercase)
        mouse_gene = gene.capitalize()
        mouse_genes.append(mouse_gene)
    return mouse_genes


def save_figure(fig, output_path: Path):
    """Save figure in multiple formats."""
    for fmt in SAVE_FORMATS:
        out_file = output_path.with_suffix(f'.{fmt}')
        fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {output_path.stem} (.png, .svg, .pdf)")


def main():
    print("=" * 70)
    print("GSEA-style Enrichment Scoring on scVI UMAP")
    print("=" * 70)
    
    # Load data
    print(f"\n[1] Loading: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"    Cells: {adata.n_obs:,}")
    print(f"    Genes: {adata.n_vars:,}")
    
    # Parse gene sets
    print(f"\n[2] Parsing gene sets...")
    
    # PPARGC1A gene set
    ppargc1a_sets = parse_gmt(GMT_PPARGC1A)
    ppargc1a_genes_human = list(ppargc1a_sets.values())[0]
    ppargc1a_genes_mouse = human_to_mouse_genes(ppargc1a_genes_human)
    print(f"    PPARGC1A_TARGET_GENES: {len(ppargc1a_genes_human)} human genes")
    
    # MOOTHA_PGC gene set
    mootha_sets = parse_gmt(GMT_MOOTHA)
    mootha_genes_human = list(mootha_sets.values())[0]
    mootha_genes_mouse = human_to_mouse_genes(mootha_genes_human)
    print(f"    MOOTHA_PGC: {len(mootha_genes_human)} human genes")
    
    # Check which genes are present in the dataset
    print(f"\n[3] Checking gene overlap with dataset...")
    genes_in_data = set(adata.var_names)
    
    ppargc1a_found = [g for g in ppargc1a_genes_mouse if g in genes_in_data]
    mootha_found = [g for g in mootha_genes_mouse if g in genes_in_data]
    
    print(f"    PPARGC1A: {len(ppargc1a_found)}/{len(ppargc1a_genes_mouse)} mouse genes found")
    print(f"    MOOTHA_PGC: {len(mootha_found)}/{len(mootha_genes_mouse)} mouse genes found")
    
    # Save gene lists for reference
    gene_info = pd.DataFrame({
        'PPARGC1A_human': pd.Series(ppargc1a_genes_human),
        'PPARGC1A_mouse': pd.Series(ppargc1a_genes_mouse),
        'PPARGC1A_found': pd.Series(ppargc1a_found),
    })
    gene_info.to_csv(OUTPUT_DIR / "ppargc1a_genes.csv", index=False)
    
    gene_info2 = pd.DataFrame({
        'MOOTHA_human': pd.Series(mootha_genes_human),
        'MOOTHA_mouse': pd.Series(mootha_genes_mouse),
        'MOOTHA_found': pd.Series(mootha_found),
    })
    gene_info2.to_csv(OUTPUT_DIR / "mootha_pgc_genes.csv", index=False)
    print(f"    Gene lists saved to {OUTPUT_DIR}")
    
    # Score cells for each gene set
    print(f"\n[4] Scoring cells for gene sets...")
    
    # PPARGC1A score
    sc.tl.score_genes(adata, gene_list=ppargc1a_found, 
                      score_name='PPARGC1A_score', ctrl_size=100)
    print(f"    PPARGC1A_score: mean={adata.obs['PPARGC1A_score'].mean():.4f}, "
          f"std={adata.obs['PPARGC1A_score'].std():.4f}")
    
    # MOOTHA_PGC score
    sc.tl.score_genes(adata, gene_list=mootha_found, 
                      score_name='MOOTHA_PGC_score', ctrl_size=100)
    print(f"    MOOTHA_PGC_score: mean={adata.obs['MOOTHA_PGC_score'].mean():.4f}, "
          f"std={adata.obs['MOOTHA_PGC_score'].std():.4f}")
    
    # Set up plotting style
    sc.set_figure_params(dpi=150, fontsize=12, figsize=(6, 5))
    
    # ===== PLOT 1: PPARGC1A enrichment on UMAP =====
    print(f"\n[5] Generating UMAP plots...")
    
    print("    - PPARGC1A enrichment UMAP")
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color='PPARGC1A_score', ax=ax, show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='PPARGC1A Target Genes Enrichment')
    save_figure(fig, OUTPUT_DIR / "umap_PPARGC1A_enrichment")
    plt.close()
    
    # ===== PLOT 2: MOOTHA_PGC enrichment on UMAP =====
    print("    - MOOTHA_PGC enrichment UMAP")
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color='MOOTHA_PGC_score', ax=ax, show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='MOOTHA PGC Enrichment')
    save_figure(fig, OUTPUT_DIR / "umap_MOOTHA_PGC_enrichment")
    plt.close()
    
    # ===== PLOT 3: Both scores side by side =====
    print("    - Combined panel")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    sc.pl.umap(adata, color='PPARGC1A_score', ax=axes[0], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='PPARGC1A Target Genes')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score', ax=axes[1], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='MOOTHA PGC')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "umap_enrichment_combined")
    plt.close()
    
    # ===== PLOT 4: Enrichment scores with condition =====
    print("    - Enrichment with condition panel")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    sc.pl.umap(adata, color='condition', ax=axes[0, 0], show=False,
               frameon=False, title='Condition')
    
    sc.pl.umap(adata, color='leiden_1.0', ax=axes[0, 1], show=False,
               frameon=False, title='Leiden Clusters (res=1.0)')
    
    sc.pl.umap(adata, color='PPARGC1A_score', ax=axes[1, 0], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='PPARGC1A Target Genes Enrichment')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score', ax=axes[1, 1], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='MOOTHA PGC Enrichment')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "umap_enrichment_with_context")
    plt.close()
    
    # ===== PLOT 5: Split by condition =====
    print("    - Enrichment split by condition")
    conditions = sorted(adata.obs['condition'].unique())
    n_cond = len(conditions)
    
    for score_name in ['PPARGC1A_score', 'MOOTHA_PGC_score']:
        fig, axes = plt.subplots(1, n_cond, figsize=(4*n_cond, 4))
        if n_cond == 1:
            axes = [axes]
        
        # Get global vmin/vmax for consistent coloring
        vmin, vmax = adata.obs[score_name].quantile([0.01, 0.99])
        vabs = max(abs(vmin), abs(vmax))
        
        for ax, cond in zip(axes, conditions):
            mask = adata.obs['condition'] == cond
            sc.pl.umap(adata[mask], color=score_name, ax=ax, show=False,
                       cmap='RdBu_r', vmin=-vabs, vmax=vabs, frameon=False,
                       title=f'{cond}\n(n={mask.sum():,})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / f"umap_{score_name}_by_condition")
        plt.close()
    
    # ===== PLOT 6: Violin plots by condition =====
    print("    - Violin plots by condition")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sc.pl.violin(adata, 'PPARGC1A_score', groupby='condition', ax=axes[0], 
                 show=False, rotation=45)
    axes[0].set_title('PPARGC1A Target Genes Enrichment')
    axes[0].set_xlabel('')
    
    sc.pl.violin(adata, 'MOOTHA_PGC_score', groupby='condition', ax=axes[1], 
                 show=False, rotation=45)
    axes[1].set_title('MOOTHA PGC Enrichment')
    axes[1].set_xlabel('')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "violin_enrichment_by_condition")
    plt.close()
    
    # ===== PLOT 7: Violin plots by Leiden cluster =====
    print("    - Violin plots by cluster")
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    
    sc.pl.violin(adata, 'PPARGC1A_score', groupby='leiden_1.0', ax=axes[0], 
                 show=False, rotation=0)
    axes[0].set_title('PPARGC1A Target Genes Enrichment by Cluster')
    axes[0].set_xlabel('Leiden Cluster')
    
    sc.pl.violin(adata, 'MOOTHA_PGC_score', groupby='leiden_1.0', ax=axes[1], 
                 show=False, rotation=0)
    axes[1].set_title('MOOTHA PGC Enrichment by Cluster')
    axes[1].set_xlabel('Leiden Cluster')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "violin_enrichment_by_cluster")
    plt.close()
    
    # ===== PLOT 8: Correlation with CytoTRACE2 =====
    print("    - Correlation with CytoTRACE2")
    if 'cytotrace2_score' in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # PPARGC1A vs CytoTRACE2
        ax = axes[0]
        x = adata.obs['cytotrace2_score']
        y = adata.obs['PPARGC1A_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('CytoTRACE2 Score')
        ax.set_ylabel('PPARGC1A Enrichment')
        ax.set_title(f'PPARGC1A vs CytoTRACE2\n(r = {corr:.3f})')
        
        # MOOTHA vs CytoTRACE2
        ax = axes[1]
        y = adata.obs['MOOTHA_PGC_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='darkred')
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('CytoTRACE2 Score')
        ax.set_ylabel('MOOTHA PGC Enrichment')
        ax.set_title(f'MOOTHA PGC vs CytoTRACE2\n(r = {corr:.3f})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "scatter_enrichment_vs_cytotrace2")
        plt.close()
    
    # ===== PLOT 9: Comparison with OXPHOS if available =====
    if 'oxphos_score' in adata.obs.columns:
        print("    - Correlation with OXPHOS score")
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        x = adata.obs['oxphos_score']
        
        ax = axes[0]
        y = adata.obs['PPARGC1A_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('OXPHOS Score')
        ax.set_ylabel('PPARGC1A Enrichment')
        ax.set_title(f'PPARGC1A vs OXPHOS\n(r = {corr:.3f})')
        
        ax = axes[1]
        y = adata.obs['MOOTHA_PGC_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='darkred')
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('OXPHOS Score')
        ax.set_ylabel('MOOTHA PGC Enrichment')
        ax.set_title(f'MOOTHA PGC vs OXPHOS\n(r = {corr:.3f})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "scatter_enrichment_vs_oxphos")
        plt.close()
    
    # Save statistics
    print(f"\n[6] Saving statistics...")
    # Select only numeric columns for aggregation
    numeric_cols = ['PPARGC1A_score', 'MOOTHA_PGC_score']
    if 'cytotrace2_score' in adata.obs.columns:
        numeric_cols.append('cytotrace2_score')
    if 'oxphos_score' in adata.obs.columns:
        numeric_cols.append('oxphos_score')
    
    stats = adata.obs[numeric_cols].copy()
    stats['condition'] = adata.obs['condition'].astype(str)
    stats['leiden_1.0'] = adata.obs['leiden_1.0'].astype(str)
    
    # Summary by condition
    summary_cond = stats.groupby('condition')[numeric_cols].agg(['mean', 'std', 'median']).round(4)
    summary_cond.to_csv(OUTPUT_DIR / "enrichment_by_condition.csv")
    
    # Summary by cluster
    summary_clust = stats.groupby('leiden_1.0')[numeric_cols].agg(['mean', 'std', 'median']).round(4)
    summary_clust.to_csv(OUTPUT_DIR / "enrichment_by_cluster.csv")
    
    print(f"\n    Enrichment by condition:")
    print(summary_cond[['PPARGC1A_score', 'MOOTHA_PGC_score']].to_string())
    
    print("\n" + "=" * 70)
    print("✓ COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"\nFiles created:")
    print("  - umap_PPARGC1A_enrichment.png/svg/pdf")
    print("  - umap_MOOTHA_PGC_enrichment.png/svg/pdf")
    print("  - umap_enrichment_combined.png/svg/pdf")
    print("  - umap_enrichment_with_context.png/svg/pdf")
    print("  - umap_PPARGC1A_score_by_condition.png/svg/pdf")
    print("  - umap_MOOTHA_PGC_score_by_condition.png/svg/pdf")
    print("  - violin_enrichment_by_condition.png/svg/pdf")
    print("  - violin_enrichment_by_cluster.png/svg/pdf")
    print("  - scatter_enrichment_vs_cytotrace2.png/svg/pdf")
    print("  - scatter_enrichment_vs_oxphos.png/svg/pdf (if available)")
    print("  - enrichment_by_condition.csv")
    print("  - enrichment_by_cluster.csv")
    print("  - ppargc1a_genes.csv")
    print("  - mootha_pgc_genes.csv")


if __name__ == "__main__":
    main()



