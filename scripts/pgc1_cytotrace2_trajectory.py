#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PGC1A/B Expression Trajectory along CytoTRACE2 Pseudotime
==========================================================

Visualizes Ppargc1a (PGC1A) and Ppargc1b (PGC1B) expression patterns
along the CytoTRACE2 pseudotime/potency axis.

Features:
- MAGIC imputation for denoising sparse expression
- Correlation analysis for ALL genes with CytoTRACE2 score
- Trajectory visualization for genes of interest

"""

import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse
from scipy.ndimage import gaussian_filter1d

# MAGIC for imputation
try:
    import magic
    MAGIC_AVAILABLE = True
except ImportError:
    print("WARNING: MAGIC not installed. Install with: pip install magic-impute")
    MAGIC_AVAILABLE = False

warnings.filterwarnings('ignore')
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

# Input AnnData object
ADATA_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad"

# Output directory
OUTPUT_DIR = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/pgc1_cytotrace2_trajectory"

# Genes of interest (mouse gene names)
GENES_OF_INTEREST = ["Ppargc1a", "Ppargc1b"]  # PGC1A and PGC1B in mouse

# CytoTRACE2 columns
CYTOTRACE_SCORE = "cytotrace2_score"     # Lower = more differentiated
CYTOTRACE_POTENCY = "cytotrace2_potency"  # Categorical potency level

# Leiden cluster column
LEIDEN_KEY = "leiden_1.0"

# Condition column
CONDITION_KEY = "condition"

# Analysis parameters
N_BINS = 20  # Number of bins for pseudotime
LOWESS_FRAC = 0.3  # Fraction of data for LOWESS smoothing
USE_MAGIC_IMPUTATION = True  # Use MAGIC for imputation
MAGIC_KNN = 10  # k for MAGIC kNN graph
MAGIC_T = 3  # Diffusion time for MAGIC

# ══════════════════════════════════════════════════════════════════════════════
# 0. Setup
# ══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "figures"), exist_ok=True)
print(f"Output directory: {OUTPUT_DIR}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load AnnData
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nLoading AnnData from: {ADATA_PATH}")
adata = sc.read_h5ad(ADATA_PATH)
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

# Check available columns
print(f"\nAvailable obs columns: {list(adata.obs.columns)}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Check for genes of interest
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nSearching for genes: {GENES_OF_INTEREST}")

# Try exact match first, then case-insensitive
genes_found = {}
for gene in GENES_OF_INTEREST:
    if gene in adata.var_names:
        genes_found[gene] = gene
    else:
        # Try case-insensitive match
        matches = [g for g in adata.var_names if g.lower() == gene.lower()]
        if matches:
            genes_found[gene] = matches[0]
        else:
            # Try partial match
            matches = [g for g in adata.var_names if gene.lower() in g.lower()]
            if matches:
                print(f"  Partial matches for {gene}: {matches[:5]}")
                genes_found[gene] = matches[0]

print(f"Genes found: {genes_found}")

if len(genes_found) == 0:
    raise ValueError(f"None of the genes {GENES_OF_INTEREST} found in the dataset!")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Extract expression and pseudotime data
# ══════════════════════════════════════════════════════════════════════════════

print("\nExtracting expression data...")

# Get expression matrix
if issparse(adata.X):
    X = adata.X.toarray()
else:
    X = np.array(adata.X)

# Check if data needs normalization
data_max = X.max()
print(f"Data max value: {data_max:.2f}")

if data_max > 100:
    print("Normalizing data...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    if issparse(adata_norm.X):
        X = adata_norm.X.toarray()
    else:
        X = np.array(adata_norm.X)
else:
    adata_norm = adata.copy()
    if issparse(adata_norm.X):
        X = adata_norm.X.toarray()
    else:
        X = np.array(adata_norm.X)

# ══════════════════════════════════════════════════════════════════════════════
# 3b. MAGIC Imputation for denoising
# ══════════════════════════════════════════════════════════════════════════════

X_imputed = None
if USE_MAGIC_IMPUTATION and MAGIC_AVAILABLE:
    print("\nRunning MAGIC imputation...")
    print(f"  Parameters: knn={MAGIC_KNN}, t={MAGIC_T}")
    
    try:
        # Create MAGIC operator
        magic_op = magic.MAGIC(knn=MAGIC_KNN, t=MAGIC_T, verbose=False)
        
        # Run MAGIC on the expression matrix
        X_imputed = magic_op.fit_transform(X)
        
        print(f"  MAGIC imputation complete!")
        print(f"  Original data sparsity: {(X == 0).mean():.2%}")
        print(f"  Imputed data sparsity: {(X_imputed == 0).mean():.2%}")
        
    except Exception as e:
        print(f"  WARNING: MAGIC imputation failed: {e}")
        print("  Continuing with original data...")
        X_imputed = None

elif USE_MAGIC_IMPUTATION and not MAGIC_AVAILABLE:
    print("\nWARNING: MAGIC requested but not installed.")
    print("  Install with: pip install magic-impute")
    print("  Continuing with original data...")

# Use imputed data if available, otherwise use normalized data
X_final = X_imputed if X_imputed is not None else X

# Build DataFrame with expression and metadata
df = pd.DataFrame(index=adata.obs_names)

# Add CytoTRACE2 score
if CYTOTRACE_SCORE in adata.obs.columns:
    df['cytotrace2_score'] = adata.obs[CYTOTRACE_SCORE].values
    print(f"\nCytoTRACE2 score range: {df['cytotrace2_score'].min():.3f} - {df['cytotrace2_score'].max():.3f}")
else:
    raise ValueError(f"Column '{CYTOTRACE_SCORE}' not found in adata.obs")

# Add CytoTRACE2 potency if available
if CYTOTRACE_POTENCY in adata.obs.columns:
    df['cytotrace2_potency'] = adata.obs[CYTOTRACE_POTENCY].values

# Add gene expression for genes of interest
print("\nGenes of interest expression (after imputation):")
for gene_name, gene_id in genes_found.items():
    gene_idx = adata_norm.var_names.get_loc(gene_id)
    df[gene_name] = X_final[:, gene_idx]
    print(f"  {gene_name}: mean={df[gene_name].mean():.3f}, max={df[gene_name].max():.3f}")

# Add cluster and condition info
if LEIDEN_KEY in adata.obs.columns:
    df['cluster'] = adata.obs[LEIDEN_KEY].astype(str).values
elif 'leiden' in adata.obs.columns:
    df['cluster'] = adata.obs['leiden'].astype(str).values

if CONDITION_KEY in adata.obs.columns:
    df['condition'] = adata.obs[CONDITION_KEY].astype(str).values

# Create "differentiation pseudotime" (inverted CytoTRACE2 score)
# CytoTRACE2: high score = stem-like, low score = differentiated
# For trajectory: we want 0 = stem-like, 1 = differentiated
df['pseudotime'] = 1 - df['cytotrace2_score']

print(f"\nData prepared: {len(df)} cells")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Save expression data
# ══════════════════════════════════════════════════════════════════════════════

csv_path = os.path.join(OUTPUT_DIR, "pgc1_cytotrace2_expression.csv")
df.to_csv(csv_path)
print(f"\nExpression data saved to: {csv_path}")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Correlation analysis
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("CORRELATION ANALYSIS")
print("="*70)

corr_results = []
for gene_name in genes_found.keys():
    # Spearman correlation with CytoTRACE2 score
    rho_score, p_score = stats.spearmanr(df[gene_name], df['cytotrace2_score'])
    
    # Spearman correlation with pseudotime (differentiation)
    rho_pt, p_pt = stats.spearmanr(df[gene_name], df['pseudotime'])
    
    print(f"\n{gene_name}:")
    print(f"  vs CytoTRACE2 score (stemness): rho={rho_score:.4f}, p={p_score:.2e}")
    print(f"  vs Pseudotime (differentiation): rho={rho_pt:.4f}, p={p_pt:.2e}")
    
    corr_results.append({
        'gene': gene_name,
        'rho_cytotrace2_score': rho_score,
        'pval_cytotrace2_score': p_score,
        'rho_pseudotime': rho_pt,
        'pval_pseudotime': p_pt
    })

corr_df = pd.DataFrame(corr_results)
corr_df.to_csv(os.path.join(OUTPUT_DIR, "pgc1_correlation_results.csv"), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# 5b. ALL GENES Correlation Analysis with CytoTRACE2 Score
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("ALL GENES CORRELATION ANALYSIS")
print("="*70)

print(f"\nCalculating Spearman correlation for all {adata_norm.n_vars} genes...")
print("This may take a few minutes...")

cytotrace_scores = df['cytotrace2_score'].values
all_gene_correlations = []

# Process in batches for progress reporting
batch_size = 1000
n_genes = adata_norm.n_vars

for batch_start in range(0, n_genes, batch_size):
    batch_end = min(batch_start + batch_size, n_genes)
    
    for j in range(batch_start, batch_end):
        gene_name = adata_norm.var_names[j]
        gene_expr = X_final[:, j]
        
        # Calculate Spearman correlation
        rho, pval = stats.spearmanr(gene_expr, cytotrace_scores)
        
        # Handle NaN values
        if not np.isfinite(rho):
            rho = 0.0
        if not np.isfinite(pval):
            pval = 1.0
        
        # Calculate mean expression and percent expressed
        mean_expr = gene_expr.mean()
        pct_expressed = (gene_expr > 0).mean() * 100
        
        all_gene_correlations.append({
            'gene': gene_name,
            'spearman_rho': rho,
            'pvalue': pval,
            'mean_expression': mean_expr,
            'pct_cells_expressed': pct_expressed
        })
    
    # Progress report
    pct_done = (batch_end / n_genes) * 100
    print(f"  Processed {batch_end}/{n_genes} genes ({pct_done:.1f}%)")

# Create DataFrame and calculate FDR-adjusted p-values
all_corr_df = pd.DataFrame(all_gene_correlations)

# BH FDR correction
from scipy.stats import rankdata
pvals = all_corr_df['pvalue'].values
n = len(pvals)
ranks = rankdata(pvals)
fdr = np.minimum(1, pvals * n / ranks)
# Ensure monotonicity
fdr_sorted_idx = np.argsort(pvals)
fdr_sorted = fdr[fdr_sorted_idx]
fdr_monotonic = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
fdr[fdr_sorted_idx] = fdr_monotonic
all_corr_df['fdr_adjusted_pvalue'] = fdr

# Add direction column
all_corr_df['direction'] = np.where(
    all_corr_df['spearman_rho'] > 0, 
    'positive_with_stemness',
    'negative_with_stemness'
)

# Sort by absolute correlation
all_corr_df['abs_rho'] = np.abs(all_corr_df['spearman_rho'])
all_corr_df = all_corr_df.sort_values('abs_rho', ascending=False)

# Save full results (stemness-oriented)
all_corr_path = os.path.join(OUTPUT_DIR, "all_genes_cytotrace2_correlations.csv")
all_corr_df.to_csv(all_corr_path, index=False)
print(f"\nAll genes correlation (stemness) saved to: {all_corr_path}")

# ══════════════════════════════════════════════════════════════════════════════
# 5c. Create DIFFERENTIATION-oriented correlation table
# ══════════════════════════════════════════════════════════════════════════════

print("\nCreating differentiation-oriented correlation table...")

# Differentiation pseudotime = 1 - CytoTRACE2 score
# So correlation with differentiation = -1 * correlation with CytoTRACE2 score
diff_corr_df = all_corr_df.copy()
diff_corr_df['spearman_rho_differentiation'] = -diff_corr_df['spearman_rho']
diff_corr_df['direction_differentiation'] = np.where(
    diff_corr_df['spearman_rho_differentiation'] > 0,
    'increases_with_differentiation',
    'decreases_with_differentiation'
)

# Rename columns for clarity
diff_corr_df = diff_corr_df.rename(columns={
    'spearman_rho': 'spearman_rho_stemness',
    'direction': 'direction_stemness'
})

# Sort by correlation with differentiation (highest first = genes that increase with differentiation)
diff_corr_df['abs_rho_diff'] = np.abs(diff_corr_df['spearman_rho_differentiation'])
diff_corr_df = diff_corr_df.sort_values('spearman_rho_differentiation', ascending=False)

# Reorder columns for clarity
column_order = [
    'gene',
    'spearman_rho_differentiation',
    'spearman_rho_stemness', 
    'pvalue',
    'fdr_adjusted_pvalue',
    'direction_differentiation',
    'direction_stemness',
    'mean_expression',
    'pct_cells_expressed',
    'abs_rho_diff'
]
diff_corr_df = diff_corr_df[column_order]

# Save differentiation-oriented results (all genes)
diff_corr_path = os.path.join(OUTPUT_DIR, "all_genes_differentiation_correlations.csv")
diff_corr_df.to_csv(diff_corr_path, index=False)
print(f"All genes correlation (differentiation) saved to: {diff_corr_path}")

# Split into POSITIVE correlations (genes increasing with differentiation)
positive_diff_df = diff_corr_df[diff_corr_df['spearman_rho_differentiation'] > 0].copy()
positive_diff_df = positive_diff_df.sort_values('spearman_rho_differentiation', ascending=False)
positive_diff_path = os.path.join(OUTPUT_DIR, "positive_differentiation_correlations.csv")
positive_diff_df.to_csv(positive_diff_path, index=False)
print(f"Positive differentiation correlations ({len(positive_diff_df)} genes) saved to: {positive_diff_path}")

# Split into NEGATIVE correlations (genes decreasing with differentiation = stem markers)
negative_diff_df = diff_corr_df[diff_corr_df['spearman_rho_differentiation'] < 0].copy()
negative_diff_df = negative_diff_df.sort_values('spearman_rho_differentiation', ascending=True)  # Most negative first
negative_diff_path = os.path.join(OUTPUT_DIR, "negative_differentiation_correlations.csv")
negative_diff_df.to_csv(negative_diff_path, index=False)
print(f"Negative differentiation correlations ({len(negative_diff_df)} genes) saved to: {negative_diff_path}")

# Also save top 50 for quick reference
top_diff_increase = positive_diff_df.head(50)
top_diff_increase.to_csv(os.path.join(OUTPUT_DIR, "top50_genes_increasing_with_differentiation.csv"), index=False)

top_diff_decrease = negative_diff_df.head(50)
top_diff_decrease.to_csv(os.path.join(OUTPUT_DIR, "top50_genes_decreasing_with_differentiation.csv"), index=False)

# Significant differentiation genes
sig_diff_genes = diff_corr_df[diff_corr_df['fdr_adjusted_pvalue'] < 0.05]
sig_diff_increase = sig_diff_genes[sig_diff_genes['spearman_rho_differentiation'] > 0]
sig_diff_decrease = sig_diff_genes[sig_diff_genes['spearman_rho_differentiation'] < 0]

print(f"\nDifferentiation correlation summary:")
print(f"  Significant genes (FDR < 0.05): {len(sig_diff_genes)}")
print(f"    - Increase with differentiation: {len(sig_diff_increase)}")
print(f"    - Decrease with differentiation: {len(sig_diff_decrease)}")

print(f"\nTop 20 genes INCREASING with differentiation:")
for _, row in top_diff_increase.head(20).iterrows():
    print(f"  {row['gene']}: rho={row['spearman_rho_differentiation']:.4f}, FDR={row['fdr_adjusted_pvalue']:.2e}")

print(f"\nTop 20 genes DECREASING with differentiation (stem markers):")
for _, row in top_diff_decrease.head(20).iterrows():
    print(f"  {row['gene']}: rho={row['spearman_rho_differentiation']:.4f}, FDR={row['fdr_adjusted_pvalue']:.2e}")

# Summary statistics
sig_genes = all_corr_df[all_corr_df['fdr_adjusted_pvalue'] < 0.05]
pos_sig = sig_genes[sig_genes['spearman_rho'] > 0]
neg_sig = sig_genes[sig_genes['spearman_rho'] < 0]

print(f"\nSummary:")
print(f"  Total genes analyzed: {len(all_corr_df)}")
print(f"  Significant genes (FDR < 0.05): {len(sig_genes)}")
print(f"    - Positively correlated with stemness: {len(pos_sig)}")
print(f"    - Negatively correlated with stemness: {len(neg_sig)}")

# Top 20 positive and negative correlations
print(f"\nTop 20 genes POSITIVELY correlated with CytoTRACE2 (stemness):")
top_pos = all_corr_df[all_corr_df['spearman_rho'] > 0].head(20)
for _, row in top_pos.iterrows():
    print(f"  {row['gene']}: rho={row['spearman_rho']:.4f}, FDR={row['fdr_adjusted_pvalue']:.2e}")

print(f"\nTop 20 genes NEGATIVELY correlated with CytoTRACE2 (differentiation markers):")
top_neg = all_corr_df[all_corr_df['spearman_rho'] < 0].head(20)
for _, row in top_neg.iterrows():
    print(f"  {row['gene']}: rho={row['spearman_rho']:.4f}, FDR={row['fdr_adjusted_pvalue']:.2e}")

# Save top correlations separately
top_pos.to_csv(os.path.join(OUTPUT_DIR, "top_positive_stemness_genes.csv"), index=False)
top_neg.to_csv(os.path.join(OUTPUT_DIR, "top_negative_stemness_genes.csv"), index=False)

# Check where PGC1A and PGC1B rank
print(f"\nPGC1A/B rankings in correlation list:")
for gene_name, gene_id in genes_found.items():
    gene_row = all_corr_df[all_corr_df['gene'] == gene_id]
    if len(gene_row) > 0:
        rank = all_corr_df.index.get_loc(gene_row.index[0]) + 1
        rho = gene_row['spearman_rho'].values[0]
        fdr = gene_row['fdr_adjusted_pvalue'].values[0]
        print(f"  {gene_name}: rank={rank}/{len(all_corr_df)}, rho={rho:.4f}, FDR={fdr:.2e}")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Visualization: Scatter plots with trend lines
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating visualizations...")

fig_dir = os.path.join(OUTPUT_DIR, "figures")

# --- Plot 1: Individual scatter plots for each gene ---
for gene_name in genes_found.keys():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: vs CytoTRACE2 score
    ax = axes[0]
    scatter = ax.scatter(
        df['cytotrace2_score'], 
        df[gene_name],
        c=df['pseudotime'],
        cmap='viridis',
        alpha=0.3,
        s=5,
        rasterized=True
    )
    
    # Add trend line (binned means)
    bins = np.linspace(df['cytotrace2_score'].min(), df['cytotrace2_score'].max(), N_BINS + 1)
    df['score_bin'] = pd.cut(df['cytotrace2_score'], bins=bins, labels=False)
    bin_means = df.groupby('score_bin')[gene_name].mean()
    bin_centers = [(bins[i] + bins[i+1])/2 for i in range(N_BINS)]
    
    # Smooth the trend - use bin indices that exist in bin_means
    valid_bins = bin_means.dropna()
    if len(valid_bins) > 3:
        x_valid = [bin_centers[int(b)] for b in valid_bins.index]
        y_smooth = gaussian_filter1d(valid_bins.values, sigma=1.5)
        ax.plot(x_valid, y_smooth, 'r-', linewidth=3, label='Smoothed trend')
    
    ax.set_xlabel('CytoTRACE2 Score (→ more stem-like)', fontsize=12)
    ax.set_ylabel(f'{gene_name} Expression', fontsize=12)
    ax.set_title(f'{gene_name} vs CytoTRACE2 Score', fontsize=14)
    plt.colorbar(scatter, ax=ax, label='Pseudotime')
    ax.legend()
    
    # Right: vs Pseudotime (differentiation)
    ax = axes[1]
    scatter = ax.scatter(
        df['pseudotime'], 
        df[gene_name],
        c=df['cytotrace2_score'],
        cmap='viridis_r',
        alpha=0.3,
        s=5,
        rasterized=True
    )
    
    # Add trend line
    bins = np.linspace(0, 1, N_BINS + 1)
    df['pt_bin'] = pd.cut(df['pseudotime'], bins=bins, labels=False)
    bin_means = df.groupby('pt_bin')[gene_name].mean()
    bin_centers = [(bins[i] + bins[i+1])/2 for i in range(N_BINS)]
    
    valid_bins = bin_means.dropna()
    if len(valid_bins) > 3:
        x_valid = [bin_centers[int(b)] for b in valid_bins.index]
        y_smooth = gaussian_filter1d(valid_bins.values, sigma=1.5)
        ax.plot(x_valid, y_smooth, 'r-', linewidth=3, label='Smoothed trend')
    
    ax.set_xlabel('Differentiation Pseudotime (→ more differentiated)', fontsize=12)
    ax.set_ylabel(f'{gene_name} Expression', fontsize=12)
    ax.set_title(f'{gene_name} vs Differentiation Pseudotime', fontsize=14)
    plt.colorbar(scatter, ax=ax, label='CytoTRACE2 Score')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f"{gene_name}_trajectory.pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(fig_dir, f"{gene_name}_trajectory.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {gene_name}_trajectory.pdf/png")

# --- Plot 2: Combined trajectory plot (both genes) ---
if len(genes_found) >= 2:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = ['#E64B35', '#4DBBD5']  # Red for PGC1A, Blue for PGC1B
    
    bins = np.linspace(0, 1, N_BINS + 1)
    df['pt_bin'] = pd.cut(df['pseudotime'], bins=bins, labels=False)
    bin_centers = [(bins[i] + bins[i+1])/2 for i in range(N_BINS)]
    
    for idx, gene_name in enumerate(genes_found.keys()):
        # Scatter (light)
        ax.scatter(
            df['pseudotime'], 
            df[gene_name],
            c=colors[idx],
            alpha=0.1,
            s=3,
            rasterized=True,
            label=f'{gene_name} (cells)'
        )
        
        # Trend line
        bin_means = df.groupby('pt_bin')[gene_name].mean()
        bin_sems = df.groupby('pt_bin')[gene_name].sem()
        
        valid_bins = bin_means.dropna()
        if len(valid_bins) > 3:
            x_valid = [bin_centers[int(b)] for b in valid_bins.index]
            y_valid = valid_bins.values
            sem_valid = bin_sems.loc[valid_bins.index].values
            
            y_smooth = gaussian_filter1d(y_valid, sigma=1.5)
            ax.plot(x_valid, y_smooth, '-', color=colors[idx], linewidth=3, 
                    label=f'{gene_name} (trend)')
            ax.fill_between(x_valid, y_smooth - sem_valid, y_smooth + sem_valid,
                           color=colors[idx], alpha=0.2)
    
    ax.set_xlabel('Differentiation Pseudotime (CytoTRACE2)', fontsize=14)
    ax.set_ylabel('Expression (log-normalized)', fontsize=14)
    ax.set_title('PGC1A & PGC1B Expression along Differentiation', fontsize=16)
    ax.legend(loc='best', fontsize=10)
    ax.set_xlim(0, 1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "pgc1a_pgc1b_combined_trajectory.pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(fig_dir, "pgc1a_pgc1b_combined_trajectory.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: pgc1a_pgc1b_combined_trajectory.pdf/png")

# --- Plot 3: Heatmap across pseudotime bins ---
if len(genes_found) >= 1:
    fig, ax = plt.subplots(figsize=(12, 3))
    
    bins = np.linspace(0, 1, N_BINS + 1)
    df['pt_bin'] = pd.cut(df['pseudotime'], bins=bins, labels=False)
    
    heatmap_data = []
    for gene_name in genes_found.keys():
        bin_means = df.groupby('pt_bin')[gene_name].mean()
        # Reindex to ensure all bins are present (fill missing with NaN then interpolate)
        bin_means = bin_means.reindex(range(N_BINS))
        bin_means = bin_means.interpolate(method='linear', limit_direction='both')
        # Z-score normalize
        z_scores = (bin_means - bin_means.mean()) / (bin_means.std() + 1e-8)
        heatmap_data.append(z_scores.values)
    
    heatmap_array = np.array(heatmap_data)
    
    im = ax.imshow(heatmap_array, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_yticks(range(len(genes_found)))
    ax.set_yticklabels(list(genes_found.keys()), fontsize=12)
    ax.set_xlabel('Differentiation Pseudotime Bins', fontsize=12)
    ax.set_title('PGC1A/B Expression (z-scored) along Differentiation', fontsize=14)
    
    # Add bin labels
    ax.set_xticks([0, N_BINS//2, N_BINS-1])
    ax.set_xticklabels(['Stem-like', 'Intermediate', 'Differentiated'])
    
    plt.colorbar(im, ax=ax, label='Z-score', shrink=0.8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "pgc1_heatmap_trajectory.pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(fig_dir, "pgc1_heatmap_trajectory.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: pgc1_heatmap_trajectory.pdf/png")

# --- Plot 4: Expression by cluster along pseudotime ---
if 'cluster' in df.columns:
    for gene_name in genes_found.keys():
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Sort clusters by mean pseudotime
        cluster_order = df.groupby('cluster')['pseudotime'].mean().sort_values().index.tolist()
        
        palette = sns.color_palette("husl", n_colors=len(cluster_order))
        
        for idx, cluster in enumerate(cluster_order):
            cluster_df = df[df['cluster'] == cluster]
            ax.scatter(
                cluster_df['pseudotime'],
                cluster_df[gene_name],
                c=[palette[idx]],
                alpha=0.3,
                s=10,
                label=f'Cluster {cluster}',
                rasterized=True
            )
        
        ax.set_xlabel('Differentiation Pseudotime', fontsize=14)
        ax.set_ylabel(f'{gene_name} Expression', fontsize=14)
        ax.set_title(f'{gene_name} Expression by Cluster along Pseudotime', fontsize=16)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_cluster.pdf"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_cluster.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {gene_name}_by_cluster.pdf/png")

# --- Plot 5: Violin plot by potency category ---
if 'cytotrace2_potency' in df.columns:
    for gene_name in genes_found.keys():
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Order potency categories
        potency_order = ['Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent']
        potency_order = [p for p in potency_order if p in df['cytotrace2_potency'].unique()]
        
        if len(potency_order) == 0:
            potency_order = df['cytotrace2_potency'].unique().tolist()
        
        sns.violinplot(
            data=df,
            x='cytotrace2_potency',
            y=gene_name,
            order=potency_order,
            palette='viridis',
            ax=ax
        )
        
        ax.set_xlabel('CytoTRACE2 Potency', fontsize=14)
        ax.set_ylabel(f'{gene_name} Expression', fontsize=14)
        ax.set_title(f'{gene_name} Expression by Potency Category', fontsize=16)
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_potency.pdf"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_potency.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {gene_name}_by_potency.pdf/png")

# --- Plot 6: By condition if available ---
if 'condition' in df.columns and df['condition'].nunique() > 1:
    for gene_name in genes_found.keys():
        fig, ax = plt.subplots(figsize=(12, 6))
        
        conditions = df['condition'].unique()
        colors = sns.color_palette("Set2", n_colors=len(conditions))
        
        bins = np.linspace(0, 1, N_BINS + 1)
        bin_centers = [(bins[i] + bins[i+1])/2 for i in range(N_BINS)]
        
        for idx, cond in enumerate(conditions):
            cond_df = df[df['condition'] == cond].copy()
            cond_df['pt_bin'] = pd.cut(cond_df['pseudotime'], bins=bins, labels=False)
            
            bin_means = cond_df.groupby('pt_bin')[gene_name].mean()
            bin_sems = cond_df.groupby('pt_bin')[gene_name].sem()
            
            # Get valid bins that have data
            valid_bins = bin_means.dropna()
            if len(valid_bins) > 3:
                x_valid = [bin_centers[int(b)] for b in valid_bins.index]
                y_valid = valid_bins.values
                sem_valid = bin_sems.loc[valid_bins.index].fillna(0).values
                
                y_smooth = gaussian_filter1d(y_valid, sigma=1.5)
                ax.plot(x_valid, y_smooth, '-', color=colors[idx], linewidth=2.5, label=cond)
                ax.fill_between(x_valid, 
                               np.array(y_smooth) - np.array(sem_valid), 
                               np.array(y_smooth) + np.array(sem_valid),
                               color=colors[idx], alpha=0.2)
        
        ax.set_xlabel('Differentiation Pseudotime', fontsize=14)
        ax.set_ylabel(f'{gene_name} Expression', fontsize=14)
        ax.set_title(f'{gene_name} Trajectory by Condition', fontsize=16)
        ax.legend(loc='best', fontsize=10)
        ax.set_xlim(0, 1)
        
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_condition.pdf"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(fig_dir, f"{gene_name}_by_condition.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {gene_name}_by_condition.pdf/png")

# --- Plot 7: Volcano plot for all genes correlation ---
print("\nGenerating all-genes correlation volcano plot...")

fig, ax = plt.subplots(figsize=(12, 8))

# Calculate -log10(FDR)
all_corr_df['neg_log10_fdr'] = -np.log10(all_corr_df['fdr_adjusted_pvalue'].clip(lower=1e-300))

# Color by significance
colors = []
for _, row in all_corr_df.iterrows():
    if row['fdr_adjusted_pvalue'] < 0.05:
        if row['spearman_rho'] > 0:
            colors.append('#E64B35')  # Red for positive (stemness)
        else:
            colors.append('#4DBBD5')  # Blue for negative (differentiation)
    else:
        colors.append('#CCCCCC')  # Grey for non-significant

ax.scatter(
    all_corr_df['spearman_rho'],
    all_corr_df['neg_log10_fdr'],
    c=colors,
    alpha=0.5,
    s=10,
    rasterized=True
)

# Highlight PGC1A/B
for gene_name, gene_id in genes_found.items():
    gene_row = all_corr_df[all_corr_df['gene'] == gene_id]
    if len(gene_row) > 0:
        ax.scatter(
            gene_row['spearman_rho'].values[0],
            gene_row['neg_log10_fdr'].values[0],
            c='gold',
            s=150,
            marker='*',
            edgecolors='black',
            linewidths=1,
            zorder=10,
            label=gene_name
        )
        ax.annotate(
            gene_name,
            (gene_row['spearman_rho'].values[0], gene_row['neg_log10_fdr'].values[0]),
            xytext=(10, 10),
            textcoords='offset points',
            fontsize=12,
            fontweight='bold'
        )

ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.7, label='FDR = 0.05')
ax.axvline(0, color='grey', linestyle='-', alpha=0.5)

ax.set_xlabel('Spearman Correlation (ρ) with CytoTRACE2 Score', fontsize=14)
ax.set_ylabel('-log₁₀(FDR)', fontsize=14)
ax.set_title('Gene Correlation with Stemness (CytoTRACE2 Score)', fontsize=16)

# Add annotations for directions
ax.text(0.7, 0.95, '← Differentiation | Stemness →', transform=ax.transAxes, 
        fontsize=10, ha='center', color='grey')

ax.legend(loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "all_genes_correlation_volcano.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, "all_genes_correlation_volcano.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: all_genes_correlation_volcano.pdf/png")

# --- Plot 8: Correlation distribution histogram ---
fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(all_corr_df['spearman_rho'], bins=100, color='steelblue', alpha=0.7, edgecolor='white')

# Mark PGC1A/B positions
for gene_name, gene_id in genes_found.items():
    gene_row = all_corr_df[all_corr_df['gene'] == gene_id]
    if len(gene_row) > 0:
        rho = gene_row['spearman_rho'].values[0]
        ax.axvline(rho, color='red', linestyle='--', linewidth=2, label=f'{gene_name} (ρ={rho:.3f})')

ax.axvline(0, color='black', linestyle='-', alpha=0.5)
ax.set_xlabel('Spearman Correlation (ρ) with CytoTRACE2 Score', fontsize=14)
ax.set_ylabel('Number of Genes', fontsize=14)
ax.set_title('Distribution of Gene Correlations with Stemness', fontsize=16)
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "correlation_distribution.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_dir, "correlation_distribution.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: correlation_distribution.pdf/png")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Summary
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"Cells analyzed: {len(df)}")
print(f"Total genes analyzed: {len(all_corr_df)}")
print(f"Genes of interest: {list(genes_found.keys())}")
print(f"MAGIC imputation: {'Applied' if X_imputed is not None else 'Not applied'}")

print(f"\nPGC1A/B Correlation with CytoTRACE2 score (stemness):")
for _, row in corr_df.iterrows():
    direction = "↑ with stemness" if row['rho_cytotrace2_score'] > 0 else "↓ with stemness"
    print(f"  {row['gene']}: rho={row['rho_cytotrace2_score']:.4f} ({direction})")

print(f"\nAll genes significant correlations (FDR < 0.05): {len(sig_genes)}")
print(f"  Stemness perspective:")
print(f"    - Increase with stemness: {len(pos_sig)}")
print(f"    - Decrease with stemness: {len(neg_sig)}")
print(f"  Differentiation perspective:")
print(f"    - Increase with differentiation: {len(sig_diff_increase)}")
print(f"    - Decrease with differentiation: {len(sig_diff_decrease)}")

print(f"\nOutput saved to: {OUTPUT_DIR}")
print("\nFiles generated:")
print(f"  - pgc1_cytotrace2_expression.csv (PGC1 expression data)")
print(f"  - pgc1_correlation_results.csv (PGC1 statistics)")
print(f"  Stemness-oriented:")
print(f"    - all_genes_cytotrace2_correlations.csv (ALL genes vs stemness)")
print(f"    - top_positive_stemness_genes.csv (top stemness markers)")
print(f"    - top_negative_stemness_genes.csv (top differentiation markers)")
print(f"  Differentiation-oriented:")
print(f"    - all_genes_differentiation_correlations.csv (ALL genes)")
print(f"    - positive_differentiation_correlations.csv (ALL positive, {len(positive_diff_df)} genes)")
print(f"    - negative_differentiation_correlations.csv (ALL negative, {len(negative_diff_df)} genes)")
print(f"    - top50_genes_increasing_with_differentiation.csv")
print(f"    - top50_genes_decreasing_with_differentiation.csv")
print(f"  - figures/ (all visualizations)")

print("\n✓ Analysis complete!")



