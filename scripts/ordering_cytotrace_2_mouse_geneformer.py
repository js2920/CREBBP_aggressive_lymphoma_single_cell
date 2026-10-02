#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CytoTRACE2-ordered gene programs with CLUSTERED y-axis (shape + amplitude)
Size-capped heatmaps; MAGIC smoothing; Viridis colormap + Okabe-Ito palette;
Condition-aware selection; per-program Enrichr + GSEA(prerank) with NES vs p-value plots;
Top annotation bar (dominant category per CT2 bin).

Adapted for: Mouse Geneformer predictions h5ad

I/O:
- Input  H5AD: /home/gusti/CREBBP_aggressive_lymphoma_single_cell/Geneformer/mouse_with_geneformer_predictions.h5ad
- Output dir : /home/gusti/CREBBP_aggressive_lymphoma_single_cell/heatmap_programs

Author: J
Updated: 2025-12-02
"""

# =========================== HARD-CODED SETTINGS ==============================

from pathlib import Path
import sys
import warnings

# Non-interactive backend (avoid X11 on servers)
import matplotlib
matplotlib.use('Agg')

# ---- Project paths ----
INPUT_H5AD   = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/Geneformer/mouse_with_geneformer_predictions.h5ad")
OUTPUT_DIR   = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/heatmap_programs")

# Data layer to use
LAYER_NAME  = "counts"   # raw counts layer

# MAGIC imputation (visualization only)
# Key parameters:
#   t     = diffusion time (higher = smoother; 3-7 typical, auto=~5)
#   knn   = neighbors for graph (higher = smoother global trends)
#   decay = affinity decay (1=standard, 15=more local)
#   n_pca = PCA dims for initial graph
USE_MAGIC   = True
MAGIC_N_PCA = 100       # 100 is good for 10k+ cells
MAGIC_KNN   = 30        # High knn for smooth visualization
MAGIC_DECAY = 1         # Standard decay
MAGIC_T     = 12        # Higher t for stronger smoothing
N_THREADS   = 20

# Condition discrimination
CONDITION_COLUMN = "condition"
USE_CONDITION_DISCRIMINATION = True
CONDITION_INTERACTION_WEIGHT = 0.5  # Balanced

# ---------------- Heatmap size control & aesthetics ---------------------------
HEATMAP_GENE_MIN = 200
HEATMAP_GENE_MAX = 300            # Top 300 discriminatory genes
FIG_HEIGHT_PER_GENE = 0.09        # Larger rows for 300 genes
GENE_LABEL_FONTSIZE = 4           # Larger font for readability

# ---------------- GSEA Filtering (disabled - run on all genes) ----------------
GSEA_MIN_ABS_CORRELATION = 0.0    # No filtering - use all genes
GSEA_MIN_GENES_PER_PROGRAM = 5    # Very low minimum

# Upstream candidate pool (keep larger than cap so you have room to filter)
# N_GENES: Initial pool of genes to consider (before final filtering)
#   - Larger (2000-5000): More comprehensive, slower
#   - Smaller (500-1000): Faster, may miss subtle programs
N_GENES = 3000                    # larger initial pool for better program detection

# FILTER_METHOD: How to select top genes for heatmap
#   - "variance": Genes with highest variance (captures dynamic genes)
#   - "max_expr": Genes with highest peak expression
#   - "mean_expr": Genes with highest mean expression
#   - "combined": Variance × max expression (balanced)
FILTER_METHOD = "combined"        # balanced selection for diverse programs

# Binning / orientation
# CT2_BINS: Number of bins along pseudotime
#   - Set to None or -1 to use INDIVIDUAL CELLS (no binning)
#   - Or set to a very large number (e.g., n_cells)
CT2_BINS         = -1              # Use individual cells (no binning)
ORDER_ASCENDING  = True            # low CT2 (differentiated) → high CT2 (stem-like)

# Visualization Smoothing (Rolling Average)
# When plotting individual cells, data can be noisy.
# Use a rolling window to smooth the visualization (but keep cell ordering).
VISUALIZATION_SMOOTHING_WINDOW = 100  # Increased for smoother dichotomy visualization

# ================== Y-AXIS (PROGRAMS) =========================================
Y_AXIS_STRATEGY       = "cluster"   # "cluster" (RECOMMENDED) or "trajectory"

# --- CLUSTERING OPTIMIZATION PARAMETERS ---
# Settings for 2-program dichotomy (Early vs Late)

# ALPHA_SHAPE: Balance between temporal pattern shape vs amplitude
ALPHA_SHAPE           = 0.5         # Balanced

# SMOOTH_BIN_SIGMA: Gaussian smoothing
SMOOTH_BIN_SIGMA      = 3.0         # Smooth for clear dichotomy

# CLUSTER_METHOD: Hierarchical linkage method
CLUSTER_METHOD        = "ward"      # Ward = compact clusters

DENDRO_THRESHOLD      = 0.8         # High threshold

# AUTO K SELECTION: Disabled - forcing 2 programs
AUTO_CHOOSE_K         = False      
AUTO_K_MIN, AUTO_K_MAX= 2, 2
MIN_GENES_PER_PROGRAM = 100         # ~250 genes per program

# FORCE_N_PROGRAMS: Force exactly 2 programs (Early vs Late)
FORCE_N_PROGRAMS      = 2

# Program stacking & separators
PROGRAM_SORT_BY         = "center_of_mass"  # "center_of_mass" | "peak_bin" | "dendrogram"
DRAW_PROGRAM_SEPARATORS = True

# Expression retention thresholds
# These control which genes pass filtering - adjust based on your data quality
NORM_TARGET_SUM          = 1e4      # counts per 10k normalization

# Stage 1: Global gene filtering (balanced for pathway detection)
MIN_FRAC_CELLS_EXPRESSED = 0.03     # 3% detection rate - captures variable genes like OXPHOS
MIN_MEAN_NORM_EXPR       = 0.05     # 0.05 CP10k - keeps mitochondrial genes

# Stage 2: Trajectory-specific filtering
BIN_MIN_MEAN_NORM        = 0.05     # Minimum expression in bins
MIN_BINS_WITH_SIGNAL     = 5        # Gene must be active in at least 5 bins

# Label coloring by correlation sign + heatmap filtering
POS_CORR_THR, NEG_CORR_THR = 0.15, -0.15  # Lowered to capture more genes
FILTER_HEATMAP_BY_CORRELATION = True  # Only show genes with |r| >= threshold

# Enrichment (ORA via Enrichr)
RUN_ENRICHR     = True
SPECIES         = "mouse"  # or "human"
TOP_ENR_TERMS   = 15

# Global + program-centric GSEA (prerank)
RUN_GSEA_PRERANK_GLOBAL   = True
RUN_GSEA_PRERANK_PER_PROG = True
GSEA_PERMUTATIONS         = 200
GSEA_TOP_LABELS           = 30

# Program scores (optional)
ADD_PROGRAM_SCORES = True

# Figure aesthetics - VIRIDIS COLOR SCHEME
FIG_CMAP       = "viridis"        # Viridis colormap
TITLE_FONTSIZE = 24
AXIS_FONTSIZE  = 14
TICK_FONTSIZE  = 12
MAX_FIG_HEIGHT = 100.0

# --------- Top annotation bar (dominant category per CT2 bin) -----------------
USE_TOP_ANNOTATION_BAR = True
TOP_BAR_COLUMN         = "condition"   # use condition for mouse data
TOP_BAR_EMPTY_COLOR    = (0.9, 0.9, 0.9, 1.0)   # for NA bins
TOP_BAR_REL_HEIGHT     = 0.06

# =============================================================================

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from scipy import sparse
from scipy.ndimage import gaussian_filter1d
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, cophenet
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from statsmodels.stats.multitest import multipletests
except Exception:
    multipletests = None

warnings.filterwarnings("ignore", category=FutureWarning, module="anndata.*")


# ============================== UTILITIES =====================================

def _get_X(adata):
    X = adata.layers[LAYER_NAME] if (LAYER_NAME in adata.layers) else adata.X
    return X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)

def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def _row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)

def _ct2_order_and_bins(ct2: np.ndarray, n_bins: int, ascending: bool = True):
    valid = ~np.isnan(ct2)
    s = ct2[valid]
    
    # Order cells
    order = np.argsort(s)
    order = order if ascending else order[::-1]
    
    # If n_bins is negative or None, use individual cells (no binning)
    if n_bins is None or n_bins < 0:
        n_cells = len(s)
        bins = np.arange(n_cells)  # Each cell is its own bin
        full_bins = np.full_like(ct2, fill_value=-1, dtype=int)
        full_bins[valid] = -1  # Not used for binning, just ordering
        
        # Re-order the full array indices
        full_order = np.where(valid)[0][order]
        uniq_bins = bins
        return full_order, full_bins, uniq_bins
        
    q = pd.qcut(s, q=n_bins, labels=False, duplicates="drop")
    bins = np.asarray(q).astype(int)
    full_bins = np.full_like(ct2, fill_value=-1, dtype=int)
    full_bins[valid] = bins
    full_order = np.where(valid)[0][order]
    uniq_bins = np.unique(bins)
    return full_order, full_bins, uniq_bins

def _mean_norm_per_gene(X, target_sum=NORM_TARGET_SUM):
    n_cells = X.shape[0]
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib = np.maximum(lib, 1e-12)
    inv_lib = (target_sum / lib)
    if sparse.issparse(X):
        v = X.T.dot(inv_lib) / n_cells
        return np.asarray(v).ravel()
    return (X.T @ inv_lib) / n_cells

def _expression_filter_basic(adata):
    X = _get_X(adata)
    n_cells = X.shape[0]
    if sparse.issparse(X):
        det_counts = np.array(X.getnnz(axis=0)).ravel()
    else:
        det_counts = (X.toarray() > 0).sum(axis=0)
    det_frac = det_counts / float(n_cells)
    mean_norm = _mean_norm_per_gene(X, target_sum=NORM_TARGET_SUM)
    keep = (det_frac >= MIN_FRAC_CELLS_EXPRESSED) & (mean_norm >= MIN_MEAN_NORM_EXPR)
    df = pd.DataFrame({
        "gene": adata.var_names.values,
        "detected_frac": det_frac,
        "mean_norm_CP10k": mean_norm,
        "keep_stage1": keep
    }).set_index("gene")
    return keep, df

def _filter_genes_by_importance(adata, genes: list, max_genes: int, method: str = "variance"):
    if len(genes) <= max_genes:
        print(f"  Already have {len(genes)} ≤ {max_genes} genes. No additional filtering needed.")
        return genes
    print(f"  Filtering from {len(genes)} to {max_genes} genes using '{method}' method...")

    gene_idx = adata.var_names.get_indexer(genes)
    ok = gene_idx >= 0
    if not np.all(ok):
        print(f"  Warning: {np.sum(~ok)} requested genes not found in var_names; dropping them.")
        genes = [g for g, k in zip(genes, ok) if k]
        gene_idx = gene_idx[ok]

    Xall = _get_X(adata)
    lib = np.asarray(Xall.sum(axis=1)).ravel()
    lib = np.maximum(lib, 1e-12)
    scale = (NORM_TARGET_SUM / lib)

    X = Xall[:, gene_idx]
    if sparse.issparse(X):
        Xn = (sparse.diags(scale).dot(X)).tocsr()
        Xlog = Xn.copy(); Xlog.data = np.log1p(Xlog.data)
        X_dense = Xlog.toarray()
    else:
        Xn = (X * scale[:, None])
        X_dense = np.log1p(Xn)

    G_var  = np.var(X_dense, axis=0)
    G_max  = np.max(X_dense, axis=0)
    G_mean = np.mean(X_dense, axis=0)

    if method == "variance":   scores = G_var
    elif method == "max_expr": scores = G_max
    elif method == "mean_expr":scores = G_mean
    elif method == "combined":
        v = (G_var - G_var.min()) / (G_var.max() - G_var.min() + 1e-12)
        m = (G_max - G_max.min()) / (G_max.max() - G_max.min() + 1e-12)
        scores = v * m
    else:
        raise ValueError(f"Unknown FILTER_METHOD: {method}")

    top_idx = np.argsort(scores)[-max_genes:]
    selected_genes = [genes[i] for i in top_idx]
    print(f"  Selected {len(selected_genes)} high-importance genes.")
    return selected_genes

def _build_gene_bin_matrices(adata, genes: list, bin_ids: np.ndarray, uniq_bins: np.ndarray, cell_order=None):
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) == 0:
        raise RuntimeError("None of the requested genes are in adata.var_names")

    gene_idx = adata.var_names.get_indexer(genes)
    Xall = _get_X(adata)

    lib = np.asarray(Xall.sum(axis=1)).ravel()
    lib = np.maximum(lib, 1e-12)
    scale = (NORM_TARGET_SUM / lib)

    X = Xall[:, gene_idx]
    if sparse.issparse(X):
        Xn  = (sparse.diags(scale).dot(X)).tocsr()
        Xlog = Xn.copy(); Xlog.data = np.log1p(Xlog.data)
    else:
        Xn  = (X * scale[:, None])
        Xlog = np.log1p(Xn)

    # If using individual cells (bin_ids is -1 or unique bins match cell count)
    if (bin_ids is not None and (bin_ids == -1).all()) or (cell_order is not None):
        if cell_order is None:
            raise ValueError("cell_order must be provided for individual cell mode")
            
        # Reorder cells directly
        Xn_ordered = Xn[cell_order, :]
        Xlog_ordered = Xlog[cell_order, :]
        
        M_norm = Xn_ordered.toarray() if sparse.issparse(Xn_ordered) else Xn_ordered
        M_log = Xlog_ordered.toarray() if sparse.issparse(Xlog_ordered) else Xlog_ordered
        
        # Apply rolling window smoothing for visualization
        if VISUALIZATION_SMOOTHING_WINDOW > 1:
            print(f"  Applying rolling smoothing (window={VISUALIZATION_SMOOTHING_WINDOW}) for visualization...")
            M_norm = pd.DataFrame(M_norm).rolling(window=VISUALIZATION_SMOOTHING_WINDOW, center=True, min_periods=1).mean().values
            M_log = pd.DataFrame(M_log).rolling(window=VISUALIZATION_SMOOTHING_WINDOW, center=True, min_periods=1).mean().values
            
        mat_norm = pd.DataFrame(M_norm.T, index=genes, columns=range(len(cell_order)))
        mat_log  = pd.DataFrame(M_log.T,  index=genes, columns=range(len(cell_order)))
        return mat_norm, mat_log

    def _mean_over_cells(mat, mask):
        if sparse.issparse(mat):
            mat = mat.tocsr()
            idx = np.flatnonzero(mask)
            if idx.size == 0:
                return np.zeros(mat.shape[1], dtype=float)
            return np.asarray(mat[idx, :].mean(axis=0)).ravel()
        else:
            return mat[mask, :].mean(axis=0)

    M_norm, M_log = [], []
    for b in uniq_bins:
        m = (bin_ids == b)
        M_norm.append(_mean_over_cells(Xn,  m))
        M_log.append( _mean_over_cells(Xlog, m))

    M_norm = np.stack(M_norm, axis=1)
    M_log  = np.stack(M_log,  axis=1)
    mat_norm = pd.DataFrame(M_norm, index=genes, columns=range(len(uniq_bins)))
    mat_log  = pd.DataFrame(M_log,  index=genes, columns=range(len(uniq_bins)))
    return mat_norm, mat_log

def _smooth_bins_df(mat: pd.DataFrame, sigma):
    if not sigma or sigma <= 0:
        return mat
    V = gaussian_filter1d(mat.values, sigma=float(sigma), axis=1, mode="nearest")
    return pd.DataFrame(V, index=mat.index, columns=mat.columns)

def _gene_positions(mat: pd.DataFrame):
    V = mat.values
    if (V < 0).any():
        V = V - V.min()
    rs = V.sum(axis=1, keepdims=True); rs[rs == 0] = 1.0
    prob = V / rs
    bins = np.arange(V.shape[1])
    com  = (prob * bins).sum(axis=1)
    peak = V.argmax(axis=1)
    return pd.DataFrame({"peak_bin": peak, "ct2_center_of_mass": com}, index=mat.index)

def _compute_single_gene_corr(i, X, s):
    col = X[:, i]
    v = np.asarray(col.toarray()).ravel() if sparse.issparse(col) else np.asarray(col).ravel()
    v = np.log1p(v)
    if np.std(v) == 0: return i, 0.0, 1.0
    r, p = stats.spearmanr(v, s)
    return i, (0.0 if np.isnan(r) else r), (1.0 if np.isnan(p) else p)

def _compute_gene_correlations(adata, gene_indices=None):
    """Spearman gene–CT2 on CP10k-normalized values."""
    X = _get_X(adata)
    s = adata.obs["ct2_score"].values.astype(float)
    valid = ~np.isnan(s)
    s = s[valid]; X = X[valid, :]

    lib = np.asarray(X.sum(axis=1)).ravel()
    scale = (NORM_TARGET_SUM / np.maximum(lib, 1e-12))
    if sparse.issparse(X):
        X = (sparse.diags(scale).dot(X)).tocsr()
    else:
        X = X * scale[:, None]

    if gene_indices is None:
        Xg = X; var_names = adata.var_names.values; n_genes = adata.n_vars
    else:
        Xg = X[:, gene_indices]; var_names = adata.var_names.values[gene_indices]; n_genes = len(gene_indices)

    corrs = np.zeros(n_genes, dtype=float)
    pvals = np.ones(n_genes, dtype=float)
    print(f"  Computing correlations for {n_genes} genes using {N_THREADS} threads...")
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = {executor.submit(_compute_single_gene_corr, i, Xg, s): i for i in range(n_genes)}
        completed = 0
        for future in as_completed(futures):
            i, r, p = future.result()
            corrs[i] = r; pvals[i] = p
            completed += 1
            if completed % 1000 == 0:
                print(f"  [corr] {completed}/{n_genes}")
    print(f"  [corr] {n_genes}/{n_genes} - Done!")

    df = pd.DataFrame({"gene": var_names, "correlation": corrs, "pvalue": pvals})
    df["abs_correlation"] = np.abs(df["correlation"])
    if multipletests is not None:
        try:
            df["pvalue_adjusted"] = multipletests(df["pvalue"], method="fdr_bh")[1]
            df["significant"] = df["pvalue_adjusted"] < 0.05
        except Exception:
            df["pvalue_adjusted"] = np.nan; df["significant"] = False
    else:
        df["pvalue_adjusted"] = np.nan; df["significant"] = False
    return df

def _label_colors_from_corr(corr_map: pd.Series) -> dict:
    out = {}
    for g, r in corr_map.items():
        if pd.isna(r): out[g] = "black"
        elif r >= POS_CORR_THR: out[g] = "darkred"      # Changed for blues theme
        elif r <= NEG_CORR_THR: out[g] = "darkblue"     # Changed for blues theme
        else: out[g] = "black"
    return out


# ===================== MAGIC IMPUTATION =======================================

def _apply_magic_imputation(adata, layer_name=LAYER_NAME):
    try:
        import magic
    except ImportError:
        print("WARNING: magic-impute not installed. Run: pip install magic-impute")
        print("Skipping MAGIC imputation.")
        return None

    print(f"  Applying MAGIC (t={MAGIC_T}, knn={MAGIC_KNN}, n_pca={MAGIC_N_PCA}, decay={MAGIC_DECAY})...")
    X = _get_X(adata)
    ad_tmp = sc.AnnData(X=X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    sc.pp.normalize_total(ad_tmp, target_sum=NORM_TARGET_SUM)
    X_norm = ad_tmp.X.toarray() if sparse.issparse(ad_tmp.X) else ad_tmp.X

    magic_op = magic.MAGIC(n_pca=MAGIC_N_PCA, knn=MAGIC_KNN, decay=MAGIC_DECAY,
                           t=MAGIC_T, n_jobs=N_THREADS, random_state=42, verbose=1)
    X_magic = magic_op.fit_transform(X_norm)
    
    # Report auto-determined t if applicable
    actual_t = magic_op.t if hasattr(magic_op, 't') else MAGIC_T
    print(f"  MAGIC complete. Actual t={actual_t}")

    adata_magic = sc.AnnData(X=sparse.csr_matrix(X_magic * NORM_TARGET_SUM),
                             obs=adata.obs.copy(), var=adata.var.copy())
    return adata_magic


# ================ CONDITION-AWARE GENE SELECTION ==============================

def _compute_condition_interaction_scores(adata, genes, condition_col=CONDITION_COLUMN):
    if condition_col not in adata.obs.columns:
        print(f"WARNING: '{condition_col}' not found in adata.obs. Skipping condition discrimination.")
        return None

    print(f"  Computing condition interaction scores for {len(genes)} genes...")
    gene_idx = adata.var_names.get_indexer(genes)
    ok = gene_idx >= 0
    if not np.all(ok):
        print(f"  Warning: {np.sum(~ok)} genes not found for interaction scoring; dropping them.")
        genes = [g for g in genes if g in adata.var_names]
        gene_idx = adata.var_names.get_indexer(genes)

    X = _get_X(adata)
    ct2 = adata.obs["ct2_score"].values.astype(float)
    conditions = adata.obs[condition_col].values

    valid = ~np.isnan(ct2)
    X = X[valid, :]; ct2 = ct2[valid]; conditions = conditions[valid]

    # CP10k normalize
    lib = np.asarray(X.sum(axis=1)).ravel()
    scale = (NORM_TARGET_SUM / np.maximum(lib, 1e-12))
    if sparse.issparse(X):
        X = (sparse.diags(scale).dot(X)).tocsr()
    else:
        X = X * scale[:, None]

    unique_conds = np.unique(conditions)
    if len(unique_conds) < 2:
        print(f"  Only {len(unique_conds)} condition found. Need ≥2 for discrimination.")
        return None

    results = []
    for i, g_idx in enumerate(gene_idx):
        if i % 500 == 0:
            print(f"  [interaction] {i}/{len(genes)}")

        col = X[:, g_idx]
        expr = np.asarray(col.toarray()).ravel() if sparse.issparse(col) else np.asarray(col).ravel()
        expr = np.log1p(expr)

        cond_corrs = {}
        for cond in unique_conds:
            mask = (conditions == cond)
            if mask.sum() < 10:
                continue
            expr_cond = expr[mask]; ct2_cond = ct2[mask]
            if np.std(expr_cond) > 0:
                r, _ = stats.spearmanr(expr_cond, ct2_cond)
                cond_corrs[cond] = 0.0 if np.isnan(r) else r

        if len(cond_corrs) < 2:
            interaction_score = 0.0
            overall_corr = 0.0
        else:
            corr_values = list(cond_corrs.values())
            interaction_score = float(np.std(corr_values))
            overall_corr, _ = stats.spearmanr(expr, ct2) if np.std(expr) > 0 else (0.0, 1.0)
            overall_corr = 0.0 if np.isnan(overall_corr) else overall_corr

        results.append({
            "gene": genes[i],
            "interaction_score": interaction_score,
            "overall_abs_corr": abs(overall_corr),
            **{f"corr_{cond}": cond_corrs.get(cond, np.nan) for cond in unique_conds}
        })

    print(f"  [interaction] {len(genes)}/{len(genes)} - Done!")
    return pd.DataFrame(results)

def _select_discriminative_genes(adata, keep_mask, n_genes=N_GENES,
                                 condition_col=CONDITION_COLUMN,
                                 interaction_weight=CONDITION_INTERACTION_WEIGHT):
    idx_stage1 = np.where(keep_mask)[0]
    corr_df = _compute_gene_correlations(adata, gene_indices=idx_stage1)

    if not USE_CONDITION_DISCRIMINATION or (condition_col not in adata.obs.columns):
        print("  Using standard correlation-based selection (no condition discrimination).")
        top = corr_df.sort_values("abs_correlation", ascending=False).head(n_genes)
        return top["gene"].tolist(), corr_df

    pool = max(n_genes * 3, 3000)
    candidate_genes = corr_df.sort_values("abs_correlation", ascending=False).head(pool)["gene"].tolist()
    interact_df = _compute_condition_interaction_scores(adata, candidate_genes, condition_col)
    if interact_df is None:
        print("  Falling back to standard correlation-based selection.")
        top = corr_df.sort_values("abs_correlation", ascending=False).head(n_genes)
        return top["gene"].tolist(), corr_df

    merged = corr_df.merge(interact_df, on="gene", how="inner")
    ac = merged["abs_correlation"]; iscore = merged["interaction_score"]
    merged["norm_abs_corr"]    = (ac - ac.min()) / (ac.max() - ac.min() + 1e-12)
    merged["norm_interaction"] = (iscore - iscore.min()) / (iscore.max() - iscore.min() + 1e-12)
    merged["combined_score"]   = (1 - interaction_weight) * merged["norm_abs_corr"] + \
                                  interaction_weight * merged["norm_interaction"]
    top = merged.sort_values("combined_score", ascending=False).head(n_genes)
    print(f"  Selected {len(top)} genes using combined score "
          f"({(1-interaction_weight):.0%} CT2 corr + {interaction_weight:.0%} interaction)")
    return top["gene"].tolist(), merged


# ===================== Y-AXIS STRATEGIES ======================================

def _combined_distance(mat_log: pd.DataFrame, alpha_shape=ALPHA_SHAPE, smooth_sigma=SMOOTH_BIN_SIGMA):
    M = _smooth_bins_df(mat_log, smooth_sigma)
    Z = _row_zscore(M)
    d_shape = pdist(Z.values, metric="correlation")  # [0,2]
    d_amp   = pdist(M.values, metric="euclidean")
    scale   = np.percentile(d_amp, 95) + 1e-12
    d_amp_s = np.clip(2.0 * (d_amp / scale), 0.0, 2.0)
    return alpha_shape * d_shape + (1.0 - alpha_shape) * d_amp_s

def _auto_select_k(Z, dist, k_min=5, k_max=20, min_genes=30):
    """
    Automatically select optimal number of clusters using multiple metrics.
    Uses silhouette score as primary metric, with Calinski-Harabasz as tiebreaker.
    """
    try:
        from sklearn.metrics import silhouette_score, calinski_harabasz_score
    except Exception as e:
        print(f"[autoK] sklearn not available ({e}); using distance threshold.")
        return None
    
    D = squareform(dist)
    results = []
    
    print(f"[autoK] Evaluating K from {k_min} to {k_max}...")
    for k in range(int(k_min), int(k_max) + 1):
        labels_k = fcluster(Z, t=int(k), criterion="maxclust")
        n_clusters = len(np.unique(labels_k))
        sizes = np.array([(labels_k == c).sum() for c in np.unique(labels_k)])
        
        # Skip if clusters are too small or too few
        if sizes.min() < int(min_genes) or n_clusters < 2:
            continue
        
        try:
            sil = silhouette_score(D, labels_k, metric="precomputed")
            # Calinski-Harabasz needs the actual data, use distance matrix as proxy
            ch = calinski_harabasz_score(D, labels_k)
            
            # Penalize extreme cluster size imbalance
            size_ratio = sizes.max() / sizes.min()
            balance_penalty = 1.0 / (1.0 + 0.1 * (size_ratio - 1))  # mild penalty for imbalance
            
            adjusted_sil = sil * balance_penalty
            results.append({
                'k': k,
                'labels': labels_k,
                'silhouette': sil,
                'adjusted_sil': adjusted_sil,
                'calinski_harabasz': ch,
                'n_clusters': n_clusters,
                'min_size': sizes.min(),
                'max_size': sizes.max(),
                'size_ratio': size_ratio
            })
            print(f"    K={k}: silhouette={sil:.3f}, CH={ch:.1f}, sizes={sizes.min()}-{sizes.max()}")
        except Exception as e:
            print(f"[autoK] k={k} skipped ({e})")
            continue
    
    if not results:
        print("[autoK] no valid k found; falling back to threshold cut.")
        return None
    
    # Sort by adjusted silhouette (accounts for balance)
    results.sort(key=lambda x: x['adjusted_sil'], reverse=True)
    best = results[0]
    
    print(f"[autoK] Selected K={best['k']} (silhouette={best['silhouette']:.3f}, "
          f"n_clusters={best['n_clusters']}, sizes={best['min_size']}-{best['max_size']})")
    
    return best['labels']


def _detect_temporal_patterns(mat_log: pd.DataFrame, n_patterns: int = 6):
    """
    Alternative: Use NMF or PCA to detect temporal expression patterns.
    Returns pattern assignments for genes.
    """
    try:
        from sklearn.decomposition import NMF, PCA
        from sklearn.cluster import KMeans
    except ImportError:
        print("[patterns] sklearn not available, skipping pattern detection")
        return None
    
    # Z-score normalize rows
    Z = _row_zscore(mat_log)
    
    # Use NMF on non-negative shifted data
    Z_shifted = Z - Z.min().min() + 0.01
    
    try:
        nmf = NMF(n_components=n_patterns, init='nndsvd', random_state=42, max_iter=500)
        W = nmf.fit_transform(Z_shifted.values)  # genes × patterns
        H = nmf.components_                       # patterns × bins
        
        # Assign each gene to its dominant pattern
        assignments = W.argmax(axis=1)
        
        # Calculate pattern centers (peak position along pseudotime)
        pattern_centers = H.argmax(axis=1)
        
        return {
            'assignments': assignments,
            'W': W,
            'H': H,
            'pattern_centers': pattern_centers,
            'reconstruction_error': nmf.reconstruction_err_
        }
    except Exception as e:
        print(f"[patterns] NMF failed: {e}")
        return None

def _cluster_and_order_genes(mat_log: pd.DataFrame, mat_norm: pd.DataFrame):
    dist = _combined_distance(mat_log, alpha_shape=ALPHA_SHAPE, smooth_sigma=SMOOTH_BIN_SIGMA)
    Z = linkage(dist, method=CLUSTER_METHOD, optimal_ordering=True)

    if FORCE_N_PROGRAMS is not None:
        labels = fcluster(Z, t=int(FORCE_N_PROGRAMS), criterion="maxclust")
        print(f"[cluster] forcing number of programs: K={int(FORCE_N_PROGRAMS)}")
        # IMPORTANT: keep K as requested — do NOT merge small clusters here.
        # Just warn if any are below MIN_GENES_PER_PROGRAM.
        sizes = pd.Series(labels).value_counts().sort_index()
        small = sizes[sizes < int(MIN_GENES_PER_PROGRAM)]
        if len(small) > 0:
            print(f"[cluster] WARNING: {len(small)} cluster(s) smaller than MIN_GENES_PER_PROGRAM={MIN_GENES_PER_PROGRAM}: {small.to_dict()}")
    elif AUTO_CHOOSE_K:
        labels = _auto_select_k(Z, dist, AUTO_K_MIN, AUTO_K_MAX, MIN_GENES_PER_PROGRAM)
        if labels is None:
            labels = fcluster(Z, t=DENDRO_THRESHOLD, criterion="distance")
    else:
        labels = fcluster(Z, t=DENDRO_THRESHOLD, criterion="distance")

    order = leaves_list(Z)

    leaf_genes  = mat_log.index[order]
    leaf_labels = pd.Series(labels, index=mat_log.index).loc[leaf_genes].values

    seen = {}; prog_ids = []; counter = 0
    for c in leaf_labels:
        if c not in seen:
            counter += 1; seen[c] = counter
        prog_ids.append(seen[c])
    prog_ids = np.array(prog_ids)

    pos_df = _gene_positions(mat_norm.loc[leaf_genes])
    assign = pd.DataFrame({
        "gene": leaf_genes,
        "program_id": prog_ids,
        "leaf_order": np.arange(len(leaf_genes)),
        "peak_bin": pos_df.loc[leaf_genes, "peak_bin"].values,
        "ct2_center_of_mass": pos_df.loc[leaf_genes, "ct2_center_of_mass"].values
    })

    if PROGRAM_SORT_BY == "center_of_mass":
        program_order = assign.groupby("program_id")["ct2_center_of_mass"].median().sort_values().index.tolist()
    elif PROGRAM_SORT_BY == "peak_bin":
        program_order = assign.groupby("program_id")["peak_bin"].median().sort_values().index.tolist()
    else:
        program_order = pd.unique(assign["program_id"]).tolist()

    program_rank = {pid: i for i, pid in enumerate(program_order)}
    assign["program_rank"] = assign["program_id"].map(program_rank)
    assign_sorted = assign.sort_values(["program_rank", "leaf_order"], kind="stable")
    final_genes = assign_sorted["gene"].tolist()

    sizes_by_pid = assign_sorted.groupby("program_id", sort=False).size().tolist()
    boundaries = np.asarray(np.cumsum(sizes_by_pid)[:-1], dtype=float)

    return final_genes, assign_sorted["program_id"].values, Z, labels, dist, boundaries


# ======================= TOP ANNOTATION (SAMPLE ORIGIN) =======================

def _get_colors_from_cmap(name: str, n: int):
    """Return a list of n RGBA tuples sampled from a colormap (handles categorical & continuous)."""
    cmap = matplotlib.colormaps.get_cmap(name)
    if hasattr(cmap, "colors") and isinstance(cmap.colors, (list, tuple)) and len(cmap.colors) >= n:
        return [tuple(cmap.colors[i]) for i in range(n)]
    # otherwise sample evenly
    if n <= 1:
        return [cmap(0.0)]
    xs = np.linspace(0, 1, n, endpoint=False)
    return [tuple(cmap(x)) for x in xs]

def _make_categorical_palette(categories):
    # categories should be a list of distinct labels (not NaN), in display order
    cats = [str(x) for x in categories if pd.notna(x)]
    n = len(cats)
    if n == 0:
        return {}
    
    # Custom color palette for condition progression (user-specified)
    CONDITION_COLORS = {
        "WT_B_cells":        "#98CAE1",  # light blue
        "Crebbp_B_cells":    "#C2E4EF",  # very light cyan
        "Pre_malignant":     "#EAECCC",  # light yellow-green
        "Malignant":         "#FEDA8B",  # light orange
        "Matched_malignant": "#FDB366",  # orange
    }
    
    # Check if all categories are in our custom palette
    cats_set = set(cats)
    if cats_set.issubset(set(CONDITION_COLORS.keys())):
        return {c: mcolors.to_rgba(CONDITION_COLORS[c]) for c in cats}
    
    # Fallback: generate colors for unknown categories
    from matplotlib.cm import viridis
    colors = [viridis(0.2 + 0.6 * i / max(n - 1, 1)) for i in range(n)]
    return {cats[i]: colors[i] for i in range(n)}

def _build_top_bar(adata, bin_ids, uniq_bins, column, ascending=True, cell_order=None):
    if column not in adata.obs.columns:
        return None, None, None
    vals = adata.obs[column]
    categories = sorted(pd.unique(vals[vals.notna()].astype(str)))
    palette = _make_categorical_palette(categories)

    vals_str = vals.astype(str)
    
    # If individual cells (no binning)
    if cell_order is not None:
        bin_labels = vals_str.iloc[cell_order].values
    else:
        bin_labels = []
        for b in uniq_bins:
            mask = (bin_ids == b)
            if mask.sum() == 0:
                bin_labels.append(None)
                continue
            s = vals_str[mask]
            lab = str(s.value_counts().idxmax()) if s.notna().any() else None
            bin_labels.append(lab)
        if not ascending:
            bin_labels = bin_labels[::-1]

    colors = []
    for lab in bin_labels:
        if lab is None or pd.isna(lab):
            colors.append(TOP_BAR_EMPTY_COLOR)
        else:
            colors.append(palette.get(lab, TOP_BAR_EMPTY_COLOR))
    colors = np.asarray(colors)
    return colors, palette, bin_labels


# ================ CONDITION-SPECIFIC VISUALIZATIONS ===========================

def _plot_condition_split_heatmaps(adata, genes, bin_ids, uniq_bins, condition_col,
                                   order, title_prefix, outfile, row_boundaries=None,
                                   ascending=True, label_colors=None):
    if condition_col not in adata.obs.columns:
        print(f"  Condition column '{condition_col}' not found. Skipping condition-split heatmaps.")
        return
    conditions = adata.obs[condition_col].values
    unique_conds = sorted([c for c in np.unique(conditions) if pd.notna(c)])
    if len(unique_conds) < 2:
        print(f"  Only {len(unique_conds)} condition. Skipping condition-split heatmaps.")
        return

    print(f"  Creating condition-split heatmaps for {len(unique_conds)} conditions...")
    n_conds = len(unique_conds)
    ng = len(genes)
    fig_h = max(6.0, min(MAX_FIG_HEIGHT, ng * FIG_HEIGHT_PER_GENE))
    fig, axes = plt.subplots(1, n_conds, figsize=(7*n_conds, fig_h), dpi=300,
                             sharey=True, constrained_layout=True)

    last_quad = None
    for ax, cond in zip(axes, unique_conds):
        cond_mask = (conditions == cond)
        adata_cond = adata[cond_mask].copy()
        bin_ids_cond = bin_ids[cond_mask]
        try:
            _, mat_log_c = _build_gene_bin_matrices(adata_cond, genes, bin_ids_cond, uniq_bins)
            mat_ordered = mat_log_c.iloc[order]
            if not ascending:
                mat_ordered = mat_ordered.iloc[:, ::-1]
            last_quad = ax.pcolormesh(mat_ordered.values, cmap=FIG_CMAP, shading="auto")
            try: last_quad.set_rasterized(True)
            except Exception: pass
            ax.set_title(f"{cond}", fontsize=AXIS_FONTSIZE)
            ax.set_xlabel("CT2 bins", fontsize=AXIS_FONTSIZE-2)
            if (row_boundaries is not None) and (len(row_boundaries) > 0):
                for y in row_boundaries:
                    ax.hlines(y, 0, mat_ordered.shape[1], colors="white", linewidth=0.4, alpha=0.7)
        except Exception as e:
            print(f"  Warning: Could not create heatmap for '{cond}': {e}")
            continue

    ax0 = axes[0]
    ax0.set_yticks(np.arange(ng) + 0.5)
    ax0.set_yticklabels(genes, fontsize=GENE_LABEL_FONTSIZE)
    if label_colors:
        for t in ax0.get_yticklabels():
            col = label_colors.get(t.get_text())
            if col: t.set_color(col)
    ax0.set_ylabel("Genes (clustered order)", fontsize=AXIS_FONTSIZE)

    if last_quad is not None:
        cb = fig.colorbar(last_quad, ax=axes, fraction=0.03, pad=0.02)
        cb.ax.tick_params(labelsize=TICK_FONTSIZE)
        cb.set_label("logₑ(1 + CP10k)", rotation=270, labelpad=15, fontsize=AXIS_FONTSIZE)

    fig.suptitle(title_prefix, fontsize=TITLE_FONTSIZE, y=1.02)
    fig.savefig(outfile, bbox_inches="tight")
    fig.savefig(outfile.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved condition-split heatmap → {outfile}")


# ===================== HEATMAP WITH TOP BAR ===================================

def _plot_heat(mat: pd.DataFrame, order: np.ndarray, title: str, cbar_label: str,
               outfile: Path, vmin=None, vmax=None, label_colors: dict = None,
               x_direction: str = "low→high", row_boundaries=None,
               top_bar_colors: np.ndarray = None, top_bar_label: str = None,
               top_bar_palette: dict = None):
    mat_ordered = mat.iloc[order]
    nb, ng = mat_ordered.shape[1], mat_ordered.shape[0]
    fig_h = max(6.0, min(MAX_FIG_HEIGHT, ng * FIG_HEIGHT_PER_GENE))

    if top_bar_colors is not None:
        # Add space for legend on the right
        fig = plt.figure(figsize=(20, fig_h), dpi=300, constrained_layout=True)
        gs = fig.add_gridspec(nrows=2, ncols=2,
                              height_ratios=[TOP_BAR_REL_HEIGHT, 1.0],
                              width_ratios=[1.0, 0.08])
        ax_top  = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[1, 0], sharex=ax_top)
        ax_legend = fig.add_subplot(gs[0, 1])

        bar_img = top_bar_colors[np.newaxis, :, :]
        ax_top.imshow(bar_img, aspect='auto', interpolation='nearest',
                      extent=[0, nb, 0, 1])
        ax_top.set_yticks([]); ax_top.set_xticks([])
        if top_bar_label:
            ax_top.set_title(f"{top_bar_label}", fontsize=AXIS_FONTSIZE-2, pad=1.5)
        
        # Add legend for conditions
        if top_bar_palette:
            from matplotlib.patches import Patch
            legend_handles = [Patch(facecolor=mcolors.to_hex(c), edgecolor='white', 
                                   label=k.replace("_", " ")) 
                             for k, c in top_bar_palette.items()]
            ax_legend.legend(handles=legend_handles, loc='center left', 
                           frameon=False, fontsize=8, title="Condition", title_fontsize=9)
        ax_legend.axis('off')
        ax = ax_heat
    else:
        fig, ax = plt.subplots(figsize=(18, fig_h), dpi=300, constrained_layout=True)

    quad = ax.pcolormesh(mat_ordered.values, cmap=FIG_CMAP, vmin=vmin, vmax=vmax, shading="auto")
    try: quad.set_rasterized(True)
    except Exception: pass

    n_xticks = min(nb, 11)
    xticks = np.linspace(0, nb-1, n_xticks).astype(int)
    ax.set_xticks(xticks + 0.5)
    ax.set_xticklabels([f"{int(100*t/(nb-1)):d}%" for t in xticks] if nb>1 else ["100%"],
                       rotation=45, ha="right", fontsize=TICK_FONTSIZE)
    ax.set_xlabel(f"Cells ordered by CytoTRACE2 ({nb} bins; {x_direction})", fontsize=AXIS_FONTSIZE)

    ax.set_yticks(np.arange(ng) + 0.5)
    ax.set_yticklabels(mat_ordered.index, fontsize=GENE_LABEL_FONTSIZE)
    if label_colors:
        for t in ax.get_yticklabels():
            col = label_colors.get(t.get_text())
            if col: t.set_color(col)
    ax.set_ylabel("Genes (clustered by pattern & level)", fontsize=AXIS_FONTSIZE)

    if DRAW_PROGRAM_SEPARATORS and (row_boundaries is not None) and len(row_boundaries) > 0:
        for y in np.asarray(row_boundaries, dtype=float):
            ax.hlines(float(y), 0, nb, colors="white", linewidth=0.4, alpha=0.7)

    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=8)
    cb = fig.colorbar(quad, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(cbar_label, rotation=270, labelpad=15, fontsize=AXIS_FONTSIZE)
    cb.ax.tick_params(labelsize=TICK_FONTSIZE)

    fig.savefig(outfile, bbox_inches="tight"); fig.savefig(outfile.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print("Saved →", outfile.resolve())


# ======================= GSEA SCATTER PLOTS ===================================

def _plot_gsea_scatter_for_library(res_df: pd.DataFrame, lib_name: str, out_png: Path,
                                   nes_col: str = "NES", term_col: str = "Term",
                                   p_col_candidates=("FDR q-val", "FDR q-value", "fdr", "FDR", 
                                                     "NOM p-val", "pval", "pvalue", "P-value",
                                                     "Adjusted P-value", "FWER p-val"),
                                   label_top=GSEA_TOP_LABELS):
    sub = res_df[res_df["library"] == lib_name].copy()
    if sub.empty:
        return
    
    # Find p-value column (case-insensitive search)
    p_col = None
    for candidate in p_col_candidates:
        if candidate in sub.columns:
            p_col = candidate
            break
    
    # Fallback: case-insensitive search
    if p_col is None:
        for col in sub.columns:
            if any(p.lower() in col.lower() for p in ['fdr', 'pval', 'p-val', 'pvalue']):
                p_col = col
                break
    
    if p_col is None:
        print(f"[GSEA-plot] No p-value column in results for {lib_name}. Columns: {list(sub.columns)}")
        return

    xvals = sub[p_col].astype(float).clip(lower=1e-300)
    x = -np.log10(xvals)
    y = sub[nes_col].astype(float)
    labels = sub[term_col].astype(str) if term_col in sub.columns else sub.index.astype(str)
    
    # Clean up term names for display
    labels_clean = labels.str.replace("_", " ").str.replace(r"\s*\(GO:\d+\)", "", regex=True)
    labels_clean = labels_clean.str[:55]  # Truncate long names

    # ===== ENHANCED AESTHETICS =====
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    
    # Background styling
    ax.set_facecolor('#FAFAFA')
    fig.patch.set_facecolor('white')
    
    # Point size based on significance
    sizes = 30 + 120 * (x / x.max()).clip(0, 1)  # Larger points for more significant
    
    # Color based on NES (blue = negative, red = positive)
    colors = np.where(y > 0, '#D62728', '#1F77B4')  # Red for positive NES, blue for negative
    
    # Add subtle grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.axhline(0.0, color="black", lw=1, alpha=0.7, zorder=1)
    ax.axvline(-np.log10(0.05), color="gray", lw=1, ls=':', alpha=0.5, label="FDR=0.05")
    
    # Main scatter plot
    scatter = ax.scatter(x, y, s=sizes, alpha=0.7, c=colors, edgecolors='white', linewidths=0.5, zorder=2)
    
    # Labels with adjustText if available, otherwise manual positioning
    try:
        from adjustText import adjust_text
        order = np.argsort(-x.values)
        texts = []
        for idx in order[:int(label_top)]:
            txt = ax.annotate(labels_clean.iloc[idx], (x.iloc[idx], y.iloc[idx]),
                             fontsize=7, fontweight='medium', color='#333333',
                             ha="left", va="center")
            texts.append(txt)
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5),
                   expand_points=(1.5, 1.5), force_text=(0.5, 0.5))
    except ImportError:
        # Fallback: simple annotation
        order = np.argsort(-x.values)
        for idx in order[:int(label_top)]:
            ax.annotate(labels_clean.iloc[idx], (x.iloc[idx] + 0.1, y.iloc[idx]),
                       fontsize=6, fontweight='medium', color='#333333',
                       ha="left", va="center", alpha=0.9)
    
    # Axis labels and title
    ax.set_xlabel(f"−log₁₀({p_col})", fontsize=12, fontweight='medium')
    ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=12, fontweight='medium')
    
    # Clean up library name for title
    lib_display = lib_name.replace("_", " ").replace("2023", "").replace("2022", "").replace("2020", "").strip()
    ax.set_title(f"GSEA: {lib_display}", fontsize=14, fontweight='bold', pad=15)
    
    # Legend for NES direction
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#D62728', markersize=10, label='Positive NES (↑)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1F77B4', markersize=10, label='Negative NES (↓)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', framealpha=0.9, fontsize=9)
    
    # Clean up spines
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#CCCCCC')
    
    ax.tick_params(axis='both', labelsize=10)
    
    plt.tight_layout()
    fig.savefig(out_png, bbox_inches="tight", facecolor='white', edgecolor='none')
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", facecolor='white', edgecolor='none')
    plt.close(fig)
    print("Saved GSEA scatter →", out_png.resolve())

def _plot_custom_diamond_dot(res_df: pd.DataFrame, lib_name: str, out_png: Path, top_n=25):
    """
    Custom plot: Blue Diamonds (NES<0, Differentiated) vs Red Dots (NES>0, Stem-like)
    Both on same axis; point size = number of genes in term.
    """
    sub = res_df[res_df["library"] == lib_name].copy()
    if sub.empty: return

    # Identify p-value column
    p_col = None
    for c in ["FDR q-val", "FDR q-value", "fdr", "Adj P-value", "Adjusted P-value", "NOM p-val"]:
        if c in sub.columns: p_col = c; break
    if not p_col: return

    # Identify gene count column
    gene_col = None
    for c in ["Lead_genes", "lead_genes", "Genes", "genes", "Gene_set_size", "matched size"]:
        if c in sub.columns: 
            gene_col = c
            break
    
    # Data prep
    sub["logp"] = -np.log10(sub[p_col].astype(float).clip(lower=1e-300))
    
    # Get gene counts (if available, otherwise use default)
    if gene_col and gene_col in sub.columns:
        # Count genes (may be semicolon-separated list or integer)
        def count_genes(x):
            if pd.isna(x): return 10
            if isinstance(x, (int, float)): return int(x)
            return len(str(x).split(";"))
        sub["n_genes"] = sub[gene_col].apply(count_genes)
    else:
        sub["n_genes"] = 20  # Default size
    
    # Sort by significance and take top from each direction
    up = sub[sub["NES"] > 0].nlargest(top_n, "logp")
    down = sub[sub["NES"] < 0].nlargest(top_n, "logp")
    plot_data = pd.concat([up, down])
    
    if plot_data.empty: return
    
    # ===== FIGURE =====
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    ax.set_facecolor('#FAFAFA')
    fig.patch.set_facecolor('white')
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5, color='gray')
    
    # Normalize sizes based on gene count
    max_genes = plot_data["n_genes"].max()
    min_size, max_size = 40, 300
    
    # 1. Plot Red Dots (Stem-like / High CT2 / Positive NES)
    stem = plot_data[plot_data["NES"] > 0].copy()
    if not stem.empty:
        s_sizes = min_size + (max_size - min_size) * (stem["n_genes"] / max_genes)
        ax.scatter(stem["logp"], stem["NES"], s=s_sizes, c='#E63946',  # Vivid red
                   marker='o', alpha=0.85, edgecolors='white', linewidth=1.0,
                   label='Stem-like (High CT2)', zorder=3)
               
    # 2. Plot Blue Diamonds (Differentiated / Low CT2 / Negative NES)
    diff = plot_data[plot_data["NES"] < 0].copy()
    if not diff.empty:
        d_sizes = min_size + (max_size - min_size) * (diff["n_genes"] / max_genes)
        ax.scatter(diff["logp"], diff["NES"], s=d_sizes, c='#457B9D',  # Steel blue
                   marker='D', alpha=0.85, edgecolors='white', linewidth=1.0,
                   label='Differentiated (Low CT2)', zorder=3)

    # Reference lines
    ax.axhline(0, color='black', lw=1.2, alpha=0.9, zorder=2)
    ax.axvline(-np.log10(0.05), color='#666666', ls='--', lw=1, alpha=0.6, zorder=1)
    
    # Labels with adjustText
    texts = []
    for _, row in plot_data.iterrows():
        label = str(row["Term"]).replace("_", " ")
        label = label.split("(GO")[0].strip()[:45]  # Remove GO ID, truncate
        txt = ax.text(row["logp"], row["NES"], label, 
                     fontsize=7, fontweight='medium', color='#222222', zorder=4)
        texts.append(txt)
        
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, 
                   arrowprops=dict(arrowstyle='-', color='#888888', alpha=0.6, lw=0.5),
                   expand_points=(1.4, 1.4), force_text=(0.6, 0.6))
    except: pass
    
    # Styling
    lib_clean = lib_name.replace("_", " ").replace("2023", "").replace("2022", "").replace("2020", "").strip()
    ax.set_title(f"{lib_clean}\nGene Set Enrichment by CytoTRACE2 Direction", 
                fontsize=16, fontweight='bold', pad=18)
    ax.set_xlabel(f"Significance (−log₁₀ {p_col})", fontsize=13, fontweight='medium')
    ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=13, fontweight='medium')
    
    # Legend with size indicator
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E63946', 
               markersize=12, label='Stem-like (Positive NES)', markeredgecolor='white'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#457B9D', 
               markersize=10, label='Differentiated (Negative NES)', markeredgecolor='white'),
    ]
    leg = ax.legend(handles=legend_elements, loc='upper left', frameon=True, 
                   framealpha=0.95, fontsize=10, title="Direction", title_fontsize=11)
    leg.get_frame().set_edgecolor('#CCCCCC')
    
    # Add note about size
    ax.text(0.99, 0.01, "Point size ∝ # genes", transform=ax.transAxes, 
           fontsize=8, ha='right', va='bottom', style='italic', color='#666666')
    
    # Clean spines
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#AAAAAA')
        ax.spines[spine].set_linewidth(1)
    
    ax.tick_params(axis='both', labelsize=10)
    
    # Save
    out_custom = out_png.parent / f"custom_diamond_dot_{lib_name}.png"
    plt.tight_layout()
    fig.savefig(out_custom, bbox_inches="tight", facecolor='white', edgecolor='none')
    fig.savefig(out_custom.with_suffix(".svg"), bbox_inches="tight", facecolor='white', edgecolor='none')
    print(f"Saved Custom Diamond/Dot Plot → {out_custom}")
    plt.close(fig)


# ============================== MAIN ==========================================

def _ensure_ct2_score_alias(adata):
    """Create/alias obs['ct2_score'] from common CytoTRACE2 columns."""
    candidates = ["ct2_score", "cytotrace2_score", "cytotrace_score", "score", "CytoTRACE2_score", "CytoTRACE_score"]
    for c in candidates:
        if c in adata.obs.columns:
            adata.obs["ct2_score"] = pd.to_numeric(adata.obs[c], errors="coerce").values
            if c != "ct2_score":
                print(f"[ct2] Using obs['{c}'] as ct2_score")
            return
    raise KeyError("No CT2 score column found (ct2_score / cytotrace2_score / score). "
                   "Please run cytotrace2_run_FIXED_dropin.py first.")

def _gseapy_libs(species: str):
    if species.lower() == "mouse":
        return ['MSigDB_Hallmark_2020','Reactome_2022','KEGG_2019_Mouse',
                'GO_Biological_Process_2023','GO_Molecular_Function_2023','GO_Cellular_Component_2023'], "Mouse"
    else:
        return ['MSigDB_Hallmark_2020','Reactome_2022','KEGG_2021_Human',
                'GO_Biological_Process_2023','GO_Molecular_Function_2023','GO_Cellular_Component_2023'], "Human"

def _run_gsea_prerank_global(selection_df: pd.DataFrame, species: str, outdir: Path, figdir: Path):
    """Global GSEA(prerank): gene ranking by correlation vs CT2 (one run per library)."""
    try:
        import gseapy as gp
    except Exception as e:
        print(f"[GSEA-global] gseapy not installed → skipping. ({e})")
        return

    if selection_df is None or "gene" not in selection_df.columns or "correlation" not in selection_df.columns:
        print("[GSEA-global] No ranking available (need columns: gene, correlation). Skipping.")
        return

    libs, _ = _gseapy_libs(species)
    rnk_df = selection_df[["gene", "correlation"]].dropna().drop_duplicates()
    rnk_df.columns = ["gene", "score"]
    rnk_df = rnk_df.sort_values("score", ascending=False)

    outdir = outdir / "gsea_prerank_global"
    _ensure_dir(outdir)
    _ensure_dir(figdir)

    all_res = []
    for lib in libs:
        print(f"[GSEA-global] prerank vs {lib} ...")
        try:
            pr = gp.prerank(
                rnk=rnk_df,
                gene_sets=lib,
                threads=N_THREADS,
                permutation_num=int(GSEA_PERMUTATIONS),
                outdir=str(outdir / lib.replace(" ", "_")),
                seed=42,
                no_plot=True
            )
            res = None
            if hasattr(pr, "res2d") and isinstance(pr.res2d, pd.DataFrame):
                res = pr.res2d.copy()
            elif hasattr(pr, "results") and isinstance(pr.results, pd.DataFrame):
                res = pr.results.copy()
            if res is None or res.empty:
                print(f"[GSEA-global] No results table for {lib}.")
                continue

            if "NES" not in res.columns and "nes" in res.columns:
                res.rename(columns={"nes": "NES"}, inplace=True)
            if "P-value" not in res.columns and "pval" in res.columns:
                res.rename(columns={"pval": "P-value"}, inplace=True)
            if "FDR q-value" not in res.columns and "fdr" in res.columns:
                res.rename(columns={"fdr": "FDR q-value"}, inplace=True)
            if "Term" not in res.columns:
                res["Term"] = res.index.astype(str)

            res["library"] = lib
            res.to_csv(outdir / f"gsea_prerank_{lib}.csv", index=False)
            all_res.append(res)
        except Exception as e:
            print(f"[GSEA-global] {lib} failed: {e}")

    if not all_res:
        print("[GSEA-global] No libraries produced results.")
        return

    comb = pd.concat(all_res, ignore_index=True)
    comb.to_csv(outdir / "gsea_prerank_global_all.csv", index=False)

    for lib in sorted(comb["library"].unique()):
        _plot_gsea_scatter_for_library(
            res_df=comb, lib_name=lib,
            out_png=figdir / f"gsea_prerank_global_{lib}.png",
            nes_col="NES", term_col="Term"
        )
        # Create the custom Diamond/Dot plot for GO BP
        if "Biological_Process" in lib:
            _plot_custom_diamond_dot(comb, lib, figdir / f"gsea_prerank_global_{lib}.png")
    print(f"[GSEA-global] Saved combined results → {outdir / 'gsea_prerank_global_all.csv'}")

def main():
    print("="*80)
    print("CytoTRACE2-ordered gene programs – MOUSE GENEFORMER DATA")
    print("Blues color scheme")
    print("="*80)
    print(f"Input h5ad  : {INPUT_H5AD}")
    print(f"Output dir  : {OUTPUT_DIR}")
    print(f"Layer       : {LAYER_NAME}")
    print(f"Color scheme: {FIG_CMAP}")
    print(f"MAGIC       : {USE_MAGIC} (t={MAGIC_T}, knn={MAGIC_KNN}, decay={MAGIC_DECAY}) [for visualization]")
    print(f"Condition   : {USE_CONDITION_DISCRIMINATION} (obs['{CONDITION_COLUMN}'])")
    print(f"Heatmap cap : ≤{HEATMAP_GENE_MAX} via '{FILTER_METHOD}'")
    print(f"CT2 bins    : {CT2_BINS} ; ascending={ORDER_ASCENDING}")
    print(f"Force K     : {FORCE_N_PROGRAMS} programs")
    print(f"Row height  : FIG_HEIGHT_PER_GENE={FIG_HEIGHT_PER_GENE}")
    print("="*80)

    if not INPUT_H5AD.exists():
        print(f"ERROR: input not found: {INPUT_H5AD}"); sys.exit(1)
    _ensure_dir(OUTPUT_DIR)
    figdir = OUTPUT_DIR / "figures"; _ensure_dir(figdir)
    progdir = OUTPUT_DIR / "programs"; _ensure_dir(progdir)
    gsea_figdir = figdir / "gsea_scatter"; _ensure_dir(gsea_figdir)

    adata = sc.read_h5ad(INPUT_H5AD)

    # Ensure ct2_score alias exists
    # Your upstream file uses 'cytotrace2_score' — this maps it to 'ct2_score'.
    candidates = ["ct2_score", "cytotrace2_score", "cytotrace_score", "score", "CytoTRACE2_score", "CytoTRACE_score"]
    found = False
    for c in candidates:
        if c in adata.obs.columns:
            adata.obs["ct2_score"] = pd.to_numeric(adata.obs[c], errors="coerce").values
            if c != "ct2_score":
                print(f"[ct2] Using obs['{c}'] as ct2_score")
            found = True
            break
    if not found:
        raise KeyError("No CT2 score column found. Please run cytotrace2_run_FIXED_dropin.py first.")

    if not adata.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    # [0.5] MAGIC (for visualization only)
    adata_for_heatmaps = adata
    if USE_MAGIC:
        print("\n[0.5/10] Applying MAGIC imputation...")
        adata_magic = _apply_magic_imputation(adata)
        if adata_magic is not None:
            adata_for_heatmaps = adata_magic
            adata_magic.write(OUTPUT_DIR / "adata_magic_imputed.h5ad")
            print(f"  Saved MAGIC-imputed data → {OUTPUT_DIR / 'adata_magic_imputed.h5ad'}")
        else:
            print("  Using original data for heatmaps.")

    # [1] Stage-1 filter (on original data)
    print("\n[1/10] Stage-1 gene filter (global detection & mean CP10k)...")
    keep_mask, expr_df = _expression_filter_basic(adata)
    expr_df.to_csv(OUTPUT_DIR / "stage1_expression_metrics.csv")
    print(f"  Kept {keep_mask.sum()}/{adata.n_vars} genes after stage-1 filter.")

    # [2] Discriminative selection (on original normalized data)
    print("\n[2/10] Selecting discriminative genes (CT2 corr + condition interaction) on non-imputed data...")
    selected_genes, selection_df = _select_discriminative_genes(
        adata, keep_mask,
        n_genes=N_GENES, condition_col=CONDITION_COLUMN,
        interaction_weight=CONDITION_INTERACTION_WEIGHT
    )
    selection_df.to_csv(OUTPUT_DIR / "gene_selection_scores.csv", index=False)
    pd.Series(selected_genes).to_csv(OUTPUT_DIR / f"selected_discriminative_genes_{N_GENES}.txt",
                                     index=False, header=False)
    print(f"  Selected {len(selected_genes)} genes.")

    # [3] CT2 ordering & bins
    print(f"\n[3/10] Ordering cells by ct2_score...")
    if CT2_BINS is None or CT2_BINS < 0:
        print(f"  Using INDIVIDUAL CELL ordering (no binning).")
        order_idx, bin_ids, uniq_bins = _ct2_order_and_bins(
            adata.obs["ct2_score"].astype(float).values, -1, ORDER_ASCENDING
        )
        # bin_ids are all -1, uniq_bins are cell indices 0..N-1
        print(f"  Ordered {len(uniq_bins)} cells.")
    else:
        print(f"  Binning into {CT2_BINS} quantile bins...")
        order_idx, bin_ids, uniq_bins = _ct2_order_and_bins(
            adata.obs["ct2_score"].astype(float).values, CT2_BINS, ORDER_ASCENDING
        )
        print(f"  Created {len(uniq_bins)} unique bins.")

    # Save CT2 bin/cell stats
    ct2 = adata.obs["ct2_score"].values
    if CT2_BINS is None or CT2_BINS < 0:
        # Save individual cell stats
        pd.DataFrame({
            "cell_rank": range(len(order_idx)),
            "ct2_score": ct2[order_idx]
        }).to_csv(OUTPUT_DIR / "ct2_cell_stats.csv", index=False)
    else:
        pd.DataFrame([
            (int(b), np.median(ct2[bin_ids == b]) if (bin_ids == b).any() else np.nan)
            for b in uniq_bins
        ], columns=["bin", "ct2_median"]).to_csv(OUTPUT_DIR / "ct2_bin_stats.csv", index=False)

    # [4] Gene × bin/cell matrices (built on heatmap backend)
    print("\n[4/10] Building gene × CT2 mean-expression matrices (CP10k + log1p)...")
    if CT2_BINS is None or CT2_BINS < 0:
        mat_norm, mat_log = _build_gene_bin_matrices(
            adata_for_heatmaps, selected_genes, bin_ids, uniq_bins, cell_order=order_idx
        )
    else:
        mat_norm, mat_log = _build_gene_bin_matrices(adata_for_heatmaps, selected_genes, bin_ids, uniq_bins)
        if not ORDER_ASCENDING:
            mat_norm = mat_norm.iloc[:, ::-1]; mat_log = mat_log.iloc[:, ::-1]

    # [5] Stage-2 filter along trajectory
    print("\n[5/10] Stage-2 gene filter (bin-level signal)...")
    bins_with_signal = (mat_norm >= BIN_MIN_MEAN_NORM).sum(axis=1)
    keep_stage2 = bins_with_signal >= MIN_BINS_WITH_SIGNAL
    mat_norm = mat_norm.loc[keep_stage2]; mat_log = mat_log.loc[keep_stage2]
    pd.DataFrame({"gene": bins_with_signal.index, "bins_with_signal": bins_with_signal.values,
                  "keep_stage2": keep_stage2.values}).to_csv(OUTPUT_DIR / "stage2_bin_signal_metrics.csv", index=False)
    print(f"  Kept {mat_norm.shape[0]}/{len(selected_genes)} genes after stage-2 filter.")

    # Preserve copies before capping
    mat_norm_all = mat_norm.copy()
    mat_log_all  = mat_log.copy()

    # [6] Importance cap to requested size (for plotting)
    print(f"\n[6/10] Enforcing heatmap size ≤{HEATMAP_GENE_MAX} (min target {HEATMAP_GENE_MIN})...")
    if mat_norm.shape[0] > HEATMAP_GENE_MAX:
        final_gene_list = _filter_genes_by_importance(
            adata, mat_norm.index.tolist(), HEATMAP_GENE_MAX, FILTER_METHOD
        )
        mat_norm = mat_norm.loc[final_gene_list]; mat_log = mat_log.loc[final_gene_list]
    else:
        final_gene_list = mat_norm.index.tolist()

    if len(final_gene_list) < HEATMAP_GENE_MIN:
        print(f"  WARNING: Only {len(final_gene_list)} genes remain (< HEATMAP_GENE_MIN={HEATMAP_GENE_MIN}).")
    else:
        print(f"  Using {len(final_gene_list)} genes on the heatmaps (target {HEATMAP_GENE_MIN}–{HEATMAP_GENE_MAX}).")

    # [6b] Filter by correlation threshold (only keep strongly correlated genes)
    if FILTER_HEATMAP_BY_CORRELATION:
        print(f"\n[6b/10] Filtering genes by |correlation| >= {POS_CORR_THR}...")
        corr_lookup = selection_df.set_index("gene")["correlation"]
        keep_corr = []
        for g in final_gene_list:
            r = corr_lookup.get(g, 0)
            if r >= POS_CORR_THR or r <= NEG_CORR_THR:
                keep_corr.append(g)
        n_removed = len(final_gene_list) - len(keep_corr)
        print(f"  Removed {n_removed} weakly correlated genes, keeping {len(keep_corr)}.")
        final_gene_list = keep_corr
        mat_norm = mat_norm.loc[final_gene_list]
        mat_log = mat_log.loc[final_gene_list]

    # Derivative panels
    ref = mat_norm.iloc[:, 0].replace(0, np.nan)
    mat_pct  = (mat_norm.div(ref, axis=0) * 100.0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    eps = 1e-6
    mat_l2fc = np.log2((mat_norm + eps).div(ref + eps, axis=0)).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # [7] Y-axis ordering (programs)
    print(f"\n[7/10] Building y-axis order using strategy: {Y_AXIS_STRATEGY} ...")
    if Y_AXIS_STRATEGY == "cluster":
        final_genes, program_ids, Z, labels, dist, boundaries = _cluster_and_order_genes(mat_log, mat_norm)
        try:
            coph, _ = cophenet(Z, dist); print(f"[cluster] cophenetic_corr={coph:.3f} ; n_clusters={len(np.unique(program_ids))}")
        except Exception as e:
            print(f"[cluster] validation skipped: {e}")
        pos_df = _gene_positions(mat_norm.loc[final_genes])
    else:
        raise NotImplementedError("Set Y_AXIS_STRATEGY='cluster' (default).")

    print(f"\n  Identified {len(np.unique(program_ids))} programs (y-axis).")

    # [8] Assignments & colors
    print("\n[8/10] Saving assignments...")
    if "correlation" in selection_df.columns:
        corr_map = selection_df.set_index("gene")["correlation"].reindex(final_genes)
    else:
        gi = adata.var_names.get_indexer(final_genes)
        corr_small = _compute_gene_correlations(adata, gene_indices=gi)
        corr_map = corr_small.set_index("gene")["correlation"].reindex(final_genes)
    label_colors = _label_colors_from_corr(corr_map)

    assign = pd.DataFrame({
        "gene": final_genes,
        "program_id": program_ids,
        "order_index": np.arange(len(final_genes)),
        "ct2_corr": corr_map.values,
        "peak_bin": pos_df.loc[final_genes, "peak_bin"].values,
        "ct2_center_of_mass": pos_df.loc[final_genes, "ct2_center_of_mass"].values
    })
    progdir = OUTPUT_DIR / "programs"; _ensure_dir(progdir)
    assign.to_csv(progdir / "programs_assignment.csv", index=False)
    for pid, sub in assign.groupby("program_id", sort=False):
        (progdir / f"program_P{int(pid):02d}_genes.txt").write_text("\n".join(sub["gene"]), encoding="utf-8")

    # Top annotation (dominant category per CT2 bin)
    top_bar_info = None
    if USE_TOP_ANNOTATION_BAR:
        col_choice = TOP_BAR_COLUMN if TOP_BAR_COLUMN in adata.obs.columns else (
            CONDITION_COLUMN if CONDITION_COLUMN in adata.obs.columns else None
        )
        if col_choice is None:
            print(f"[top-bar] Column '{TOP_BAR_COLUMN}' not found and no '{CONDITION_COLUMN}' fallback. Skipping top bar.")
        else:
            # Pass cell_order if using individual cells
            cell_ord = order_idx if (CT2_BINS is None or CT2_BINS < 0) else None
            colors_bar, pal, labels_bar = _build_top_bar(
                adata, bin_ids, uniq_bins, col_choice, 
                ascending=ORDER_ASCENDING, cell_order=cell_ord
            )
            if colors_bar is not None:
                top_bar_info = {"colors": colors_bar, "label": col_choice, "palette": pal}
                pal_df = pd.DataFrame({"category": list(pal.keys()),
                                       "color_hex": [mcolors.to_hex(pal[k]) for k in pal.keys()]})
                pal_df.to_csv(figdir / f"top_bar_palette_{col_choice}.tsv", sep="\t", index=False)
                print(f"[top-bar] Saved palette → {figdir / f'top_bar_palette_{col_choice}.tsv'}")

    # [9] Heatmaps (with top annotation bar)
    print("\n[9/10] Plotting heatmaps (clustered y-axis; compact rows; top annotation)...")
    idx_order = np.arange(len(final_genes))
    def _re(mat): return mat.reindex(index=final_genes)
    xdir_text = "low→high" if ORDER_ASCENDING else "high→low"
    vmax_abs = np.nanpercentile(mat_norm.values, 99)

    figdir = OUTPUT_DIR / "figures"; _ensure_dir(figdir)

    _plot_heat(_re(mat_norm), idx_order,
               f"{len(final_genes)} genes – normalized mean (CP10k)",
               "CP10k", figdir / "ct2_heat_abs.png", vmin=0, vmax=vmax_abs,
               label_colors=label_colors, x_direction=xdir_text,
               row_boundaries=boundaries,
               top_bar_colors=(top_bar_info["colors"] if top_bar_info else None),
               top_bar_label=(top_bar_info["label"] if top_bar_info else None),
               top_bar_palette=(top_bar_info["palette"] if top_bar_info else None))

    _plot_heat(_re(mat_log), idx_order,
               f"{len(final_genes)} genes – logₑ(1 + CP10k)",
               "logₑ(1 + CP10k)", figdir / "ct2_heat_log1p.png",
               label_colors=label_colors, x_direction=xdir_text,
               row_boundaries=boundaries,
               top_bar_colors=(top_bar_info["colors"] if top_bar_info else None),
               top_bar_label=(top_bar_info["label"] if top_bar_info else None),
               top_bar_palette=(top_bar_info["palette"] if top_bar_info else None))

    _plot_heat(_re(mat_pct), idx_order,
               f"{len(final_genes)} genes – % of baseline (leftmost bin)",
               "% of baseline", figdir / "ct2_heat_pct_of_leftmost.png", vmin=0, vmax=120,
               label_colors=label_colors, x_direction=xdir_text,
               row_boundaries=boundaries,
               top_bar_colors=(top_bar_info["colors"] if top_bar_info else None),
               top_bar_label=(top_bar_info["label"] if top_bar_info else None),
               top_bar_palette=(top_bar_info["palette"] if top_bar_info else None))

    _plot_heat(_re(mat_l2fc), idx_order,
               f"{len(final_genes)} genes – log₂ FC vs baseline (leftmost bin)",
               "log₂ FC", figdir / "ct2_heat_log2fc_vs_leftmost.png", vmin=-3, vmax=3,
               label_colors=label_colors, x_direction=xdir_text,
               row_boundaries=boundaries,
               top_bar_colors=(top_bar_info["colors"] if top_bar_info else None),
               top_bar_label=(top_bar_info["label"] if top_bar_info else None),
               top_bar_palette=(top_bar_info["palette"] if top_bar_info else None))

    # [10] Condition-specific panels
    if USE_CONDITION_DISCRIMINATION:
        print("\n[10/10] Creating condition-specific heatmaps...")
        _plot_condition_split_heatmaps(
            adata_for_heatmaps, final_genes, bin_ids, uniq_bins,
            CONDITION_COLUMN, idx_order,
            f"Condition-specific programs ({len(final_genes)} genes)",
            figdir / "ct2_heat_condition_split.png",
            row_boundaries=boundaries,
            ascending=ORDER_ASCENDING,
            label_colors=label_colors
        )

    # ---------------- Enrichr per program + Program-centric GSEA ---------------
    if RUN_ENRICHR or RUN_GSEA_PRERANK_PER_PROG:
        try:
            import gseapy as gp
        except Exception as e:
            print(f"[Enrichment] gseapy not installed → skipping enrichment. ({e})")
            gp = None

        if gp is not None:
            libs, org = _gseapy_libs(SPECIES)

            # Program-centric base: z-scored log matrices across bins
            Z_all = _row_zscore(mat_log_all)  # genes × bins
            L = Z_all.shape[1]

            all_enr_rows = []
            for pid, sub in assign.groupby("program_id", sort=False):
                genes_prog = [g for g in sub["gene"] if g in Z_all.index]
                print(f"[Program {int(pid):02d}] n_genes={len(genes_prog)}")

                # ---- Enrichr ORA per program ----
                if RUN_ENRICHR:
                    try:
                        enr = gp.enrichr(gene_list=genes_prog, gene_sets=libs,
                                         organism=org,
                                         outdir=str(progdir / f"gsea_P{int(pid):02d}"),
                                         cutoff=0.05, verbose=False)
                        res = None
                        if hasattr(enr, "results") and isinstance(enr.results, pd.DataFrame):
                            res = enr.results.copy()
                        elif hasattr(enr, "res2d") and isinstance(enr.res2d, pd.DataFrame):
                            res = enr.res2d.copy()
                        if res is not None and not res.empty:
                            cn = {c.lower().strip(): c for c in res.columns}
                            term_col = cn.get("term") or cn.get("name") or "Term"
                            adjp_col = cn.get("adjusted p-value") or cn.get("fdr q-value") or cn.get("adj p") or "Adjusted P-value"
                            p_col   = cn.get("p-value") or "P-value"
                            comb_col= cn.get("combined score") or "Combined Score"
                            gs_col  = cn.get("gene_set") or "Gene_set"
                            if gs_col not in res.columns: res[gs_col] = "unknown"
                            res_small = res.sort_values(adjp_col, ascending=True).head(TOP_ENR_TERMS).copy()
                            res_small["program_id"] = int(pid)
                            res_small.rename(columns={
                                term_col: "Term", adjp_col: "Adjusted P-value",
                                p_col: "P-value", comb_col: "Combined Score",
                                gs_col: "library"
                            }, inplace=True)
                            res_small.to_csv(progdir / f"gsea_enrichr_P{int(pid):02d}.csv", index=False)
                            all_enr_rows.append(res_small)
                        else:
                            print(f"[Enrichr] Program P{int(pid):02d}: no terms retained at cutoff.")
                    except Exception as e:
                        print(f"[Enrichr] Program P{int(pid):02d} failed: {e}")

                # ---- Program-centric GSEA(prerank) ----
                if RUN_GSEA_PRERANK_PER_PROG:
                    try:
                        centroid = Z_all.loc[genes_prog].mean(axis=0).values  # length L
                        centroid = (centroid - centroid.mean()) / (centroid.std() + 1e-12)
                        scores = (Z_all.values @ centroid) / float(L)
                        rnk = pd.Series(scores, index=Z_all.index).sort_values(ascending=False)
                        rnk_df = rnk.reset_index()
                        rnk_df.columns = ["gene", "score"]

                        gsea_combined = []
                        for gs in libs:
                            try:
                                pr = gp.prerank(rnk=rnk_df, gene_sets=gs,
                                                threads=N_THREADS,
                                                permutation_num=int(GSEA_PERMUTATIONS),
                                                outdir=str(progdir / f"gsea_prerank_P{int(pid):02d}_{gs}"),
                                                seed=42, no_plot=True)
                                res = None
                                if hasattr(pr, "res2d") and isinstance(pr.res2d, pd.DataFrame):
                                    res = pr.res2d.copy()
                                elif hasattr(pr, "results") and isinstance(pr.results, pd.DataFrame):
                                    res = pr.results.copy()
                                if res is None or res.empty:
                                    continue
                                if "NES" not in res.columns and "nes" in res.columns:
                                    res.rename(columns={"nes":"NES"}, inplace=True)
                                if "P-value" not in res.columns and "pval" in res.columns:
                                    res.rename(columns={"pval":"P-value"}, inplace=True)
                                if "FDR q-value" not in res.columns and "fdr" in res.columns:
                                    res.rename(columns={"fdr":"FDR q-value"}, inplace=True)
                                if "Term" not in res.columns:
                                    res["Term"] = res.index.astype(str)
                                res["library"] = gs
                                res["program_id"] = int(pid)
                                gsea_combined.append(res)
                            except Exception as e:
                                print(f"[GSEA-prerank] P{int(pid):02d}, lib {gs} failed: {e}")

                        if gsea_combined:
                            gsea_df = pd.concat(gsea_combined, ignore_index=True)
                            gsea_df.to_csv(progdir / f"gsea_prerank_P{int(pid):02d}_combined.csv", index=False)
                            # Plot per library
                            for lib in sorted(gsea_df["library"].unique()):
                                _plot_gsea_scatter_for_library(
                                    res_df=gsea_df, lib_name=lib,
                                    out_png=gsea_figdir / f"P{int(pid):02d}_gsea_{lib}.png",
                                    nes_col="NES", term_col="Term"
                                )
                        else:
                            print(f"[GSEA-prerank] Program P{int(pid):02d}: no results across libraries.")
                    except Exception as e:
                        print(f"[Ranking/Plot] Program P{int(pid):02d} ranking/plot failed: {e}")

            if all_enr_rows:
                df_enr_all = pd.concat(all_enr_rows, ignore_index=True)
                df_enr_all.to_csv(progdir / "gsea_enrichr_top_terms_all_programs.csv", index=False)
                print("Saved Enrichr summary →", (progdir / "gsea_enrichr_top_terms_all_programs.csv").resolve())

    # Global GSEA (prerank) + scatter plots
    if RUN_GSEA_PRERANK_GLOBAL:
        print("\n[+] Running GLOBAL GSEA (prerank) and plotting NES vs p-value (per library) ...")
        _run_gsea_prerank_global(selection_df=selection_df, species=SPECIES,
                                 outdir=progdir, figdir=gsea_figdir)

    # Correlation-subset enrichment (positive vs negative) - GSEA prerank + Enrichr
    if FILTER_HEATMAP_BY_CORRELATION:
        print("\n[+] Running GSEA prerank on correlation subsets (positive vs negative)...")
        try:
            import gseapy as gp
            libs, org = _gseapy_libs(SPECIES)
            
            # Get all genes with correlations (not just heatmap genes)
            corr_df = selection_df[["gene", "correlation"]].dropna().copy()
            
            # Positive correlation genes (stem-like) - use ALL genes with r > 0
            pos_df = corr_df[corr_df["correlation"] > 0].copy()
            pos_df = pos_df.sort_values("correlation", ascending=False)
            pos_df.columns = ["gene", "score"]
            
            # Negative correlation genes (differentiated) - use ALL genes with r < 0, flip sign for ranking
            neg_df = corr_df[corr_df["correlation"] < 0].copy()
            neg_df["correlation"] = neg_df["correlation"].abs()  # Use absolute value for ranking
            neg_df = neg_df.sort_values("correlation", ascending=False)
            neg_df.columns = ["gene", "score"]
            
            print(f"  Positive correlation subset: {len(pos_df)} genes (stem-like)")
            print(f"  Negative correlation subset: {len(neg_df)} genes (differentiated)")
            
            for subset_name, rnk_df in [("positive_corr_stemlike", pos_df), 
                                         ("negative_corr_differentiated", neg_df)]:
                if len(rnk_df) < 15:
                    print(f"  Skipping {subset_name}: too few genes ({len(rnk_df)})")
                    continue
                
                subset_dir = progdir / f"gsea_prerank_{subset_name}"
                _ensure_dir(subset_dir)
                
                all_res = []
                for lib in libs:
                    try:
                        pre_res = gp.prerank(rnk=rnk_df, gene_sets=lib, threads=N_THREADS,
                                            permutation_num=int(GSEA_PERMUTATIONS),
                                            outdir=str(subset_dir / lib.replace(" ", "_")),
                                            seed=42, no_plot=True)
                        res = pre_res.res2d
                        if res is not None and not res.empty:
                            res["library"] = lib
                            res["subset"] = subset_name
                            all_res.append(res)
                    except Exception as e:
                        print(f"    {subset_name} {lib}: {e}")
                
                if all_res:
                    combined = pd.concat(all_res, ignore_index=True)
                    combined.to_csv(subset_dir / f"gsea_prerank_{subset_name}_all.csv", index=False)
                    print(f"  Saved → gsea_prerank_{subset_name}_all.csv ({len(combined)} terms)")
                    
                    # Plot scatter for each library
                    for lib in sorted(combined["library"].unique()):
                        _plot_gsea_scatter_for_library(
                            res_df=combined, lib_name=lib,
                            out_png=gsea_figdir / f"gsea_prerank_{subset_name}_{lib}.png",
                            nes_col="NES", term_col="Term"
                        )
                        
        except Exception as e:
            print(f"  Correlation subset GSEA failed: {e}")

    # Program scores
    if ADD_PROGRAM_SCORES:
        print("\n[+] Computing per-cell program scores (mean log1p CP10k) ...")
        X = _get_X(adata)
        ad_tmp = sc.AnnData(X=X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
        sc.pp.normalize_total(ad_tmp, target_sum=NORM_TARGET_SUM); sc.pp.log1p(ad_tmp)
        for pid, sub in assign.groupby("program_id", sort=False):
            genes_in = [g for g in sub["gene"] if g in ad_tmp.var_names]
            Xi = ad_tmp[:, genes_in].X
            mean_vec = np.asarray(Xi.mean(axis=1)).ravel() if sparse.issparse(Xi) else Xi.mean(axis=1)
            adata.obs[f"program_P{int(pid):02d}_score"] = mean_vec
        out_h5ad = OUTPUT_DIR / "adata_with_program_scores.h5ad"
        adata.write(out_h5ad); print("Saved →", out_h5ad.resolve())

    print("\n" + "="*80)
    print("✓ ANALYSIS COMPLETE")
    print("="*80)
    print(f"Results → {OUTPUT_DIR.resolve()}")
    print("  Figures : ct2_heat_[abs|log1p|pct_of_leftmost|log2fc_vs_leftmost].png/.svg (with top annotation bar)")
    print("            ct2_heat_condition_split.png")
    print("            figures/gsea_scatter/* (NES vs p-value plots)")
    print("  Tables  : programs_assignment.csv + per-program .txt")
    print("            gene_selection_scores.csv, stage1_expression_metrics.csv, stage2_bin_signal_metrics.csv, ct2_bin_stats.csv")
    print("            top_bar_palette_*.tsv")
    if RUN_ENRICHR:
        print("  Enrichr : programs/gsea_enrichr_Pxx.csv, gsea_enrichr_top_terms_all_programs.csv")
    if FILTER_HEATMAP_BY_CORRELATION:
        print("  Subsets : programs/gsea_prerank_positive_corr_stemlike/, gsea_prerank_negative_corr_differentiated/")
    if RUN_GSEA_PRERANK_GLOBAL or RUN_GSEA_PRERANK_PER_PROG:
        print("  GSEA    : programs/gsea_prerank_*_combined.csv and per-library CSVs")
    if USE_MAGIC:          print("  H5AD    : adata_magic_imputed.h5ad")
    if ADD_PROGRAM_SCORES: print("  H5AD    : adata_with_program_scores.h5ad")

if __name__ == "__main__":
    main()



