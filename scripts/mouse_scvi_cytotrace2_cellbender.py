#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mouse scVI + CytoTRACE2 Analysis (CellBender Filtered Data)
============================================================

Loads CellBender-filtered h5 files, performs QC (MT% <= 10%), 
doublet removal via Scrublet, scVI integration, and CytoTRACE2.

Author: J
Date: 2025-12-01
"""

# ============================== SETUP ========================================
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

import re
import shlex
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

import torch
import scvi

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
scvi.settings.seed = 0
np.random.seed(0)
try:
    torch.set_float32_matmul_precision("high")
except Exception:
    pass

# ============================== PATHS ========================================
CELLBENDER_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/cellbender_filtered")
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
FIGDIR_CT2 = OUTDIR / "figures_ct2"
CT2_WORKDIR = OUTDIR / "ct2_io"
CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"
CT2_OUTDIR = OUTDIR / "cytotrace2_results"

for p in (OUTDIR, FIGDIR, FIGDIR_CT2, CT2_WORKDIR, CT2_OUTDIR):
    p.mkdir(parents=True, exist_ok=True)

# CellBender filtered h5 files with sample metadata
SAMPLE_FILES = {
    "SIGAA3_Matched_malignant_R1": {
        "path": CELLBENDER_DIR / "SIGAA3_Matched_malignant_R1_GEX_cellbender_filtered.h5",
        "condition": "Matched_malignant",
        "replicate": "R1"
    },
    "SIGAA4_Matched_malignant_R2": {
        "path": CELLBENDER_DIR / "SIGAA4_Matched_malignant_R2_GEX_cellbender_filtered.h5",
        "condition": "Matched_malignant",
        "replicate": "R2"
    },
    "SIGAA6_Crebbp_B_cells_R1": {
        "path": CELLBENDER_DIR / "SIGAA6_Crebbp_B_cells_R1_GEX_cellbender_filtered.h5",
        "condition": "Crebbp_B_cells",
        "replicate": "R1"
    },
    "SIGAB6_WT_B_cells_R2": {
        "path": CELLBENDER_DIR / "SIGAB6_WT_B_cells_R2_GEX_cellbender_filtered.h5",
        "condition": "WT_B_cells",
        "replicate": "R2"
    },
    "SIGAC2_Pre_malignant_R2": {
        "path": CELLBENDER_DIR / "SIGAC2_Pre_malignant_R2_GEX_cellbender_filtered.h5",
        "condition": "Pre_malignant",
        "replicate": "R2"
    },
    "SIGAC6_WT_B_cells_R1": {
        "path": CELLBENDER_DIR / "SIGAC6_WT_B_cells_R1_GEX_cellbender_filtered.h5",
        "condition": "WT_B_cells",
        "replicate": "R1"
    },
    "SIGAD5_Malignant_R2": {
        "path": CELLBENDER_DIR / "SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R2"
    },
    "SIGAD6_Crebbp_B_cells_R1_2": {
        "path": CELLBENDER_DIR / "SIGAD6_Crebbp_B_cells_R1_GEX_cellbender_filtered.h5",
        "condition": "Crebbp_B_cells",
        "replicate": "R1_2"
    },
    "SIGAF2_Pre_malignant_R1": {
        "path": CELLBENDER_DIR / "SIGAF2_Pre_malignant_R1_GEX_cellbender_filtered.h5",
        "condition": "Pre_malignant",
        "replicate": "R1"
    },
    "SIGAH1_Malignant_R1": {
        "path": CELLBENDER_DIR / "SIGAH1_Malignant_R1_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R1"
    },
}

# QC thresholds
MT_THRESHOLD = 10.0  # Maximum mitochondrial gene percentage
MIN_GENES = 200
MIN_CELLS = 3

print("=" * 84)
print("MOUSE scVI + CytoTRACE2 ANALYSIS (CellBender Filtered)")
print("=" * 84)
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"Output dir    : {OUTDIR}")
print(f"MT threshold  : {MT_THRESHOLD}%\n")

# ============================== GENE SETS ====================================
OXPHOS_GENES = [
    'Cox4i1','Cox5a','Cox5b','Cox6a1','Cox6b1','Cox6c','Cox7a2','Cox7b','Cox7c','Cox8a',
    'Cyc1','Cycs','Ndufa1','Ndufa2','Ndufa3','Ndufa4','Ndufa5','Ndufa6','Ndufa7','Ndufa8','Ndufa9',
    'Ndufa10','Ndufa11','Ndufa12','Ndufa13','Ndufab1','Ndufb1','Ndufb2','Ndufb3','Ndufb4','Ndufb5',
    'Ndufb6','Ndufb7','Ndufb8','Ndufb9','Ndufb10','Ndufb11','Ndufc1','Ndufc2','Ndufs1','Ndufs2',
    'Ndufs3','Ndufs4','Ndufs5','Ndufs6','Ndufs7','Ndufs8','Ndufv1','Ndufv2','Ndufv3','Sdha','Sdhb',
    'Sdhc','Sdhd','Uqcr10','Uqcr11','Uqcrb','Uqcrc1','Uqcrc2','Uqcrfs1','Uqcrh','Uqcrq',
    'Atp5f1a','Atp5f1b','Atp5f1c','Atp5f1d','Atp5f1e','Atp5mc1','Atp5mc2','Atp5mc3','Atp5me','Atp5mf',
    'Atp5mg','Atp5pb','Atp5pd','Atp5pf','Atp5po'
]
BCR_GENES = [
    'Cd79a','Cd79b','Cd19','Cd22','Cd72','Cr2','Fcrl1','Fcrl2','Fcrl3','Fcrl4','Fcrl5','Ms4a1',
    'Ighm','Ighd','Igha','Ighg1','Ighg2a','Ighg2b','Ighg2c','Ighg3','Ighe','Btk','Lyn','Syk','Blk',
    'Blnk','Pik3cd','Pik3ap1','Plcg2','Prkcb','Nfkb1','Nfkb2','Rel','Rela','Nfatc1','Nfatc2',
    'Bcl10','Card11','Malt1','Map3k7','Ikbkb','Ikbkg','Chuk','Ptpn6','Ptprc','Vav1','Vav2',
    'Vav3','Grb2','Sos1','Sos2','Hras','Kras','Nras','Raf1','Map2k1','Map2k2','Mapk1','Mapk3'
]

# ============================== HELPER FUNCTIONS =============================
def _is_intlike(mat, n_check=200000) -> bool:
    """Check if matrix contains integer-like values."""
    if sp.issparse(mat):
        data = mat.data[:min(n_check, mat.data.size)]
    else:
        flat = np.ravel(mat)
        data = flat[:min(n_check, flat.size)]
    return data.size > 0 and np.all((data >= 0) & np.isclose(data, np.round(data)))


def load_cellbender_h5(h5_path: Path, sample_name: str) -> ad.AnnData:
    """Load CellBender filtered h5 file."""
    print(f"  Loading: {h5_path.name}")
    
    # Try reading as 10x h5
    try:
        adata = sc.read_10x_h5(str(h5_path))
    except Exception as e1:
        print(f"    read_10x_h5 failed: {e1}, trying read_h5ad...")
        try:
            adata = sc.read_h5ad(str(h5_path))
        except Exception as e2:
            # Try generic h5 read
            import h5py
            with h5py.File(h5_path, 'r') as f:
                print(f"    H5 keys: {list(f.keys())}")
            raise RuntimeError(f"Could not read {h5_path}: {e1}, {e2}")
    
    # Make var_names unique
    adata.var_names_make_unique()
    
    # Add sample metadata
    adata.obs["sample_id"] = sample_name
    
    # Ensure we have a proper barcode index
    if not adata.obs_names.str.contains("-").any():
        adata.obs_names = [f"{bc}-{sample_name}" for bc in adata.obs_names]
    else:
        adata.obs_names = [f"{bc.split('-')[0]}-{sample_name}" for bc in adata.obs_names]
    
    print(f"    Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


def run_scrublet_safe(adata: ad.AnnData, sample_name: str) -> pd.Series:
    """
    Run Scrublet for doublet detection with proper error handling.
    Returns a boolean Series indicating predicted doublets.
    """
    import scrublet as scr
    
    print(f"    Running Scrublet on {sample_name}...")
    
    # Get count matrix
    if sp.issparse(adata.X):
        counts = adata.X.tocsr()
    else:
        counts = sp.csr_matrix(adata.X)
    
    # Ensure counts are non-negative integers
    counts.data = np.clip(counts.data, 0, None)
    counts.data = np.round(counts.data)
    
    try:
        # Initialize Scrublet
        scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
        
        # Run doublet detection
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            verbose=False
        )
        
        # If automatic threshold fails, use manual threshold
        if predicted_doublets is None or np.all(~predicted_doublets):
            threshold = 0.25
            predicted_doublets = doublet_scores > threshold
            print(f"      Using manual threshold {threshold}")
        
        n_doublets = predicted_doublets.sum()
        print(f"      Detected {n_doublets:,} doublets ({100*n_doublets/len(predicted_doublets):.1f}%)")
        
        return pd.Series(predicted_doublets, index=adata.obs_names)
        
    except Exception as e:
        print(f"      Scrublet failed: {e}")
        print(f"      Marking all cells as singlets for {sample_name}")
        return pd.Series(False, index=adata.obs_names)


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Compute QC metrics including mitochondrial gene percentage."""
    # Identify mitochondrial genes (mouse: mt-)
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt"], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    return adata


def filter_cells(adata: ad.AnnData, mt_threshold: float = 10.0, 
                 min_genes: int = 200) -> ad.AnnData:
    """Filter cells based on QC metrics."""
    n_before = adata.n_obs
    
    # Filter by minimum genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter by MT percentage
    adata = adata[adata.obs["pct_counts_mt"] <= mt_threshold].copy()
    
    n_after = adata.n_obs
    print(f"    Filtered: {n_before:,} → {n_after:,} cells "
          f"({n_before - n_after:,} removed, {100*(n_before-n_after)/n_before:.1f}%)")
    
    return adata


def confounder_mask(varnames: pd.Index) -> pd.Series:
    """Create mask for confounder genes (MT, ribosomal, IG, TCR)."""
    v = pd.Index([str(g) for g in varnames])
    is_mt = v.str.lower().str.startswith("mt-")
    is_ribo = v.str.match(r"^(Rps|Rpl|RPS|RPL)", na=False)
    is_ig = v.str.match(r"^(Igh|Igk|Igl|IGH|IGK|IGL)[vdjc]", case=False, na=False)
    tcr_prefixes = ("Trav", "Trbv", "Trgv", "Trdv", "Traj", "Trbj", "Trgj", "Trdj",
                    "Trac", "Trbc", "Trgc", "Trdc")
    is_tcr = v.str.lower().str.startswith(tuple(p.lower() for p in tcr_prefixes))
    return is_mt | is_ribo | is_ig | is_tcr


# ============================== CT2 FUNCTIONS ================================
def stream_counts_to_ct2_txt(adata: ad.AnnData, dest: Path, 
                              collapse_duplicates: bool = False) -> None:
    """Write CytoTRACE2-compatible matrix: genes × cells."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    
    if "counts" not in adata.layers:
        raise RuntimeError("[ct2-io] layers['counts'] missing for CT2 export.")
    
    cells = adata.obs_names.astype(str).tolist()
    X = adata.layers["counts"]
    is_sparse = sp.issparse(X)
    n_vars = adata.n_vars
    
    na_like = {"", "N/A", "NA", "NULL", "NONE", "NAN", "<NA>", "<na>"}
    
    print(f"[ct2-io] → {dest}  (genes: {n_vars:,}; cells: {len(cells):,})")
    with open(dest, "w", buffering=1024*1024) as fh:
        fh.write("gene\t" + "\t".join(cells) + "\n")
        
        wrote = 0
        for j in range(n_vars):
            g = str(adata.var_names[j]).strip()
            if g.upper() in na_like:
                continue
            col = X[:, j]
            arr = col.toarray().ravel() if is_sparse else np.asarray(col).ravel()
            if np.all(np.isfinite(arr)) and np.all(np.isclose(arr, np.round(arr))):
                fh.write(g + "\t" + "\t".join(map(lambda v: str(int(v)), arr)) + "\n")
            else:
                fh.write(g + "\t" + "\t".join(map(lambda v: f"{float(v):.6g}", arr)) + "\n")
            wrote += 1
            if (wrote % 1000 == 0) or (j + 1 == n_vars):
                print(f"  [ct2-io] {wrote:,} genes written", end="\r")
    print()


def try_import_cytotrace2():
    """Try to import CytoTRACE2 Python API."""
    try:
        from cytotrace2_py.cytotrace2_py import cytotrace2 as fn
        return fn
    except Exception:
        try:
            from cytotrace2_py import cytotrace2 as fn
            return fn
        except Exception:
            return None


def run_ct2_python(input_txt: Path, species: str, outdir: Path):
    """Run CytoTRACE2 via Python API."""
    fn = try_import_cytotrace2()
    if fn is None:
        print("[ct2] Python API not found; trying CLI fallback.")
        return None
    print("[ct2] Running via Python API cytotrace2(...)")
    try:
        return fn(str(input_txt), species=species, output_dir=str(outdir))
    except TypeError:
        return fn(str(input_txt))


def _standardize_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize column names."""
    df = df.copy()
    df.columns = [re.sub(r"\s+", "_", c.strip().lower()) for c in df.columns]
    return df


def _coerce_float_series(s: pd.Series) -> pd.Series:
    """Coerce series to float."""
    if pd.api.types.is_numeric_dtype(s):
        return s.astype(float)
    t = s.astype(str).str.strip().str.replace(",", ".", regex=False)
    t = t.str.replace(r"[^0-9eE\.\+\-]+", "", regex=True)
    return pd.to_numeric(t, errors="coerce")


def parse_ct2_scores(outdir: Path, adata: ad.AnnData) -> pd.DataFrame:
    """Parse CytoTRACE2 output scores."""
    cands: List[Path] = []
    for ext in ("*.csv", "*.tsv", "*.txt"):
        cands += list(outdir.rglob(ext))
    
    best = None
    for f in cands:
        try:
            df = pd.read_csv(f, sep=None, engine="python", dtype=str)
        except Exception:
            continue
        df = _standardize_cols(df)
        if df.empty:
            continue
        
        id_col = next((c for c in ("cell", "cell_id", "barcode", "barcodes", 
                                    "cellname", "cell_id_or_barcode")
                       if c in df.columns), df.columns[0])
        score_cols = [c for c in df.columns if ("cytotrace2" in c and "score" in c)] \
                  or [c for c in df.columns if c in ("score", "cytotrace_score", "ct2_score")]
        
        if not score_cols:
            continue
        
        tmp = pd.DataFrame(index=df[id_col].astype(str).values)
        tmp["cytotrace2_score"] = _coerce_float_series(df[score_cols[0]]).values
        
        pot_cols = [c for c in df.columns if "potency" in c]
        if pot_cols:
            tmp["cytotrace2_potency"] = df[pot_cols[0]].astype(str).values
        
        rel_cols = [c for c in df.columns if "relative" in c]
        if rel_cols:
            tmp["cytotrace2_relative"] = _coerce_float_series(df[rel_cols[0]]).values
        
        n_match = len(set(tmp.index) & set(map(str, adata.obs_names)))
        if best is None or n_match > best[0]:
            best = (n_match, f, tmp)
    
    if best is None:
        raise FileNotFoundError(f"[ct2] No results table found under {outdir}")
    
    print(f"[ct2] Using: {best[1]} (matched {best[0]:,} cells)")
    df_best = best[2]
    df_best = df_best[~df_best.index.duplicated(keep="first")]
    return df_best.reindex(adata.obs_names.astype(str))


# ============================== PLOTTING FUNCTIONS ===========================
def _save_umap(adata, color, fname, title=None, palette=None, cmap=None):
    """Save UMAP plot."""
    fig, ax = plt.subplots(figsize=(10, 9))
    sc.pl.umap(adata, color=color, title=(title or color), palette=palette, cmap=cmap,
               frameon=False, legend_loc="right margin", ax=ax, show=False, s=25)
    plt.tight_layout()
    plt.savefig(FIGDIR / fname, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  ✓ {fname}")


def save_umaps_by_sample(adata: ad.AnnData, figdir: Path, group_key: str = "sample_id"):
    """Save individual UMAPs highlighting each sample (with background in gray)."""
    if group_key not in adata.obs:
        print(f"  (warn) '{group_key}' not in adata.obs; skipping per-group UMAPs.")
        return
    
    labels = adata.obs[group_key].astype(str)
    outdir = figdir / f"umap_by_{group_key}"
    outdir.mkdir(parents=True, exist_ok=True)
    
    background_label = "__rest__"
    
    for grp in sorted(labels.unique()):
        adata.obs["_highlight"] = background_label
        adata.obs.loc[labels == grp, "_highlight"] = grp
        adata.obs["_highlight"] = pd.Categorical(
            adata.obs["_highlight"], categories=[grp, background_label]
        )
        palette = {grp: "#d62728", background_label: "#d3d3d3"}
        
        fig, ax = plt.subplots(figsize=(9, 8))
        sc.pl.umap(
            adata,
            color="_highlight",
            palette=palette,
            frameon=False,
            legend_loc=None,
            show=False,
            s=25,
            ax=ax,
            title=f"{group_key} → {grp}"
        )
        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        adata.obs.drop(columns="_highlight", inplace=True)
    
    print(f"  ✓ Per-group UMAPs saved to {outdir}")


def save_umaps_individual_samples(adata: ad.AnnData, figdir: Path, 
                                   group_key: str = "sample_id",
                                   color_by: Optional[List[str]] = None):
    """
    Save UMAPs showing ONLY cells from each sample (no background cells).
    
    Parameters
    ----------
    adata : AnnData
        Full dataset with UMAP coordinates
    figdir : Path
        Output directory for figures
    group_key : str
        Column in obs to group by (e.g., 'sample_id', 'condition')
    color_by : list of str, optional
        Additional columns to color by for each sample subset.
        Default: ['leiden_0.5', 'condition']
    """
    if group_key not in adata.obs:
        print(f"  (warn) '{group_key}' not in adata.obs; skipping individual sample UMAPs.")
        return
    
    if color_by is None:
        color_by = ['leiden_0.5']
        if 'condition' in adata.obs.columns and group_key != 'condition':
            color_by.append('condition')
        if 'cytotrace2_score' in adata.obs.columns:
            color_by.append('cytotrace2_score')
        if 'phase' in adata.obs.columns:
            color_by.append('phase')
    
    outdir = figdir / f"umap_individual_{group_key}"
    outdir.mkdir(parents=True, exist_ok=True)
    
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    
    # Get global UMAP limits for consistent axes across all plots
    umap_coords = adata.obsm["X_umap"]
    x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
    y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    
    # Color palettes
    condition_colors = {
        "WT_B_cells": "#2ecc71",
        "Crebbp_B_cells": "#3498db", 
        "Pre_malignant": "#f39c12",
        "Matched_malignant": "#e74c3c",
        "Malignant": "#8e44ad"
    }
    
    phase_colors = {
        "G1": "#3498db",
        "S": "#e74c3c", 
        "G2M": "#2ecc71"
    }
    
    for grp in unique_groups:
        # Subset to only this group's cells
        mask = labels == grp
        adata_sub = adata[mask].copy()
        n_cells = adata_sub.n_obs
        
        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)
        
        # Plot 1: Just the cells colored by a single color
        fig, ax = plt.subplots(figsize=(9, 8))
        ax.scatter(
            adata_sub.obsm["X_umap"][:, 0],
            adata_sub.obsm["X_umap"][:, 1],
            c="#d62728",
            s=15,
            alpha=0.7,
            rasterized=True
        )
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_title(f"{grp} (n={n_cells:,})", fontsize=14)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_aspect('equal', adjustable='box')
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        
        # Additional plots colored by different variables
        for col in color_by:
            if col not in adata_sub.obs.columns:
                continue
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            if col == 'cytotrace2_score' or adata_sub.obs[col].dtype in [np.float64, np.float32, float]:
                # Continuous variable
                scatter = ax.scatter(
                    adata_sub.obsm["X_umap"][:, 0],
                    adata_sub.obsm["X_umap"][:, 1],
                    c=adata_sub.obs[col].astype(float),
                    cmap="viridis",
                    s=15,
                    alpha=0.8,
                    rasterized=True
                )
                plt.colorbar(scatter, ax=ax, label=col, shrink=0.8)
            else:
                # Categorical variable
                categories = adata_sub.obs[col].astype(str).unique()
                
                # Use appropriate color palette
                if col == 'condition':
                    palette = condition_colors
                elif col == 'phase':
                    palette = phase_colors
                else:
                    # Generate colors for categories
                    from matplotlib import cm
                    cmap_cat = cm.get_cmap('tab20', len(categories))
                    palette = {cat: cmap_cat(i) for i, cat in enumerate(sorted(categories))}
                
                for cat in sorted(categories):
                    cat_mask = adata_sub.obs[col].astype(str) == cat
                    color = palette.get(cat, "#999999")
                    ax.scatter(
                        adata_sub.obsm["X_umap"][cat_mask, 0],
                        adata_sub.obsm["X_umap"][cat_mask, 1],
                        c=[color],
                        s=15,
                        alpha=0.7,
                        label=cat,
                        rasterized=True
                    )
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            
            ax.set_xlim(x_min - x_margin, x_max + x_margin)
            ax.set_ylim(y_min - y_margin, y_max + y_margin)
            ax.set_title(f"{grp} — {col} (n={n_cells:,})", fontsize=14)
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            ax.set_aspect('equal', adjustable='box')
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
            
            safe_col = re.sub(r"[^A-Za-z0-9._-]+", "_", col)
            fig.savefig(outdir / f"umap_{safe_grp}_{safe_col}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
    
    print(f"  ✓ Individual sample UMAPs saved to {outdir}")


def add_ct2_umap_and_plots(adata: ad.AnnData, figdir: Path, title_suffix=""):
    """Add CytoTRACE2 UMAP and histogram plots."""
    if "X_umap" not in adata.obsm:
        rep = "X_scvi" if "X_scvi" in adata.obsm else None
        if rep is None:
            sc.pp.pca(adata, n_comps=50, use_highly_variable=False)
            rep = "X_pca"
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=30)
        sc.tl.umap(adata, min_dist=0.2, spread=1.5)
    
    figdir.mkdir(parents=True, exist_ok=True)
    
    if "cytotrace2_score" in adata.obs:
        fig, ax = plt.subplots(figsize=(10, 9))
        sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=30, frameon=False,
                   cmap="viridis", title=f"CytoTRACE2 score {title_suffix}")
        fig.savefig(figdir / "umap_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        
        fig, ax = plt.subplots(figsize=(7, 5))
        adata.obs["cytotrace2_score"].astype(float).plot(kind="hist", bins=50, ax=ax)
        ax.set_xlabel("CytoTRACE2 score")
        ax.set_ylabel("Cell count")
        ax.set_title("Distribution of CT2 scores")
        fig.savefig(figdir / "hist_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


# ============================== MAIN PIPELINE ================================
print("=" * 84)
print("STEP 1 — Load and QC CellBender filtered samples")
print("=" * 84)

adata_list = []
qc_stats = []

for sample_name, sample_info in SAMPLE_FILES.items():
    h5_path = sample_info["path"]
    
    if not h5_path.exists():
        print(f"  WARNING: {h5_path} not found, skipping...")
        continue
    
    # Load sample
    adata = load_cellbender_h5(h5_path, sample_name)
    
    # Add metadata
    adata.obs["condition"] = sample_info["condition"]
    adata.obs["replicate"] = sample_info["replicate"]
    
    # Compute QC metrics
    adata = compute_qc_metrics(adata)
    
    # Store pre-filter stats
    n_pre = adata.n_obs
    
    # Run Scrublet for doublet detection
    doublet_mask = run_scrublet_safe(adata, sample_name)
    adata.obs["predicted_doublet"] = doublet_mask.values
    
    # Filter cells (MT threshold = 10%)
    adata = filter_cells(adata, mt_threshold=MT_THRESHOLD, min_genes=MIN_GENES)
    
    # Remove doublets
    n_pre_doublet = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_doublets_removed = n_pre_doublet - adata.n_obs
    print(f"    Removed {n_doublets_removed:,} doublets → {adata.n_obs:,} cells")
    
    # Store counts in layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    
    # Record QC stats
    qc_stats.append({
        "sample": sample_name,
        "condition": sample_info["condition"],
        "cells_raw": n_pre,
        "cells_after_qc": adata.n_obs,
        "pct_removed": 100 * (n_pre - adata.n_obs) / n_pre
    })
    
    adata_list.append(adata)
    print(f"    Final: {adata.n_obs:,} cells\n")

# Save QC stats
qc_df = pd.DataFrame(qc_stats)
qc_df.to_csv(OUTDIR / "qc_statistics.csv", index=False)
print(f"\n  QC stats saved to: {OUTDIR / 'qc_statistics.csv'}")
print(qc_df.to_string())

print("\n" + "=" * 84)
print("STEP 2 — Concatenate samples")
print("=" * 84)

adata = ad.concat(adata_list, join="outer", index_unique=None)
print(f"  Combined: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Ensure counts layer exists after concat
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
print(f"  After gene filter (min_cells={MIN_CELLS}): {adata.n_vars:,} genes")

print("\n" + "=" * 84)
print("STEP 3 — Normalize and find HVGs")
print("=" * 84)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["normalized"] = adata.X.copy()

# Find HVGs (using counts layer)
sc.pp.highly_variable_genes(
    adata, 
    layer="counts",
    n_top_genes=5000, 
    flavor="seurat_v3",
    batch_key="sample_id",
    subset=False
)

# Remove confounder genes from HVGs
hvg_mask = adata.var["highly_variable"].copy()
confounder = confounder_mask(adata.var_names)
hvg_mask = hvg_mask & ~confounder
adata.var["highly_variable"] = hvg_mask
n_hvgs = hvg_mask.sum()
print(f"  HVGs after removing confounders: {n_hvgs:,}")

print("\n" + "=" * 84)
print("STEP 4 — Cell cycle scoring (BEFORE scVI for regression)")
print("=" * 84)

# Cell cycle genes (mouse)
s_genes = ['Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
           'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Mlf1ip', 'Hells', 'Rfc2',
           'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
           'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1',
           'Tipin', 'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b',
           'Brip1', 'E2f8']
g2m_genes = ['Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
             'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
             'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
             'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
             'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
             'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
             'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa']

# Filter to genes present in data
s_genes_present = [g for g in s_genes if g in adata.var_names]
g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]

cell_cycle_computed = False
if s_genes_present and g2m_genes_present:
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
    cell_cycle_computed = True
    print(f"  Cell cycle scoring: {len(s_genes_present)} S genes, {len(g2m_genes_present)} G2M genes")
    print(f"  Phase distribution:")
    print(adata.obs['phase'].value_counts().to_string())
    
    # Calculate cell cycle difference score for regression
    # This is a common approach: regress out the difference between G2M and S scores
    adata.obs['cc_difference'] = adata.obs['G2M_score'] - adata.obs['S_score']
    print(f"  Added 'cc_difference' score for regression (G2M - S)")
else:
    print(f"  WARNING: Not enough cell cycle genes found!")
    print(f"    S genes found: {len(s_genes_present)}")
    print(f"    G2M genes found: {len(g2m_genes_present)}")

print("\n" + "=" * 84)
print("STEP 5 — scVI integration WITH cell cycle regression")
print("=" * 84)

# Subset to HVGs for scVI
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# Setup scVI
for cat in ("sample_id", "condition"):
    if cat in adata_hvg.obs:
        adata_hvg.obs[cat] = adata_hvg.obs[cat].astype("category")

# Prepare continuous covariates for cell cycle regression
continuous_covariates = []
if cell_cycle_computed:
    # Use S_score and G2M_score as continuous covariates for regression
    continuous_covariates = ["S_score", "G2M_score"]
    print(f"  Cell cycle regression enabled using: {continuous_covariates}")
else:
    print(f"  Cell cycle regression DISABLED (no scores available)")

# Setup anndata with cell cycle covariates
if continuous_covariates:
    scvi.model.SCVI.setup_anndata(
        adata_hvg, 
        layer="counts", 
        batch_key="sample_id",
        categorical_covariate_keys=["condition"],
        continuous_covariate_keys=continuous_covariates
    )
else:
    scvi.model.SCVI.setup_anndata(
        adata_hvg, 
        layer="counts", 
        batch_key="sample_id",
        categorical_covariate_keys=["condition"]
    )

model = scvi.model.SCVI(
    adata_hvg,
    n_latent=96,
    n_layers=2,
    dropout_rate=0.1,
    gene_likelihood="nb",
    dispersion="gene-batch",
    use_layer_norm="both",
    use_batch_norm="none"
)

max_epochs = 125
print(f"  Training scVI for up to {max_epochs} epochs...")

model.train(
    max_epochs=max_epochs,
    early_stopping=True,
    early_stopping_patience=20,
    check_val_every_n_epoch=5,
    plan_kwargs={"lr": 1e-3, "reduce_lr_on_plateau": True}
)

# Get latent representation
adata.obsm["X_scvi"] = model.get_latent_representation()
adata_hvg.obsm["X_scvi"] = adata.obsm["X_scvi"]

print("\n" + "=" * 84)
print("STEP 6 — UMAP and clustering")
print("=" * 84)

sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

for res in [0.3, 0.5, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
    print(f"  Leiden res={res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")

print("\n" + "=" * 84)
print("STEP 7 — Gene set scores")
print("=" * 84)

# Gene set scores (cell cycle already computed in Step 4)
avail = set(adata.var_names)
ox = [g for g in OXPHOS_GENES if g in avail]
bcr = [g for g in BCR_GENES if g in avail]

if ox:
    sc.tl.score_genes(adata, gene_list=ox, score_name='oxphos_score')
    print(f"  OXPHOS score: {len(ox)} genes")
if bcr:
    sc.tl.score_genes(adata, gene_list=bcr, score_name='bcr_score')
    print(f"  BCR score: {len(bcr)} genes")

print("\n" + "=" * 84)
print("STEP 8 — CytoTRACE2")
print("=" * 84)

# Export counts for CT2
stream_counts_to_ct2_txt(adata, CT2_INPUT_TXT, collapse_duplicates=False)

# Run CT2
ct2_obj = run_ct2_python(CT2_INPUT_TXT, species="mouse", outdir=CT2_OUTDIR)

if ct2_obj is None:
    cmd = f"cytotrace2 -f {shlex.quote(str(CT2_INPUT_TXT))} -sp mouse --output-dir {shlex.quote(str(CT2_OUTDIR))} --disable-plotting"
    print("[ct2] CLI:", cmd)
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        print(f"[ct2] WARNING: CytoTRACE2 CLI failed with exit code {ret.returncode}")
        print("[ct2] Continuing without CT2 scores...")
    else:
        try:
            df_scores = parse_ct2_scores(CT2_OUTDIR, adata)
            for col in df_scores.columns:
                adata.obs[col] = df_scores[col].values
            print("✓ CytoTRACE2 scores attached")
        except Exception as e:
            print(f"[ct2] Could not parse scores: {e}")
else:
    try:
        # CT2 Python API can return either AnnData or DataFrame
        if hasattr(ct2_obj, 'obs'):
            # AnnData object
            obs = ct2_obj.obs.copy()
        elif isinstance(ct2_obj, pd.DataFrame):
            # DataFrame returned directly
            obs = ct2_obj.copy()
        else:
            raise TypeError(f"Unexpected CT2 return type: {type(ct2_obj)}")
        
        obs.index = obs.index.astype(str)
        obs = _standardize_cols(obs)
        idx = adata.obs_names.astype(str)
        
        sname = next((c for c in obs.columns if ("cytotrace2" in c and "score" in c)), None)
        if sname:
            adata.obs["cytotrace2_score"] = _coerce_float_series(obs.reindex(idx)[sname]).values
        
        p = next((c for c in obs.columns if "potency" in c), None)
        if p:
            adata.obs["cytotrace2_potency"] = obs.reindex(idx)[p].astype(str).values
        
        print("✓ CytoTRACE2 scores attached")
    except Exception as e:
        print(f"[ct2] Could not extract scores from object: {e}")
        # Fallback to parsing output files
        try:
            df_scores = parse_ct2_scores(CT2_OUTDIR, adata)
            for col in df_scores.columns:
                adata.obs[col] = df_scores[col].values
            print("✓ CytoTRACE2 scores attached (from output files)")
        except Exception as e2:
            print(f"[ct2] Could not parse output files either: {e2}")

print("\n" + "=" * 84)
print("STEP 9 — Generate figures")
print("=" * 84)

# Color palettes
condition_colors = {
    "WT_B_cells": "#2ecc71",
    "Crebbp_B_cells": "#3498db", 
    "Pre_malignant": "#f39c12",
    "Matched_malignant": "#e74c3c",
    "Malignant": "#8e44ad"
}

phase_colors = {
    "G1": "#3498db",
    "S": "#e74c3c", 
    "G2M": "#2ecc71"
}

_save_umap(adata, "condition", "umap_condition.png", "Condition", condition_colors)
_save_umap(adata, "sample_id", "umap_sample_id.png", "Sample ID")
_save_umap(adata, "leiden_0.5", "umap_leiden_0.5.png", "Leiden (res=0.5)")
_save_umap(adata, "leiden_1.0", "umap_leiden_1.0.png", "Leiden (res=1.0)")

# Cell cycle phase plots
if "phase" in adata.obs:
    _save_umap(adata, "phase", "umap_cell_cycle_phase.png", "Cell Cycle Phase", phase_colors)
    
    # Additional cell cycle score plots
    if "S_score" in adata.obs:
        _save_umap(adata, "S_score", "umap_S_score.png", "S Phase Score", cmap="RdYlBu_r")
    if "G2M_score" in adata.obs:
        _save_umap(adata, "G2M_score", "umap_G2M_score.png", "G2M Phase Score", cmap="RdYlBu_r")
    if "cc_difference" in adata.obs:
        _save_umap(adata, "cc_difference", "umap_cc_difference.png", "Cell Cycle Difference (G2M-S)", cmap="RdBu_r")
    
    # Cell cycle phase distribution by condition
    fig, ax = plt.subplots(figsize=(10, 6))
    phase_counts = adata.obs.groupby(['condition', 'phase']).size().unstack(fill_value=0)
    phase_pct = phase_counts.div(phase_counts.sum(axis=1), axis=0) * 100
    phase_pct.plot(kind='bar', stacked=True, ax=ax, color=[phase_colors.get(p, '#999999') for p in phase_pct.columns])
    ax.set_ylabel("Percentage of cells")
    ax.set_xlabel("Condition")
    ax.set_title("Cell Cycle Phase Distribution by Condition")
    ax.legend(title="Phase", bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(FIGDIR / "cell_cycle_phase_by_condition.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  ✓ cell_cycle_phase_by_condition.png")
    
    # Cell cycle phase distribution by sample
    fig, ax = plt.subplots(figsize=(14, 6))
    phase_counts_sample = adata.obs.groupby(['sample_id', 'phase']).size().unstack(fill_value=0)
    phase_pct_sample = phase_counts_sample.div(phase_counts_sample.sum(axis=1), axis=0) * 100
    phase_pct_sample.plot(kind='bar', stacked=True, ax=ax, color=[phase_colors.get(p, '#999999') for p in phase_pct_sample.columns])
    ax.set_ylabel("Percentage of cells")
    ax.set_xlabel("Sample")
    ax.set_title("Cell Cycle Phase Distribution by Sample")
    ax.legend(title="Phase", bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(FIGDIR / "cell_cycle_phase_by_sample.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  ✓ cell_cycle_phase_by_sample.png")

if "oxphos_score" in adata.obs:
    _save_umap(adata, "oxphos_score", "umap_oxphos_score.png", "OXPHOS Score", cmap="RdYlBu_r")

if "bcr_score" in adata.obs:
    _save_umap(adata, "bcr_score", "umap_bcr_score.png", "BCR Score", cmap="RdYlBu_r")

# Per-sample UMAPs (highlighted with background in gray)
save_umaps_by_sample(adata, FIGDIR, group_key="sample_id")
save_umaps_by_sample(adata, FIGDIR, group_key="condition")

# Individual sample UMAPs (only that sample's cells, no background)
# Include phase in the color_by list
individual_color_by = ['leiden_0.5', 'condition']
if 'phase' in adata.obs.columns:
    individual_color_by.append('phase')
if 'S_score' in adata.obs.columns:
    individual_color_by.append('S_score')
if 'G2M_score' in adata.obs.columns:
    individual_color_by.append('G2M_score')
if 'cytotrace2_score' in adata.obs.columns:
    individual_color_by.append('cytotrace2_score')

save_umaps_individual_samples(adata, FIGDIR, group_key="sample_id", color_by=individual_color_by)
save_umaps_individual_samples(adata, FIGDIR, group_key="condition", color_by=individual_color_by)

# CT2 plots
if "cytotrace2_score" in adata.obs:
    add_ct2_umap_and_plots(adata, FIGDIR_CT2, title_suffix="(mouse)")

# QC violin plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sc.pl.violin(adata, keys="n_genes_by_counts", groupby="condition", ax=axes[0], show=False)
sc.pl.violin(adata, keys="total_counts", groupby="condition", ax=axes[1], show=False)
sc.pl.violin(adata, keys="pct_counts_mt", groupby="condition", ax=axes[2], show=False)
plt.tight_layout()
fig.savefig(FIGDIR / "qc_violins_by_condition.png", dpi=300, bbox_inches="tight")
plt.close()
print("  ✓ qc_violins_by_condition.png")

print("\n" + "=" * 84)
print("STEP 10 — Save outputs")
print("=" * 84)

# Save scVI model
try:
    model.save(OUTDIR / "scvi_model", overwrite=True)
    print(f"  ✓ scVI model → {OUTDIR / 'scvi_model'}")
except Exception as e:
    print(f"  (warn) Could not save scVI model: {e}")

# Save AnnData
adata.write(OUTDIR / "mouse_integrated.h5ad")
print(f"  ✓ AnnData → {OUTDIR / 'mouse_integrated.h5ad'}")

# Save cell metadata
adata.obs.to_csv(OUTDIR / "cell_metadata.csv")
print(f"  ✓ Metadata → {OUTDIR / 'cell_metadata.csv'}")

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTDIR}")
print(f"  Figures: {FIGDIR}")
print(f"  CT2 results: {CT2_OUTDIR}")
print(f"  Total cells: {adata.n_obs:,}")
print(f"  Total genes: {adata.n_vars:,}")
print("\nDONE.\n")



