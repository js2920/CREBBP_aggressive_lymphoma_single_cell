#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CytoTRACE2 Analysis for Human DLBCL + Mouse Malignant + Tonsil Integration
===========================================================================

MEMORY-OPTIMIZED VERSION: Streams data to CT2 input file incrementally

This script:
1. Loads the integrated scVI data (for cell IDs and UMAP coords)
2. Streams source data to CT2 input file (never loads all genes at once)
3. Runs CytoTRACE2
4. Generates publication-quality visualizations

Author: J
Date: 2025-12-11
"""

import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"

import gc
import re
import shlex
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Set, Optional
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
np.random.seed(0)

# ============================== CONFIGURATION ================================
# Input: Integrated scVI object (for cell IDs and UMAP coords)
INPUT_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Source data files
SOURCE_MOUSE_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/cellbender_filtered")
SOURCE_MOUSE_FILES = {
    "SIGAA3_Matched_malignant_R1": SOURCE_MOUSE_DIR / "SIGAA3_Matched_malignant_R1_GEX_cellbender_filtered.h5",
    "SIGAA4_Matched_malignant_R2": SOURCE_MOUSE_DIR / "SIGAA4_Matched_malignant_R2_GEX_cellbender_filtered.h5",
    "SIGAD5_Malignant_R2": SOURCE_MOUSE_DIR / "SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5",
    "SIGAH1_Malignant_R1": SOURCE_MOUSE_DIR / "SIGAH1_Malignant_R1_GEX_cellbender_filtered.h5",
}
SOURCE_DLBCL_FILES = [
    Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL1_raw.h5ad"),
    Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL2_raw.h5ad"),
    Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL3_raw.h5ad"),
    Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad"),
]
SOURCE_TONSIL_GC = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_GCBC_RNA.h5ad")
SOURCE_TONSIL_MBC = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad")

# Output directory
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/CytoTRACE2")
FIGDIR = OUTPUT_DIR / "figures"
CT2_WORKDIR = OUTPUT_DIR / "ct2_io"
CT2_OUTDIR = OUTPUT_DIR / "ct2_results"

for p in [OUTPUT_DIR, FIGDIR, CT2_WORKDIR, CT2_OUTDIR]:
    p.mkdir(parents=True, exist_ok=True)

CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"

# ============================== HELPER FUNCTIONS =============================
def _is_intlike(mat, n_check=200000) -> bool:
    """Check if matrix contains integer-like values."""
    if sp.issparse(mat):
        data = mat.data[:min(n_check, mat.data.size)]
    else:
        flat = np.ravel(mat)
        data = flat[:min(n_check, flat.size)]
    return data.size > 0 and np.all((data >= 0) & np.isclose(data, np.round(data)))


def get_counts_matrix(adata: ad.AnnData):
    """Get raw counts matrix from AnnData."""
    if "counts" in adata.layers and _is_intlike(adata.layers["counts"]):
        return adata.layers["counts"]
    if adata.raw is not None and _is_intlike(adata.raw.X):
        return adata.raw.X
    if _is_intlike(adata.X):
        return adata.X
    return adata.X  # Fallback


def try_import_cytotrace2():
    """Try to import CytoTRACE2."""
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
    print("[ct2] Running via Python API...")
    try:
        return fn(str(input_txt), species=species, output_dir=str(outdir))
    except TypeError:
        try:
            return fn(str(input_txt))
        except Exception as e:
            print(f"[ct2] Python API failed: {e}")
            return None


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


def parse_ct2_scores(outdir: Path, cell_names: List[str]) -> pd.DataFrame:
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
        
        id_col = next((c for c in ("cell", "cell_id", "barcode") if c in df.columns), df.columns[0])
        score_cols = [c for c in df.columns if ("cytotrace2" in c and "score" in c)]
        
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
        
        n_match = len(set(tmp.index) & set(cell_names))
        if best is None or n_match > best[0]:
            best = (n_match, f, tmp)
    
    if best is None:
        raise FileNotFoundError(f"[ct2] No results found under {outdir}")
    
    print(f"[ct2] Using: {best[1]} (matched {best[0]:,} cells)")
    df_best = best[2]
    df_best = df_best[~df_best.index.duplicated(keep="first")]
    return df_best.reindex(cell_names)


# ============================== STREAMING CT2 WRITER =========================
class StreamingCT2Writer:
    """
    Memory-efficient CT2 matrix writer.
    Collects gene expression per-cell across multiple source files,
    then writes genes × cells matrix.
    """
    
    def __init__(self, target_cells: List[str], output_path: Path):
        self.target_cells = target_cells
        self.cell_to_idx = {c: i for i, c in enumerate(target_cells)}
        self.output_path = output_path
        
        # Gene -> array of counts (sparse storage)
        self.gene_data: Dict[str, np.ndarray] = {}
        self.n_cells = len(target_cells)
        
    def add_dataset(self, adata: ad.AnnData, tag: str):
        """Add counts from an AnnData object for matching cells."""
        # Find matching cells
        adata_cells = set(adata.obs_names.astype(str))
        matching = [c for c in self.target_cells if c in adata_cells]
        
        if not matching:
            print(f"    [{tag}] No matching cells found")
            return 0
        
        print(f"    [{tag}] Processing {len(matching):,} matching cells...")
        
        # Get counts matrix
        X = get_counts_matrix(adata)
        is_sparse = sp.issparse(X)
        
        # Create cell index mapping for this dataset
        adata_cell_list = list(adata.obs_names.astype(str))
        local_idx = {c: i for i, c in enumerate(adata_cell_list)}
        
        # Gene names (uppercase)
        gene_names = [str(g).upper() for g in adata.var_names]
        
        # For each gene, accumulate counts
        for j, gene in enumerate(gene_names):
            if gene in ("", "N/A", "NA", "NULL", "NONE", "NAN"):
                continue
            
            # Get column
            col = X[:, j]
            if is_sparse:
                arr = col.toarray().ravel()
            else:
                arr = np.asarray(col).ravel()
            
            # Initialize gene array if needed
            if gene not in self.gene_data:
                self.gene_data[gene] = np.zeros(self.n_cells, dtype=np.float32)
            
            # Fill in counts for matching cells
            for cell in matching:
                target_idx = self.cell_to_idx[cell]
                local_cell_idx = local_idx[cell]
                self.gene_data[gene][target_idx] = arr[local_cell_idx]
            
            if (j + 1) % 5000 == 0:
                print(f"      [{tag}] {j+1:,}/{len(gene_names):,} genes", end="\r")
        
        print(f"    [{tag}] Done: {len(matching):,} cells, {len(gene_names):,} genes")
        return len(matching)
    
    def write(self, min_cells_expressed: int = 10):
        """Write the CT2 input file."""
        print(f"\n[ct2-io] Writing CT2 input: {self.output_path}")
        print(f"[ct2-io] Total genes collected: {len(self.gene_data):,}")
        print(f"[ct2-io] Total cells: {self.n_cells:,}")
        
        # Filter genes by expression
        genes_to_write = []
        for gene, counts in self.gene_data.items():
            n_expressed = np.sum(counts > 0)
            if n_expressed >= min_cells_expressed:
                genes_to_write.append(gene)
        
        print(f"[ct2-io] Genes expressed in ≥{min_cells_expressed} cells: {len(genes_to_write):,}")
        
        # Write file
        with open(self.output_path, "w", buffering=1024*1024) as fh:
            # Header: gene + cell names
            fh.write("gene\t" + "\t".join(self.target_cells) + "\n")
            
            wrote = 0
            for gene in sorted(genes_to_write):
                counts = self.gene_data[gene]
                
                # Format as integers if possible
                if np.all(np.isclose(counts, np.round(counts))):
                    line = gene + "\t" + "\t".join(str(int(v)) for v in counts)
                else:
                    line = gene + "\t" + "\t".join(f"{v:.6g}" for v in counts)
                
                fh.write(line + "\n")
                wrote += 1
                
                if wrote % 2000 == 0:
                    print(f"  [ct2-io] {wrote:,}/{len(genes_to_write):,} genes written", end="\r")
        
        print(f"\n[ct2-io] ✓ Wrote {wrote:,} genes × {self.n_cells:,} cells")
        
        # Clear memory
        self.gene_data.clear()
        gc.collect()
        
        return wrote


# ============================== MAIN PIPELINE ================================
print("=" * 84)
print("CYTOTRACE2 ANALYSIS: HUMAN DLBCL + MOUSE MALIGNANT + TONSIL")
print("=" * 84)
print(f"Input: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print()

# ==================== STEP 1: Load Integrated Data ====================
print("=" * 84)
print("STEP 1 — Load integrated data (for cell IDs and UMAP)")
print("=" * 84)

print(f"Loading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if "disease_state" in adata.obs:
    print(f"  Disease states: {list(adata.obs['disease_state'].unique())}")
if "species" in adata.obs:
    print(f"  Species: {list(adata.obs['species'].unique())}")

# Get target cell IDs
target_cells = list(adata.obs_names.astype(str))
print(f"  Target cells: {len(target_cells):,}")

# ==================== STEP 2: Stream Source Data to CT2 File ====================
print("\n" + "=" * 84)
print("STEP 2 — Stream source data to CT2 input file (MEMORY EFFICIENT)")
print("=" * 84)

writer = StreamingCT2Writer(target_cells, CT2_INPUT_TXT)
total_matched = 0

# Process mouse samples one at a time
print("\n  Processing mouse malignant samples...")
for sample_id, path in SOURCE_MOUSE_FILES.items():
    if not path.exists():
        print(f"    [{sample_id}] File not found, skipping")
        continue
    
    try:
        # Load
        a = sc.read_10x_h5(str(path))
        a.var_names_make_unique()
        
        # Fix cell names to match integrated object
        a.obs_names = pd.Index([f"{bc.split('-')[0]}-{sample_id}" for bc in a.obs_names])
        
        # Convert mouse genes to human (uppercase)
        a.var_names = pd.Index([str(g).upper() for g in a.var_names])
        a.var_names_make_unique()
        
        # Add to writer
        n = writer.add_dataset(a, sample_id)
        total_matched += n
        
        # Free memory
        del a
        gc.collect()
        
    except Exception as e:
        print(f"    [{sample_id}] Error: {e}")

# Process DLBCL samples one at a time
print("\n  Processing DLBCL samples...")
for path in SOURCE_DLBCL_FILES:
    if not path.exists():
        print(f"    [{path.stem}] File not found, skipping")
        continue
    
    try:
        a = sc.read_h5ad(path)
        a.var_names_make_unique()
        
        # Uppercase gene names
        a.var_names = pd.Index([str(g).upper() for g in a.var_names])
        a.var_names_make_unique()
        
        # Add to writer
        n = writer.add_dataset(a, path.stem)
        total_matched += n
        
        del a
        gc.collect()
        
    except Exception as e:
        print(f"    [{path.stem}] Error: {e}")

# Process tonsil samples one at a time (BACKED MODE)
print("\n  Processing tonsil samples (backed mode for memory efficiency)...")
for path in [SOURCE_TONSIL_GC, SOURCE_TONSIL_MBC]:
    if not path.exists():
        print(f"    [{path.stem}] File not found, skipping")
        continue
    
    try:
        # Load in backed mode
        print(f"    [{path.stem}] Loading in backed mode...")
        a_backed = sc.read_h5ad(path, backed="r")
        
        # Find matching cells
        target_set = set(target_cells)
        backed_cells = set(a_backed.obs_names.astype(str))
        matching = list(target_set & backed_cells)
        
        print(f"    [{path.stem}] Found {len(matching):,} matching cells (from {a_backed.n_obs:,} total)")
        
        if matching:
            # Load ONLY matching cells into memory
            mask = a_backed.obs_names.isin(matching)
            a = a_backed[mask].to_memory()
            
            # Close backed file
            if hasattr(a_backed, 'file') and a_backed.file is not None:
                try:
                    a_backed.file.close()
                except:
                    pass
            del a_backed
            gc.collect()
            
            # Uppercase gene names
            a.var_names = pd.Index([str(g).upper() for g in a.var_names])
            a.var_names_make_unique()
            
            # Add to writer
            n = writer.add_dataset(a, path.stem)
            total_matched += n
            
            del a
            gc.collect()
        else:
            if hasattr(a_backed, 'file') and a_backed.file is not None:
                try:
                    a_backed.file.close()
                except:
                    pass
            del a_backed
            gc.collect()
            
    except Exception as e:
        print(f"    [{path.stem}] Error: {e}")
        import traceback
        traceback.print_exc()

print(f"\n  Total matched cells across sources: {total_matched:,}")
print(f"  Target cells: {len(target_cells):,}")

# Write CT2 input file
min_cells = max(10, int(len(target_cells) * 0.001))  # At least 0.1% or 10 cells
n_genes = writer.write(min_cells_expressed=min_cells)

# ==================== STEP 3: Run CytoTRACE2 ====================
print("\n" + "=" * 84)
print("STEP 3 — Run CytoTRACE2")
print("=" * 84)

print("  NOTE: CytoTRACE2 uses species='human' since all genes are in human symbols")
print("        (Mouse genes were converted via uppercase)")
print()

ct2_obj = run_ct2_python(CT2_INPUT_TXT, species="human", outdir=CT2_OUTDIR)

if ct2_obj is None:
    print("[ct2] Trying CLI fallback...")
    cmd = f"cytotrace2 -f {shlex.quote(str(CT2_INPUT_TXT))} -sp human --output-dir {shlex.quote(str(CT2_OUTDIR))} --disable-plotting"
    print(f"[ct2] CLI command: {cmd}")
    try:
        ret = subprocess.run(cmd, shell=True, timeout=7200)
        if ret.returncode != 0:
            print(f"[ct2] WARNING: CT2 CLI failed with exit code {ret.returncode}")
    except subprocess.TimeoutExpired:
        print("[ct2] WARNING: CT2 CLI timed out after 2 hours")
    except Exception as e:
        print(f"[ct2] WARNING: CT2 CLI failed: {e}")

# ==================== STEP 4: Parse and attach scores ====================
print("\n" + "=" * 84)
print("STEP 4 — Parse CytoTRACE2 results")
print("=" * 84)

df_scores = None

# Try to parse from object
if ct2_obj is not None:
    try:
        if hasattr(ct2_obj, 'obs'):
            obs = ct2_obj.obs.copy()
        elif isinstance(ct2_obj, pd.DataFrame):
            obs = ct2_obj.copy()
        else:
            raise TypeError(f"Unexpected CT2 return type: {type(ct2_obj)}")
        
        obs.index = obs.index.astype(str)
        obs = _standardize_cols(obs)
        
        df_scores = pd.DataFrame(index=target_cells)
        sname = next((c for c in obs.columns if ("cytotrace2" in c and "score" in c)), None)
        if sname:
            df_scores["cytotrace2_score"] = _coerce_float_series(obs.reindex(target_cells)[sname]).values
        
        pname = next((c for c in obs.columns if "potency" in c), None)
        if pname:
            df_scores["cytotrace2_potency"] = obs.reindex(target_cells)[pname].astype(str).values
            
        rname = next((c for c in obs.columns if "relative" in c), None)
        if rname:
            df_scores["cytotrace2_relative"] = _coerce_float_series(obs.reindex(target_cells)[rname]).values
            
    except Exception as e:
        print(f"[ct2] Could not parse object: {e}")
        df_scores = None

# Try to parse from files
if df_scores is None:
    try:
        df_scores = parse_ct2_scores(CT2_OUTDIR, target_cells)
    except FileNotFoundError as e:
        print(f"[ct2] WARNING: {e}")
        print("[ct2] CytoTRACE2 may have failed. Check output directory.")

# Attach scores to adata
if df_scores is not None and not df_scores.empty:
    for col in df_scores.columns:
        adata.obs[col] = df_scores[col].values
    print("✓ CytoTRACE2 scores attached to adata.obs")
    
    # Save updated adata
    output_h5ad = OUTPUT_DIR / "integrated_with_cytotrace2.h5ad"
    adata.write_h5ad(output_h5ad)
    print(f"✓ Saved: {output_h5ad}")
else:
    print("⚠ CytoTRACE2 scores NOT attached")
    print("  Check CT2 output directory for errors")

# ==================== STEP 5: Generate Figures ====================
print("\n" + "=" * 84)
print("STEP 5 — Generate figures")
print("=" * 84)

if "cytotrace2_score" not in adata.obs:
    print("  Skipping figures (no CT2 scores)")
else:
    # Color palettes
    disease_colors = {
        "Mouse_Malignant": "#FF6B6B",
        "Mouse_Matched_malignant": "#FF9999",
        "DLBCL": "#8B0000",
        "Tonsil_GC_B": "#4169E1",
        "Tonsil_Normal": "#4169E1"
    }
    species_colors = {"mouse": "#98FB98", "human": "#6495ED"}
    
    def save_fig(fig, name):
        for fmt in ["png", "pdf", "svg"]:
            fig.savefig(FIGDIR / f"{name}.{fmt}", dpi=300, bbox_inches="tight")
        print(f"  ✓ {name}")
    
    # 1. UMAP by CT2 score
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=20, frameon=False,
               cmap="viridis", title="CytoTRACE2 Score\n(Higher = Less Differentiated)")
    save_fig(fig, "umap_cytotrace2_score")
    plt.close()
    
    # 2. UMAP by CT2 score (reversed colormap)
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=20, frameon=False,
               cmap="viridis_r", title="CytoTRACE2 Score (Reversed)\n(Darker = Less Differentiated)")
    save_fig(fig, "umap_cytotrace2_score_reversed")
    plt.close()
    
    # 3. CT2 score histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    scores = adata.obs["cytotrace2_score"].dropna()
    ax.hist(scores, bins=50, edgecolor="black", alpha=0.7, color="steelblue")
    ax.axvline(scores.median(), color="red", linestyle="--", label=f"Median: {scores.median():.3f}")
    ax.axvline(scores.mean(), color="green", linestyle="--", label=f"Mean: {scores.mean():.3f}")
    ax.set_xlabel("CytoTRACE2 Score", fontsize=12)
    ax.set_ylabel("Cell Count", fontsize=12)
    ax.set_title("CytoTRACE2 Score Distribution", fontsize=14)
    ax.legend()
    save_fig(fig, "histogram_cytotrace2_score")
    plt.close()
    
    # 4. CT2 by disease state (violin)
    if "disease_state" in adata.obs:
        fig, ax = plt.subplots(figsize=(12, 6))
        order = sorted(adata.obs["disease_state"].unique())
        colors = [disease_colors.get(d, "#808080") for d in order]
        sns.violinplot(data=adata.obs, x="disease_state", y="cytotrace2_score",
                       order=order, palette=colors, ax=ax)
        ax.set_xlabel("Disease State", fontsize=12)
        ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
        ax.set_title("Differentiation Potential by Disease State\n(Higher = Less Differentiated)", fontsize=14)
        plt.xticks(rotation=45, ha="right")
        save_fig(fig, "violin_cytotrace2_by_disease")
        plt.close()
    
    # 5. CT2 by species (violin)
    if "species" in adata.obs:
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.violinplot(data=adata.obs, x="species", y="cytotrace2_score",
                       palette=species_colors, ax=ax)
        ax.set_xlabel("Species", fontsize=12)
        ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
        ax.set_title("Differentiation Potential by Species", fontsize=14)
        save_fig(fig, "violin_cytotrace2_by_species")
        plt.close()
    
    # 6. CT2 by sample (boxplot)
    if "sample_batch" in adata.obs:
        fig, ax = plt.subplots(figsize=(14, 6))
        order = sorted(adata.obs["sample_batch"].unique())
        sns.boxplot(data=adata.obs, x="sample_batch", y="cytotrace2_score",
                    order=order, ax=ax)
        ax.set_xlabel("Sample", fontsize=12)
        ax.set_ylabel("CytoTRACE2 Score", fontsize=12)
        ax.set_title("Differentiation Potential by Sample", fontsize=14)
        plt.xticks(rotation=45, ha="right")
        save_fig(fig, "boxplot_cytotrace2_by_sample")
        plt.close()
    
    # 7. CT2 potency category distribution
    if "cytotrace2_potency" in adata.obs:
        fig, ax = plt.subplots(figsize=(10, 6))
        potency_counts = adata.obs["cytotrace2_potency"].value_counts()
        potency_counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="black")
        ax.set_xlabel("Potency Category", fontsize=12)
        ax.set_ylabel("Cell Count", fontsize=12)
        ax.set_title("CytoTRACE2 Potency Distribution", fontsize=14)
        plt.xticks(rotation=45, ha="right")
        save_fig(fig, "barplot_cytotrace2_potency")
        plt.close()
        
        # Potency by disease state
        if "disease_state" in adata.obs:
            crosstab = pd.crosstab(adata.obs["disease_state"], adata.obs["cytotrace2_potency"])
            crosstab_pct = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
            
            fig, ax = plt.subplots(figsize=(12, 8))
            crosstab_pct.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
            ax.set_xlabel("Disease State", fontsize=12)
            ax.set_ylabel("Percentage", fontsize=12)
            ax.set_title("Potency Distribution by Disease State", fontsize=14)
            ax.legend(title="Potency", bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            save_fig(fig, "stacked_bar_potency_by_disease")
            plt.close()
            
            # Save crosstab
            crosstab.to_csv(OUTPUT_DIR / "potency_by_disease_crosstab.csv")
            crosstab_pct.to_csv(OUTPUT_DIR / "potency_by_disease_percentage.csv")
    
    # 8. Split UMAP by species
    if "species" in adata.obs:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        umap_coords = adata.obsm["X_umap"]
        vmin, vmax = adata.obs["cytotrace2_score"].quantile([0.01, 0.99])
        
        for idx, species_name in enumerate(["human", "mouse"]):
            if species_name not in adata.obs["species"].values:
                continue
            ax = axes[idx]
            mask = adata.obs["species"] == species_name
            
            scatter = ax.scatter(
                umap_coords[mask, 0], umap_coords[mask, 1],
                c=adata.obs.loc[mask, "cytotrace2_score"],
                cmap="viridis", s=15, alpha=0.7, vmin=vmin, vmax=vmax,
                rasterized=True
            )
            ax.set_title(f"{species_name.upper()} cells: CytoTRACE2 Score\n({mask.sum():,} cells)", fontsize=12)
            ax.axis("off")
            plt.colorbar(scatter, ax=ax, label="CT2 Score")
        
        plt.tight_layout()
        save_fig(fig, "umap_cytotrace2_split_by_species")
        plt.close()

# ==================== STEP 6: Statistics ====================
print("\n" + "=" * 84)
print("STEP 6 — Generate statistics")
print("=" * 84)

if "cytotrace2_score" in adata.obs:
    scores = adata.obs["cytotrace2_score"].dropna()
    
    print(f"\nOverall CT2 statistics:")
    print(f"  Mean: {scores.mean():.4f}")
    print(f"  Median: {scores.median():.4f}")
    print(f"  Std: {scores.std():.4f}")
    print(f"  Min: {scores.min():.4f}")
    print(f"  Max: {scores.max():.4f}")
    
    # By disease state
    if "disease_state" in adata.obs:
        print("\nCT2 by disease state:")
        stats_disease = adata.obs.groupby("disease_state")["cytotrace2_score"].agg(
            ["mean", "median", "std", "count"]
        ).round(4)
        print(stats_disease.to_string())
        stats_disease.to_csv(OUTPUT_DIR / "cytotrace2_by_disease_state.csv")
    
    # By species
    if "species" in adata.obs:
        print("\nCT2 by species:")
        stats_species = adata.obs.groupby("species")["cytotrace2_score"].agg(
            ["mean", "median", "std", "count"]
        ).round(4)
        print(stats_species.to_string())
        stats_species.to_csv(OUTPUT_DIR / "cytotrace2_by_species.csv")
    
    # By sample
    if "sample_batch" in adata.obs:
        stats_sample = adata.obs.groupby("sample_batch")["cytotrace2_score"].agg(
            ["mean", "median", "std", "count"]
        ).round(4)
        stats_sample.to_csv(OUTPUT_DIR / "cytotrace2_by_sample.csv")

# ==================== DONE ====================
print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTPUT_DIR}")
print(f"  Figures: {FIGDIR}")
print(f"  CT2 results: {CT2_OUTDIR}")

if "cytotrace2_score" in adata.obs:
    print(f"\n  Key outputs:")
    print(f"    - integrated_with_cytotrace2.h5ad")
    print(f"    - cytotrace2_by_disease_state.csv")
    print(f"    - cytotrace2_by_species.csv")
    print(f"    - umap_cytotrace2_score.png/pdf/svg")
    print(f"    - violin_cytotrace2_by_disease.png/pdf/svg")

print("\nDONE.\n")


