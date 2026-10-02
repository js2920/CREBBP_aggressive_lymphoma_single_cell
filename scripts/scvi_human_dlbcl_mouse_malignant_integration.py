#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Human DLBCL + Malignant Mouse + Tonsil GC B Integration (scVI + CytoTRACE2)
===========================================================================

Integrates:
- Human DLBCL samples (Roider DLBCL1/2/3 + Alizadeh CD20+)
- Malignant mouse samples (ONLY Malignant + Matched_malignant, CellBender filtered)
- Human Tonsil GC B cells (LZ/DZ + centroblasts + centrocytes, with optional proliferation filtering)

Key features:
- Mouse→Human ortholog mapping via BioMart 1:1
- Proper Scrublet doublet detection for mouse samples
- Tonsil GC B cell filtering with proliferation removal
- scVI integration with species covariate
- CytoTRACE2 analysis (uses UNION of all genes, as recommended by CT2)
- Cross-species HVG filtering

Author: J
Date: 2025-12-10
"""

# ============================== SETUP ========================================
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

import gc
import re
import shlex
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import urllib.request
import urllib.parse

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
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration")
FIGDIR = OUTDIR / "figures"
FIGDIR_CT2 = OUTDIR / "figures_ct2"
CT2_WORKDIR = OUTDIR / "ct2_io"
CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"
CT2_OUTDIR = OUTDIR / "cytotrace2_results"

for p in (OUTDIR, FIGDIR, FIGDIR_CT2, CT2_WORKDIR, CT2_OUTDIR):
    p.mkdir(parents=True, exist_ok=True)

# Human DLBCL files
DLBCL_FILES = [
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL1_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL2_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL3_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad",
]

# Human Tonsil data files
TONSIL_GC_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_GCBC_RNA.h5ad"
TONSIL_MBC_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad"
MAX_TONSIL_CELLS = 12500  # Maximum total GC + Memory B cells to keep

# Proliferation filtering for tonsil cells
FILTER_PROLIFERATING_TONSIL = False  # Changed to False to preserve Dark Zone (Centroblasts)
PROLIFERATION_THRESHOLD = 0.20  # Relaxed threshold if enabled (was 0.10)
EXCLUDE_PROLIFERATIVE_GC_TYPES = False  # Changed to False to preserve Dark Zone annotations

# Toggle CytoTRACE2 (Set to True later to run CT2 after scVI completes)
RUN_CYTOTRACE2 = False  # Set to False to skip CT2 and save memory

# Mouse malignant samples (CellBender filtered)
CELLBENDER_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/cellbender_filtered")
MOUSE_MALIGNANT_FILES = {
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
    "SIGAD5_Malignant_R2": {
        "path": CELLBENDER_DIR / "SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R2"
    },
    "SIGAH1_Malignant_R1": {
        "path": CELLBENDER_DIR / "SIGAH1_Malignant_R1_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R1"
    },
}

# QC thresholds
MT_THRESHOLD = 10.0
MIN_GENES = 200
MIN_CELLS = 3

print("=" * 84)
print("HUMAN DLBCL + MALIGNANT MOUSE + TONSIL GC B INTEGRATION (scVI + CytoTRACE2)")
print("=" * 84)
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"Output dir    : {OUTDIR}")
print(f"MT threshold  : {MT_THRESHOLD}%")
print(f"HVGs for integration: 4000")
print(f"Tonsil B cell types: LZ/DZ + centroblasts + centrocytes + GC B cells + memory B cells")
print(f"Max tonsil cells: {MAX_TONSIL_CELLS:,}")
if FILTER_PROLIFERATING_TONSIL:
    print(f"Proliferation filtering: ON (S+G2M score threshold: {PROLIFERATION_THRESHOLD})")
    print(f"  -> WARNING: This may remove Dark Zone (Centroblast) populations.")
else:
    print(f"Proliferation filtering: OFF (Preserving Dark Zone/Centroblasts)")
if EXCLUDE_PROLIFERATIVE_GC_TYPES:
    print(f"Exclude proliferative GC types: ON")
else:
    print(f"Exclude proliferative GC types: OFF")
print()

# ============================== GENE SETS ====================================
OXPHOS_GENES = [
    'COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A2','COX7B','COX7C','COX8A',
    'CYC1','CYCS','NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9',
    'NDUFA10','NDUFA11','NDUFA12','NDUFA13','NDUFAB1','NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5',
    'NDUFB6','NDUFB7','NDUFB8','NDUFB9','NDUFB10','NDUFB11','NDUFC1','NDUFC2','NDUFS1','NDUFS2',
    'NDUFS3','NDUFS4','NDUFS5','NDUFS6','NDUFS7','NDUFS8','NDUFV1','NDUFV2','NDUFV3','SDHA','SDHB',
    'SDHC','SDHD','UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ',
    'ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5MC1','ATP5MC2','ATP5MC3','ATP5ME','ATP5MF',
    'ATP5MG','ATP5PB','ATP5PD','ATP5PF','ATP5PO'
]
BCR_GENES = [
    'CD79A','CD79B','CD19','CD22','CD72','CR2','FCRL1','FCRL2','FCRL3','FCRL4','FCRL5','MS4A1',
    'IGHM','IGHD','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHE','BTK','LYN','SYK','BLK',
    'BLNK','PIK3CD','PIK3AP1','PLCG2','PRKCB','NFKB1','NFKB2','REL','RELA','NFATC1','NFATC2',
    'BCL10','CARD11','MALT1','MAP3K7','IKBKB','IKBKG','CHUK','PTPN6','PTPRC','VAV1','VAV2',
    'VAV3','GRB2','SOS1','SOS2','HRAS','KRAS','NRAS','RAF1','MAP2K1','MAP2K2','MAPK1','MAPK3'
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


def enforce_raw_counts(adata: ad.AnnData, prefer_layers=("counts", "raw_counts"), tag=""):
    """Ensure layers['counts'] contains raw integer counts."""
    for lyr in prefer_layers:
        if lyr in adata.layers and _is_intlike(adata.layers[lyr]):
            adata.layers["counts"] = adata.layers[lyr]
            print(f"  [{tag}] Using '{lyr}' as counts.")
            return adata
    if adata.raw is not None and _is_intlike(adata.raw.X):
        adata.layers["counts"] = adata.raw.X
        print(f"  [{tag}] Copied raw.X as counts.")
        return adata
    if _is_intlike(adata.X):
        adata.layers["counts"] = adata.X.copy()
        print(f"  [{tag}] Using X as counts.")
        return adata
    raise ValueError(f"[{tag}] No integer-like counts found for scVI.")


def harmonize_preprocessing(adata, target_sum=1e4, tag=""):
    """Normalize and log-transform, preserving counts."""
    print(f"  [{tag}] normalize_total→log1p (counts preserved)")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["normalized"] = adata.X.copy()
    return adata


# ============================== ORTHOLOG MAPPING =============================
def fetch_biomart_m2h_one2one(host: str = "https://www.ensembl.org") -> Dict[str, str]:
    """BioMart mouse→human 1:1 via martservice XML."""
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="mmusculus_gene_ensembl" interface="default">
    <Filter name="with_hsapiens_homolog" excluded="0"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
    <Attribute name="hsapiens_homolog_orthology_type"/>
  </Dataset>
</Query>"""
    url = f"{host}/biomart/martservice?query=" + urllib.parse.quote(xml)
    print("  [ortholog] Querying Ensembl BioMart for 1:1 orthologs ...")
    txt = urllib.request.urlopen(url, timeout=120).read().decode("utf-8")
    m2h = {}
    for line in txt.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        mm, hs, typ = (p.strip() for p in parts)
        if typ == "ortholog_one2one" and mm and hs and mm not in m2h:
            m2h[mm] = hs.upper()
    if len(m2h) < 1000:
        raise RuntimeError("BioMart map unexpectedly small")
    print(f"  [ortholog] BioMart 1:1 pairs: {len(m2h):,}")
    return m2h


def load_one_to_one_orthologs_mgi(mgi_file: Optional[Path]) -> Dict[str, str]:
    """Load MGI ortholog mapping as fallback."""
    if mgi_file is None or not mgi_file.exists():
        return {}
    df = pd.read_csv(mgi_file, sep="\t", dtype=str)
    out: Dict[str, str] = {}
    for _, g in df.groupby("DB Class Key"):
        m = g[g["Common Organism Name"] == "mouse, laboratory"]["Symbol"].dropna().unique()
        h = g[g["Common Organism Name"] == "human"]["Symbol"].dropna().unique()
        if len(m) == 1 and len(h) == 1:
            out[m[0]] = h[0].upper()
    print(f"  [ortholog] MGI 1:1 pairs: {len(out):,}")
    return out


def get_m2h_map(prefer_biomart=True) -> Dict[str, str]:
    """Get mouse-to-human ortholog mapping."""
    m2h = {}
    if prefer_biomart:
        try:
            m2h.update(fetch_biomart_m2h_one2one())
        except Exception as e:
            print(f"  [ortholog] BioMart failed: {e}")
    rpt = OUTDIR / "HOM_MouseHumanSequence.rpt"
    if not rpt.exists():
        try:
            print("  [ortholog] Downloading MGI HOM report ...")
            urllib.request.urlretrieve(
                "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                str(rpt)
            )
        except Exception as e:
            print(f"  [ortholog] MGI download failed: {e}")
            return m2h
    mgi_map = load_one_to_one_orthologs_mgi(rpt)
    for k, v in mgi_map.items():
        m2h.setdefault(k, v)
    print(f"  [ortholog] total unique pairs: {len(m2h):,}")
    return m2h


def convert_mouse_to_human_genes(adata: ad.AnnData, m2h: Dict[str, str], tag=""):
    """Convert mouse gene symbols to human orthologs."""
    print(f"  [{tag}] mouse→human gene symbols (1:1 map → uppercase fallback)")
    mapped = [m2h.get(str(g), str(g).upper()) for g in adata.var_names]
    adata.var_names = pd.Index(mapped)
    adata.var_names_make_unique()
    return adata


def harmonize_human_genes(adata: ad.AnnData, tag=""):
    """Harmonize human gene symbols to uppercase."""
    adata.var_names = pd.Index([str(g).upper() for g in adata.var_names])
    adata.var_names_make_unique()
    return adata


# ============================== TONSIL FILTERING =============================
def find_label_column(adata: ad.AnnData) -> Optional[str]:
    """Find cell type annotation column in tonsil atlas."""
    candidates = [
        "annotation_20230508", "annotation_20220414", "cell_type", "CellType",
        "celltype", "label", "celltype_l2", "celltype_l1", "Azimuth.celltype",
        "celltype_final", "predicted.celltype", "predicted_celltype"
    ]
    for col in candidates:
        if col in adata.obs.columns:
            return col
    # Try regex search
    for col in adata.obs.columns:
        if re.search(r"(cell.?type|annotation|label)", col, re.I):
            return col
    return None


def filter_gc_b_cells(adata: ad.AnnData,
                      filter_proliferating: bool = False,
                      proliferation_threshold: float = 0.10,
                      exclude_proliferative_gc_types: bool = False) -> ad.AnnData:
    """Filter tonsil atlas for germinal center and memory B cell types, optionally filter proliferating cells.

    Includes: LZ/DZ B cells, centroblasts, centrocytes, GC B cells, memory B cells, and related populations.
    """
    # Convert to memory if backed (needed for filtering and copying)
    if adata.isbacked:
        print(f"  [tonsil] Converting backed AnnData to memory...")
        adata = adata.to_memory()

    label_col = find_label_column(adata)
    if label_col is None:
        print(f"  [tonsil] WARNING: No cell type column found!")
        print(f"  [tonsil] Available columns: {list(adata.obs.columns[:10])}")
        print(f"  [tonsil] Using all cells (no filtering)")
        return adata.copy()

    labels = adata.obs[label_col].astype(str).str.lower()

    # Comprehensive pattern matching for germinal center B cell types
    # LZ/DZ patterns
    is_gc = (
        labels.str.contains(r"\b(lz|light\s*zone)\b", case=False, na=False) |
        labels.str.contains(r"\b(dz|dark\s*zone)\b", case=False, na=False) |
        labels.str.contains(r"\bgcb.*(lz|dz)\b", case=False, na=False) |
        labels.str.contains(r"\b(lz|dz).*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"lz[-/]dz", case=False, na=False)
    )

    # Centroblast patterns
    is_gc |= (
        labels.str.contains(r"\bcentroblast", case=False, na=False) |
        labels.str.contains(r"\bcb\s*b\s*cells?\b", case=False, na=False)
    )

    # Centrocyte patterns
    is_gc |= (
        labels.str.contains(r"\bcentrocyte", case=False, na=False) |
        labels.str.contains(r"\bcc\s*b\s*cells?\b", case=False, na=False)
    )

    # General GC B cell patterns
    is_gc |= (
        labels.str.contains(r"\bgc\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bgerminal\s*center\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bgcb\b", case=False, na=False) |
        labels.str.contains(r"\bgc\s*b\s*lymphocyte", case=False, na=False)
    )

    # Memory B cell patterns
    is_gc |= (
        labels.str.contains(r"\bmemory\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bmbc\b", case=False, na=False) |
        labels.str.contains(r"\bmem\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bb\s*memory", case=False, na=False)
    )

    # Activated B cells in GC context (but exclude plasma cells)
    is_gc |= (
        (labels.str.contains(r"\bactivated\s*b\s*cells?\b", case=False, na=False) |
         labels.str.contains(r"\bact\s*b\s*cells?\b", case=False, na=False)) &
        ~labels.str.contains(r"\bplasma\b", case=False, na=False)
    )

    n_matched = is_gc.sum()
    print(f"  [tonsil] Found {n_matched:,} GC + memory B cells (from {adata.n_obs:,} total)")

    if n_matched == 0:
        print(f"  [tonsil] WARNING: No GC/memory B cells found! Available cell types:")
        for ct, count in adata.obs[label_col].value_counts().head(15).items():
            print(f"    {ct}: {count:,} cells")
        print(f"  [tonsil] Using all cells instead")
        filtered = adata.copy()
    else:
        filtered = adata[is_gc].copy()
        # Report breakdown by cell type
        if label_col:
            print(f"  [tonsil] B cell type breakdown:")
            gc_types = filtered.obs[label_col].value_counts()
            for ct, count in gc_types.head(10).items():
                pct = 100 * count / len(filtered)
                print(f"    {ct}: {count:,} cells ({pct:.1f}%)")
        print(f"  [tonsil] Filtered to {filtered.n_obs:,} GC + memory B cells")

        # Print final breakdown
        if label_col:
            print(f"  [tonsil] Final composition:")
            final_types = filtered.obs[label_col].value_counts()
            for ct, count in final_types.head(10).items():
                 pct = 100 * count / len(filtered)
                 print(f"    - {ct}: {count:,} ({pct:.1f}%)")

    # Exclude proliferative GC cell types by annotation (e.g., "DZ late Sphase", "DZ early G2Mphase")
    if exclude_proliferative_gc_types and label_col and filtered.n_obs > 0:
        labels = filtered.obs[label_col].astype(str).str.lower()
        is_prolif_type = (
            labels.str.contains(r"\bproliferative\b", case=False, na=False) |
            labels.str.contains(r"\bs\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bg2m\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bg2\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bm\s*phase\b", case=False, na=False)
        )
        n_prolif_type = is_prolif_type.sum()
        if n_prolif_type > 0:
            filtered = filtered[~is_prolif_type].copy()
            print(f"  [tonsil] Excluded {n_prolif_type:,} cells with proliferative cell type annotations")
            print(f"  [tonsil] Retained {filtered.n_obs:,} non-proliferative B cells")

    # Filter proliferating cells if requested
    if filter_proliferating and filtered.n_obs > 0:
        print(f"  [tonsil] Computing cell cycle scores to filter proliferating cells...")
        # Preserve counts if not already in layers
        if "counts" not in filtered.layers:
            if _is_intlike(filtered.X):
                filtered.layers["counts"] = filtered.X.copy()
            elif filtered.raw is not None and _is_intlike(filtered.raw.X):
                filtered.layers["counts"] = filtered.raw.X.copy()
        # Ensure we have normalized data for cell cycle scoring
        if "normalized" not in filtered.layers:
            sc.pp.normalize_total(filtered, target_sum=1e4)
            sc.pp.log1p(filtered)
            filtered.layers["normalized"] = filtered.X.copy()

        # Cell cycle genes
        cc_genes = [x.strip() for x in """
MCM5,PCNA,TYMS,FEN1,MCM2,MCM4,RRM1,UNG,GINS2,MCM6,CDCA7,DTL,PRIM1,UHRF1,MLF1IP,HELLS,RFC2,RPA2,NASP,RAD51AP1,GMNN,WDR76,SLBP,CCNE2,UBR7,POLD3,MSH2,
ATAD2,RAD51,RRM2,CDC45,CDC6,EXO1,TIPIN,DSCC1,BLM,CASP8AP2,USP1,CLSPN,POLA1,CHAF1B,BRIP1,E2F8,
HMGB2,CDK1,NUSAP1,UBE2C,BIRC5,TPX2,TOP2A,NDC80,CKS2,NUF2,CKS1B,MKI67,TMPO,CENPF,TACC3,FAM64A,SMC4,CCNB2,CKAP2L,CKAP2,AURKB,BUB1,KIF11,ANP32E,TUBB4B,
GTSE1,KIF20B,HJURP,CDCA3,HN1,CDC20,TTK,CDC25C,KIF2C,RANGAP1,NCAPD2,DLGAP5,CDCA2,CDCA8,ECT2,KIF23,HMMR,AURKA,PSRC1,ANLN,LBR,CKAP5,CENPE,CTCF,NEK2,G2E3,GAS2L3,CBX5,CENPA
""".replace("\n",",").split(",") if x.strip()]
        s_genes = cc_genes[:43]
        g2m_genes = cc_genes[43:]

        # Score cell cycle
        sc.tl.score_genes_cell_cycle(filtered, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False)

        # Compute combined S+G2M score
        s_score = filtered.obs.get("S_score", pd.Series(0, index=filtered.obs.index))
        g2m_score = filtered.obs.get("G2M_score", pd.Series(0, index=filtered.obs.index))
        combined_cc_score = s_score + g2m_score
        filtered.obs["combined_cc_score"] = combined_cc_score

        # Filter out proliferating cells
        is_low_prolif = combined_cc_score <= proliferation_threshold
        n_before = filtered.n_obs
        filtered = filtered[is_low_prolif].copy()
        n_after = filtered.n_obs
        n_removed = n_before - n_after
        print(f"  [tonsil] Removed {n_removed:,} proliferating cells (S+G2M > {proliferation_threshold})")
        print(f"  [tonsil] Retained {n_after:,} low-proliferation B cells")

    return filtered


def load_and_filter_tonsil_combined(gc_path: str, mbc_path: str, max_cells: int,
                                    filter_proliferating: bool, proliferation_threshold: float,
                                    exclude_proliferative_gc_types: bool) -> ad.AnnData:
    """Load GC and MBC datasets using backed mode to save memory, filter, concatenate, and subsample."""
    parts = []

    # Helper to load subset from backed file
    def load_subset(path, tag):
        print(f"  [tonsil] Loading {tag} from: {Path(path).name} (backed mode)")
        try:
            # Load in backed mode - only reads metadata initially
            raw = sc.read_h5ad(path, backed="r")
            print(f"    Raw {tag}: {raw.n_obs:,} cells (on disk)")

            # Identify cells to keep based on obs (cell type)
            # We use filter_gc_b_cells but specialized to just return indices/mask first if possible?
            # Actually filter_gc_b_cells expects an AnnData.
            # If we pass backed AnnData, it will read .obs which is fine.
            # But we need to ensure it doesn't try to read .X or copy the whole object.

            # Let's do the filtering manually here to ensure memory safety
            label_col = find_label_column(raw)
            if label_col is None:
                print(f"    (warn) No label column found, loading all...")
                subset = raw.to_memory()
            else:
                labels = raw.obs[label_col].astype(str).str.lower()

                # Re-use the regex logic (simplified/copied for safety or extract to helper if I could)
                # For now I will rely on the pattern matching being fast on just the obs series
                is_target = (
                    labels.str.contains(r"\b(lz|light\s*zone)\b", case=False, na=False) |
                    labels.str.contains(r"\b(dz|dark\s*zone)\b", case=False, na=False) |
                    labels.str.contains(r"\bgcb.*(lz|dz)\b", case=False, na=False) |
                    labels.str.contains(r"\b(lz|dz).*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"lz[-/]dz", case=False, na=False) |
                    labels.str.contains(r"\bcentroblast", case=False, na=False) |
                    labels.str.contains(r"\bcb\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bcentrocyte", case=False, na=False) |
                    labels.str.contains(r"\bcc\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgc\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgerminal\s*center\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgcb\b", case=False, na=False) |
                    labels.str.contains(r"\bgc\s*b\s*lymphocyte", case=False, na=False) |
                    labels.str.contains(r"\bmemory\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bmbc\b", case=False, na=False) |
                    labels.str.contains(r"\bmem\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bb\s*memory", case=False, na=False) |
                    ((labels.str.contains(r"\bactivated\s*b\s*cells?\b", case=False, na=False) |
                      labels.str.contains(r"\bact\s*b\s*cells?\b", case=False, na=False)) &
                     ~labels.str.contains(r"\bplasma\b", case=False, na=False))
                )

                n_match = is_target.sum()
                print(f"    Found {n_match:,} target cells in {tag}")

                if n_match == 0:
                    print("    (warn) No target cells found, skipping...")
                    return None

                # Subset in backed mode (lazy)
                subset_backed = raw[is_target]

                # Now load ONLY the subset into memory
                print(f"    Loading {n_match:,} cells into memory...")
                subset = subset_backed.to_memory()

            # Close file handle
            if hasattr(raw, 'file') and raw.file is not None:
                try: raw.file.close()
                except: pass

            # Post-load processing
            subset = harmonize_human_genes(subset, tag=f"tonsil_{tag}")
            subset = enforce_raw_counts(subset, tag=f"tonsil_{tag}")

            # Apply fine-grained filtering (proliferation, etc.) in memory
            # We reuse filter_gc_b_cells but now it works on a much smaller object
            filt = filter_gc_b_cells(subset,
                                     filter_proliferating=filter_proliferating,
                                     proliferation_threshold=proliferation_threshold,
                                     exclude_proliferative_gc_types=exclude_proliferative_gc_types)
            return filt

        except Exception as e:
            print(f"    (error) Failed to load {tag}: {e}")
            return None

    # 1. Load GC B cells
    gc_filt = load_subset(gc_path, "gc")
    if gc_filt is not None:
        parts.append(gc_filt)
    gc.collect()

    # 2. Load Memory B cells
    mbc_filt = load_subset(mbc_path, "mbc")
    if mbc_filt is not None:
        parts.append(mbc_filt)
    gc.collect()

    if not parts:
        return None

    # 3. Concatenate
    print("  [tonsil] Concatenating datasets...")
    combined = ad.concat(parts, join="outer", index_unique=None, fill_value=0)
    # Restore counts
    combined.layers["counts"] = combined.layers.get("counts", combined.X)

    # Clean up parts to free memory
    del parts, gc_filt, mbc_filt
    gc.collect()

    print(f"  [tonsil] Combined GC + MBC: {combined.n_obs:,} cells")

    # 4. Subsample
    if combined.n_obs > max_cells:
        np.random.seed(0)
        idx = np.random.choice(combined.n_obs, size=max_cells, replace=False)
        combined = combined[idx].copy()
        print(f"  [tonsil] Subsampled to {max_cells:,} cells")

    return combined


# ============================== SCRUBLET =====================================
def run_scrublet_safe(adata: ad.AnnData, sample_name: str) -> pd.Series:
    """Run Scrublet for doublet detection with proper error handling."""
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
        scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            verbose=False
        )

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


# ============================== QC FUNCTIONS =================================
def compute_qc_metrics(adata: ad.AnnData, species: str = "mouse") -> ad.AnnData:
    """Compute QC metrics including mitochondrial gene percentage."""
    if species == "mouse":
        adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    else:  # human
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

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
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= mt_threshold].copy()
    n_after = adata.n_obs
    print(f"    Filtered: {n_before:,} → {n_after:,} cells "
          f"({n_before - n_after:,} removed)")
    return adata


def load_cellbender_h5(h5_path: Path, sample_name: str) -> ad.AnnData:
    """Load CellBender filtered h5 file."""
    print(f"  Loading: {h5_path.name}")

    try:
        adata = sc.read_10x_h5(str(h5_path))
    except Exception as e1:
        print(f"    read_10x_h5 failed: {e1}, trying read_h5ad...")
        adata = sc.read_h5ad(str(h5_path))

    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_name

    if not adata.obs_names.str.contains("-").any():
        adata.obs_names = [f"{bc}-{sample_name}" for bc in adata.obs_names]
    else:
        adata.obs_names = [f"{bc.split('-')[0]}-{sample_name}" for bc in adata.obs_names]

    print(f"    Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


# ============================== HVG FUNCTIONS ================================
def confounder_mask(varnames: pd.Index) -> pd.Series:
    """Create mask for confounder genes."""
    v = pd.Index([str(g) for g in varnames])
    is_mt = v.str.startswith("MT-")
    is_ribo = v.str.startswith(("RPS", "RPL"))
    is_ig = v.str.match(r"^(IGH|IGK|IGL)[VDJC].*", na=False)
    tcr_prefixes = ("TRAV", "TRBV", "TRGV", "TRDV", "TRAJ", "TRBJ",
                    "TRGJ", "TRDJ", "TRAC", "TRBC", "TRGC", "TRDC")
    is_tcr = v.str.startswith(tcr_prefixes)
    return is_mt | is_ribo | is_ig | is_tcr


def build_sample_batch(adx, tag=""):
    """Build sample_batch column for HVG detection."""
    for cand in ["sample_id", "donor", "library_id", "orig.ident", "sample",
                 "Sample", "Library", "Donor", "batch"]:
        if cand in adx.obs:
            adx.obs["sample_batch"] = adx.obs[cand].astype(str)
            print(f"  [{tag}] sample_batch ← '{cand}'")
            return adx
    adx.obs["sample_batch"] = (adx.obs.get("study", "study").astype(str) + "_" +
                               adx.obs.get("sample_type", "sample").astype(str))
    print(f"  [{tag}] sample_batch ← fallback")
    return adx


def compute_shared_hvgs(adata_list, n_top=3500, min_datasets=2):
    """Compute HVGs shared across datasets."""
    hvg_counts = {}
    for i, adx in enumerate(adata_list):
        print(f"    HVGs for dataset {i+1}/{len(adata_list)} (seurat_v3 on counts)")
        bk = "sample_batch" if "sample_batch" in adx.obs.columns else None
        sc.pp.highly_variable_genes(
            adx, layer="counts", n_top_genes=n_top, flavor="seurat_v3",
            subset=False, batch_key=bk
        )
        for g in adx.var_names[adx.var["highly_variable"]]:
            hvg_counts[g] = hvg_counts.get(g, 0) + 1
    min_count = max(min_datasets, len(adata_list) // 2)
    shared = [g for g, c in hvg_counts.items() if c >= min_count]
    print(f"  Genes HVG in ≥{min_count} datasets: {len(shared):,}")
    return set(shared)


# ============================== CT2 FUNCTIONS ================================
def stream_counts_to_ct2_txt(adata: ad.AnnData, dest: Path) -> None:
    """Write CytoTRACE2-compatible matrix."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    if "counts" not in adata.layers:
        raise RuntimeError("[ct2-io] layers['counts'] missing.")

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
            g = str(adata.var_names[j]).strip().upper()
            if g in na_like:
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

        id_col = next((c for c in ("cell", "cell_id", "barcode") if c in df.columns), df.columns[0])
        score_cols = [c for c in df.columns if ("cytotrace2" in c and "score" in c)]

        if not score_cols:
            continue

        tmp = pd.DataFrame(index=df[id_col].astype(str).values)
        tmp["cytotrace2_score"] = _coerce_float_series(df[score_cols[0]]).values

        pot_cols = [c for c in df.columns if "potency" in c]
        if pot_cols:
            tmp["cytotrace2_potency"] = df[pot_cols[0]].astype(str).values

        n_match = len(set(tmp.index) & set(map(str, adata.obs_names)))
        if best is None or n_match > best[0]:
            best = (n_match, f, tmp)

    if best is None:
        raise FileNotFoundError(f"[ct2] No results found under {outdir}")

    print(f"[ct2] Using: {best[1]} (matched {best[0]:,} cells)")
    df_best = best[2]
    df_best = df_best[~df_best.index.duplicated(keep="first")]
    return df_best.reindex(adata.obs_names.astype(str))


# ============================== PLOTTING FUNCTIONS ===========================
def _save_umap(adata, color, fname, title=None, palette=None, cmap=None):
    """Save UMAP plot in PNG, PDF, and SVG formats."""
    fig, ax = plt.subplots(figsize=(10, 9))
    sc.pl.umap(adata, color=color, title=(title or color), palette=palette, cmap=cmap,
               frameon=False, legend_loc="right margin", ax=ax, show=False, s=25)
    plt.tight_layout()

    # Save in multiple formats
    base_name = fname.replace('.png', '')
    plt.savefig(FIGDIR / f"{base_name}.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIGDIR / f"{base_name}.pdf", bbox_inches="tight")
    plt.savefig(FIGDIR / f"{base_name}.svg", bbox_inches="tight")
    plt.close()
    print(f"  ✓ {base_name} (.png, .pdf, .svg)")


def save_umaps_by_sample(adata: ad.AnnData, figdir: Path, group_key: str = "sample_batch"):
    """Save individual UMAPs highlighting each sample in PNG, PDF, and SVG formats."""
    if group_key not in adata.obs:
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
            adata, color="_highlight", palette=palette, frameon=False,
            legend_loc=None, show=False, s=25, ax=ax, title=f"{group_key} → {grp}"
        )
        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.pdf", bbox_inches="tight")
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.svg", bbox_inches="tight")
        plt.close(fig)
        adata.obs.drop(columns="_highlight", inplace=True)

    print(f"  ✓ Per-group UMAPs saved to {outdir} (.png, .pdf, .svg)")


def save_individual_sample_umaps(adata: ad.AnnData, figdir: Path,
                                  group_key: str = "sample_batch"):
    """Save UMAPs showing ONLY cells from each sample in PNG, PDF, and SVG formats."""
    if group_key not in adata.obs:
        return

    outdir = figdir / f"umap_individual_{group_key}"
    outdir.mkdir(parents=True, exist_ok=True)

    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())

    # Get global UMAP limits
    umap_coords = adata.obsm["X_umap"]
    x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
    y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    for grp in unique_groups:
        mask = labels == grp
        adata_sub = adata[mask].copy()
        n_cells = adata_sub.n_obs

        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)

        fig, ax = plt.subplots(figsize=(9, 8))
        ax.scatter(
            adata_sub.obsm["X_umap"][:, 0],
            adata_sub.obsm["X_umap"][:, 1],
            c="#d62728", s=15, alpha=0.7, rasterized=True
        )
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_title(f"{grp} (n={n_cells:,})", fontsize=14)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.pdf", bbox_inches="tight")
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.svg", bbox_inches="tight")
        plt.close(fig)

    print(f"  ✓ Individual sample UMAPs saved to {outdir} (.png, .pdf, .svg)")


def add_ct2_umap_and_plots(adata: ad.AnnData, figdir: Path, title_suffix=""):
    """Add CytoTRACE2 UMAP and histogram plots in PNG, PDF, and SVG formats."""
    figdir.mkdir(parents=True, exist_ok=True)

    if "cytotrace2_score" in adata.obs:
        fig, ax = plt.subplots(figsize=(10, 9))
        sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=30, frameon=False,
                   cmap="viridis", title=f"CytoTRACE2 score {title_suffix}")
        fig.savefig(figdir / "umap_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        fig.savefig(figdir / "umap_cytotrace2_score.pdf", bbox_inches="tight")
        fig.savefig(figdir / "umap_cytotrace2_score.svg", bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        adata.obs["cytotrace2_score"].astype(float).plot(kind="hist", bins=50, ax=ax)
        ax.set_xlabel("CytoTRACE2 score")
        ax.set_ylabel("Cell count")
        ax.set_title("Distribution of CT2 scores")
        fig.savefig(figdir / "hist_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        fig.savefig(figdir / "hist_cytotrace2_score.pdf", bbox_inches="tight")
        fig.savefig(figdir / "hist_cytotrace2_score.svg", bbox_inches="tight")
        plt.close(fig)


# ============================== MAIN PIPELINE ================================
print("=" * 84)
print("STEP 0 — Ortholog mapping (BioMart preferred)")
print("=" * 84)
M2H_1TO1 = get_m2h_map(prefer_biomart=True)

# ============================== LOAD MOUSE MALIGNANT =========================
print("\n" + "=" * 84)
print("STEP 1 — Load and QC mouse malignant samples (CellBender filtered)")
print("=" * 84)

mouse_list = []
qc_stats = []

for sample_name, sample_info in MOUSE_MALIGNANT_FILES.items():
    h5_path = sample_info["path"]

    if not h5_path.exists():
        print(f"  WARNING: {h5_path} not found, skipping...")
        continue

    # Load sample
    adata = load_cellbender_h5(h5_path, sample_name)

    # Add metadata
    adata.obs["condition"] = sample_info["condition"]
    adata.obs["replicate"] = sample_info["replicate"]
    adata.obs["species"] = "mouse"
    adata.obs["study"] = "mouse_malignant"
    adata.obs["disease_state"] = f"Mouse_{sample_info['condition']}"

    # Compute QC metrics
    adata = compute_qc_metrics(adata, species="mouse")
    n_pre = adata.n_obs

    # Run Scrublet
    doublet_mask = run_scrublet_safe(adata, sample_name)
    adata.obs["predicted_doublet"] = doublet_mask.values

    # Filter cells
    adata = filter_cells(adata, mt_threshold=MT_THRESHOLD, min_genes=MIN_GENES)

    # Remove doublets
    n_pre_doublet = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_doublets_removed = n_pre_doublet - adata.n_obs
    print(f"    Removed {n_doublets_removed:,} doublets → {adata.n_obs:,} cells")

    # Store counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    qc_stats.append({
        "sample": sample_name,
        "condition": sample_info["condition"],
        "species": "mouse",
        "cells_raw": n_pre,
        "cells_after_qc": adata.n_obs
    })

    mouse_list.append(adata)
    print(f"    Final: {adata.n_obs:,} cells\n")

# Concatenate mouse samples
mouse_adata = ad.concat(mouse_list, join="outer", index_unique=None)
print(f"  Mouse combined: {mouse_adata.n_obs:,} cells × {mouse_adata.n_vars:,} genes")

# Ensure counts layer
if "counts" not in mouse_adata.layers:
    mouse_adata.layers["counts"] = mouse_adata.X.copy()

# Convert mouse genes to human
mouse_adata = convert_mouse_to_human_genes(mouse_adata, M2H_1TO1, tag="mouse")
mouse_adata = enforce_raw_counts(mouse_adata, tag="mouse")

# ============================== LOAD HUMAN DLBCL =============================
print("\n" + "=" * 84)
print("STEP 2 — Load human DLBCL samples")
print("=" * 84)

dlbcl_parts = []
for p in DLBCL_FILES:
    if not Path(p).exists():
        print(f"  WARNING: {p} not found, skipping...")
        continue

    a = sc.read_h5ad(p)
    nm = Path(p).stem
    a.obs["species"] = "human"
    a.obs["study"] = ("DLBCL_Alizadeh" if "Alizadeh" in p else "DLBCL_Roider")
    a.obs["sample_type"] = "dlbcl"
    a.obs["disease_state"] = "DLBCL"
    a = harmonize_human_genes(a, tag=nm)
    a = enforce_raw_counts(a, tag=nm)
    if "sample" not in a.obs:
        a.obs["sample"] = nm
    dlbcl_parts.append(a)
    print(f"  + {Path(p).name}: {a.n_obs:,} cells, {a.n_vars:,} genes")

    qc_stats.append({
        "sample": nm,
        "condition": "DLBCL",
        "species": "human",
        "cells_raw": a.n_obs,
        "cells_after_qc": a.n_obs
    })

dlbcl = ad.concat(dlbcl_parts, join="outer", index_unique=None)
dlbcl.layers["counts"] = dlbcl.layers.get("counts", dlbcl.X)
print(f"  DLBCL combined: {dlbcl.n_obs:,} × {dlbcl.n_vars:,}")

# ============================== LOAD HUMAN TONSIL ============================
print("\n" + "=" * 84)
print("STEP 2.5 — Load Human Tonsil GC B cells")
print("=" * 84)

# Load and filter Tonsil data (GC + MBC combined)
tonsil = load_and_filter_tonsil_combined(
    gc_path=TONSIL_GC_PATH,
    mbc_path=TONSIL_MBC_PATH,
    max_cells=MAX_TONSIL_CELLS,
    filter_proliferating=FILTER_PROLIFERATING_TONSIL,
    proliferation_threshold=PROLIFERATION_THRESHOLD,
    exclude_proliferative_gc_types=EXCLUDE_PROLIFERATIVE_GC_TYPES
)

if tonsil is None:
    raise RuntimeError("Failed to load tonsil data.")

tonsil.obs["species"] = "human"
tonsil.obs["study"] = "Tonsil_GCBC_MBC"
tonsil.obs["sample_type"] = "tonsil_b_cells"
tonsil.obs["disease_state"] = "Tonsil_Normal"
print(f"  Tonsil (GC+MBC) final: {tonsil.n_obs:,} × {tonsil.n_vars:,}")

qc_stats.append({
    "sample": "Tonsil_GC_MBC",
    "condition": "Tonsil_Normal",
    "species": "human",
    "cells_raw": "N/A (combined)",
    "cells_after_qc": tonsil.n_obs
})

# Save QC stats
qc_df = pd.DataFrame(qc_stats)
qc_df.to_csv(OUTDIR / "qc_statistics.csv", index=False)
print(f"\n  QC stats saved to: {OUTDIR / 'qc_statistics.csv'}")
print(qc_df.to_string())

# ============================== NORMALIZE ====================================
print("\n" + "=" * 84)
print("STEP 3 — Normalize datasets")
print("=" * 84)

mouse_adata = harmonize_preprocessing(mouse_adata, tag="mouse")
dlbcl = harmonize_preprocessing(dlbcl, tag="DLBCL")
tonsil = harmonize_preprocessing(tonsil, tag="tonsil")

# ============================== COMMON GENES + HVGs ==========================
print("\n" + "=" * 84)
print("STEP 4 — Common genes + HVGs")
print("=" * 84)

datasets = [mouse_adata, dlbcl, tonsil]
common_genes = sorted(list(set(datasets[0].var_names).intersection(
    *[set(x.var_names) for x in datasets[1:]])))
print(f"  Common genes across datasets: {len(common_genes):,}")

# Build sample_batch
mouse_adata = build_sample_batch(mouse_adata, "mouse")
dlbcl = build_sample_batch(dlbcl, "DLBCL")
tonsil = build_sample_batch(tonsil, "tonsil")

adata_for_hvg = [x[:, common_genes].copy() for x in datasets]
hvgs = compute_shared_hvgs(adata_for_hvg, n_top=4000, min_datasets=2)
hvgi = pd.Index(sorted(hvgs))
hvgs_clean = sorted(list(hvgi[~confounder_mask(hvgi)]))
print(f"  HVGs after removing confounders: {len(hvgs_clean):,}")

# ============================== CONCATENATE ==================================
print("\n" + "=" * 84)
print("STEP 5 — Concatenate datasets")
print("=" * 84)

# For scVI: use common genes (inner join)
adata_full = ad.concat([x[:, common_genes] for x in datasets],
                        merge="same", join="inner", index_unique=None)
if "counts" not in adata_full.layers:
    raise RuntimeError("layers['counts'] missing in adata_full.")
print(f"  Combined FULL (common genes): {adata_full.n_obs:,} × {adata_full.n_vars:,}")

# For CT2: use UNION of all genes (outer join) - CT2 explicitly asks NOT to pre-filter genes
if RUN_CYTOTRACE2:
    print("  Creating CT2 object with UNION of all genes (outer join, fill missing with zeros)...")
    all_genes_union = sorted(list(set().union(*[set(x.var_names) for x in datasets])))
    print(f"  Union of all genes: {len(all_genes_union):,} (vs {len(common_genes):,} common genes)")

    # Verify all datasets have counts layer before concatenation
    for i, ds in enumerate(datasets):
        if "counts" not in ds.layers:
            raise RuntimeError(f"[ct2] Dataset {i} missing layers['counts'].")

    # Concatenate with outer join - fills missing genes with zeros
    adata_ct2 = ad.concat(datasets, merge="same", join="outer", index_unique=None, fill_value=0)

    # Verify/reconstruct counts layer
    if "counts" not in adata_ct2.layers:
        if _is_intlike(adata_ct2.X):
            adata_ct2.layers["counts"] = adata_ct2.X.copy()
            print("  [ct2] Reconstructed counts layer from X")
        else:
            raise RuntimeError("[ct2] Could not find or reconstruct raw counts for CT2.")
    print(f"  CT2 object: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes (UNION, no pre-filtering)")
else:
    print("  Skipping CT2 object creation (RUN_CYTOTRACE2=False)")
    adata_ct2 = None

# Cross-species HVG filtering
if "species" in adata_full.obs and adata_full.obs["species"].nunique() > 1:
    print("  Ensuring HVGs detected across species (≥5% of cells each)")
    hvg_view = adata_full[:, hvgs_clean]
    counts_layer = hvg_view.layers["counts"] if "counts" in hvg_view.layers else hvg_view.X
    is_sparse = sp.issparse(counts_layer)
    if is_sparse:
        counts_layer = counts_layer.tocsr()
    species_labels = hvg_view.obs["species"].astype(str).values
    keep_mask = np.ones(len(hvgs_clean), dtype=bool)
    for sp_name in np.unique(species_labels):
        row_idx = np.where(species_labels == sp_name)[0]
        if row_idx.size == 0:
            continue
        if is_sparse:
            sub = counts_layer[row_idx, :]
            detected = np.asarray((sub > 0).mean(axis=0)).ravel()
        else:
            sub = counts_layer[row_idx, :]
            detected = (sub > 0).astype(float).mean(axis=0)
        keep_mask &= detected >= 0.05
    hvgs_filtered = [g for g, keep in zip(hvgs_clean, keep_mask) if keep]
    if hvgs_filtered:
        hvgs_clean = hvgs_filtered
        print(f"  HVGs after cross-species detection ≥5%: {len(hvgs_clean):,}")

adata_all = adata_full[:, hvgs_clean].copy()
print(f"  Combined HVG: {adata_all.n_obs:,} × {adata_all.n_vars:,}")

# ============================== scVI =========================================
print("\n" + "=" * 84)
print("STEP 6 — scVI integration (NO cell cycle regression)")
print("=" * 84)

for cat in ("sample_batch", "species"):
    if cat in adata_all.obs:
        adata_all.obs[cat] = adata_all.obs[cat].astype("category")

scvi.model.SCVI.setup_anndata(
    adata_all, layer="counts", batch_key="sample_batch",
    categorical_covariate_keys=["species"]
)

model = scvi.model.SCVI(
    adata_all,
    n_latent=96,
    n_layers=2,
    dropout_rate=0.15,
    gene_likelihood="nb",
    dispersion="gene-batch",
    use_layer_norm="both",
    use_batch_norm="none"
)

max_epochs = 150 if torch.cuda.is_available() else 600
print(f"  Training scVI for {max_epochs} epochs (early stopping patience=60)...")
model.train(
    max_epochs=max_epochs,
    plan_kwargs={"lr": 3e-4, "n_epochs_kl_warmup": 200, "reduce_lr_on_plateau": True},
    check_val_every_n_epoch=20,
    early_stopping=True,
    early_stopping_patience=60
)

adata_all.obsm["X_scvi"] = model.get_latent_representation()

# ============================== UMAP + CLUSTERING ============================
print("\n" + "=" * 84)
print("STEP 7 — UMAP and clustering")
print("=" * 84)

sc.pp.neighbors(adata_all, use_rep="X_scvi", n_neighbors=30)
sc.tl.umap(adata_all, min_dist=0.2, spread=1.5)

for res in [0.5, 1.0, 1.5]:
    sc.tl.leiden(adata_all, resolution=res, key_added=f"leiden_{res}")
    print(f"  Leiden res={res}: {adata_all.obs[f'leiden_{res}'].nunique()} clusters")

# ============================== CELL CYCLE + GENE SET SCORES =================
print("\n" + "=" * 84)
print("STEP 8 — Cell cycle scoring + gene-set scores")
print("=" * 84)

# Cell cycle scoring (for visualization only, NOT regressed out)
cc_genes = [x.strip() for x in """
MCM5,PCNA,TYMS,FEN1,MCM2,MCM4,RRM1,UNG,GINS2,MCM6,CDCA7,DTL,PRIM1,UHRF1,MLF1IP,HELLS,RFC2,RPA2,NASP,RAD51AP1,GMNN,WDR76,SLBP,CCNE2,UBR7,POLD3,MSH2,
ATAD2,RAD51,RRM2,CDC45,CDC6,EXO1,TIPIN,DSCC1,BLM,CASP8AP2,USP1,CLSPN,POLA1,CHAF1B,BRIP1,E2F8,
HMGB2,CDK1,NUSAP1,UBE2C,BIRC5,TPX2,TOP2A,NDC80,CKS2,NUF2,CKS1B,MKI67,TMPO,CENPF,TACC3,FAM64A,SMC4,CCNB2,CKAP2L,CKAP2,AURKB,BUB1,KIF11,ANP32E,TUBB4B,
GTSE1,KIF20B,HJURP,CDCA3,HN1,CDC20,TTK,CDC25C,KIF2C,RANGAP1,NCAPD2,DLGAP5,CDCA2,CDCA8,ECT2,KIF23,HMMR,AURKA,PSRC1,ANLN,LBR,CKAP5,CENPE,CTCF,NEK2,G2E3,GAS2L3,CBX5,CENPA
""".replace("\n", ",").split(",") if x.strip()]
s_genes = cc_genes[:43]
g2m_genes = cc_genes[43:]

# Filter to genes present
s_genes_present = [g for g in s_genes if g in adata_all.var_names]
g2m_genes_present = [g for g in g2m_genes if g in adata_all.var_names]

if s_genes_present and g2m_genes_present:
    sc.tl.score_genes_cell_cycle(adata_all, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
    print(f"  Cell cycle scoring: {len(s_genes_present)} S genes, {len(g2m_genes_present)} G2M genes")
    print(f"  Phase distribution:")
    print(adata_all.obs['phase'].value_counts().to_string())
else:
    print(f"  WARNING: Not enough cell cycle genes found for scoring")

# Gene set scores
avail = set(adata_all.var_names)
ox = [g for g in OXPHOS_GENES if g in avail]
bcr = [g for g in BCR_GENES if g in avail]
if ox:
    sc.tl.score_genes(adata_all, gene_list=ox, score_name='oxphos_score')
    print(f"  OXPHOS score: {len(ox)} genes")
if bcr:
    sc.tl.score_genes(adata_all, gene_list=bcr, score_name='bcr_score')
    print(f"  BCR score: {len(bcr)} genes")

# ============================== CytoTRACE2 ===================================
print("\n" + "=" * 84)
print("STEP 9 — CytoTRACE2")
print("=" * 84)

if RUN_CYTOTRACE2:
    # Use adata_ct2 (union of all genes) for CT2, as recommended by CT2 authors
    print("  Using CT2 object with UNION of all genes (no pre-filtering, as recommended by CT2)")
    print(f"  CT2 input size: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes")

    # Filter to genes expressed in at least some cells to reduce file size/memory
    if adata_ct2.n_vars > 30000:
        print(f"  [ct2] Filtering to genes expressed in ≥0.1% of cells to reduce memory...")
        X_counts = adata_ct2.layers["counts"] if "counts" in adata_ct2.layers else adata_ct2.X
        if sp.issparse(X_counts):
            n_expressed = np.array((X_counts > 0).sum(axis=0)).ravel()
        else:
            n_expressed = (X_counts > 0).sum(axis=0)
        min_cells = max(1, int(adata_ct2.n_obs * 0.001))  # At least 0.1% of cells
        keep_genes = n_expressed >= min_cells
        n_keep = keep_genes.sum()
        print(f"  [ct2] Keeping {n_keep:,} genes (expressed in ≥{min_cells} cells) out of {adata_ct2.n_vars:,}")
        adata_ct2 = adata_ct2[:, keep_genes].copy()
        print(f"  [ct2] Filtered CT2 input: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes")

    stream_counts_to_ct2_txt(adata_ct2, CT2_INPUT_TXT)

    ct2_obj = run_ct2_python(CT2_INPUT_TXT, species="human", outdir=CT2_OUTDIR)

    if ct2_obj is None:
        cmd = f"cytotrace2 -f {shlex.quote(str(CT2_INPUT_TXT))} -sp human --output-dir {shlex.quote(str(CT2_OUTDIR))} --disable-plotting"
        print("[ct2] CLI:", cmd)
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode != 0:
            print(f"[ct2] WARNING: CytoTRACE2 CLI failed")

    # Parse CT2 scores
    df_scores = None
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
            idx = adata_ct2.obs_names.astype(str)

            df_scores = pd.DataFrame(index=idx)
            sname = next((c for c in obs.columns if ("cytotrace2" in c and "score" in c)), None)
            if sname:
                df_scores["cytotrace2_score"] = _coerce_float_series(obs.reindex(idx)[sname]).values

            p = next((c for c in obs.columns if "potency" in c), None)
            if p:
                df_scores["cytotrace2_potency"] = obs.reindex(idx)[p].astype(str).values
        except Exception as e:
            print(f"[ct2] Could not parse object: {e}")

    if df_scores is None:
        try:
            df_scores = parse_ct2_scores(CT2_OUTDIR, adata_ct2)
        except FileNotFoundError:
            print("[ct2] WARNING: CT2 results not found")

    if df_scores is not None and not df_scores.empty:
        for col in df_scores.columns:
            ser_full = df_scores[col].reindex(adata_full.obs_names.astype(str))
            ser_all = df_scores[col].reindex(adata_all.obs_names.astype(str))
            adata_full.obs[col] = ser_full
            adata_all.obs[col] = ser_all
        print("✓ CytoTRACE2 scores attached")
else:
    print("  Skipping CytoTRACE2 execution (RUN_CYTOTRACE2=False)")


# ============================== FIGURES ======================================
print("\n" + "=" * 84)
print("STEP 10 — Generate figures")
print("=" * 84)

disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000",
    "Tonsil_GC_B": "#4169E1",
    "Tonsil_Normal": "#4169E1"  # Same color as Tonsil_GC_B (royal blue)
}
species_colors = {"mouse": "#98FB98", "human": "#6495ED"}
phase_colors = {"G1": "#3498db", "S": "#e74c3c", "G2M": "#2ecc71"}

_save_umap(adata_all, "disease_state", "umap_disease_state.png", "Disease State", disease_colors)
_save_umap(adata_all, "species", "umap_species.png", "Species", species_colors)
_save_umap(adata_all, "sample_batch", "umap_sample_batch.png", "Sample Batch")
_save_umap(adata_all, "leiden_1.0", "umap_leiden_1.0.png", "Leiden (res=1.0)")
_save_umap(adata_all, "phase", "umap_cell_cycle.png", "Cell Cycle Phase", phase_colors)

if "oxphos_score" in adata_all.obs:
    _save_umap(adata_all, "oxphos_score", "umap_oxphos_score.png", "OXPHOS Score", cmap="RdYlBu_r")
if "bcr_score" in adata_all.obs:
    _save_umap(adata_all, "bcr_score", "umap_bcr_score.png", "BCR Score", cmap="RdYlBu_r")

save_umaps_by_sample(adata_all, FIGDIR, group_key="sample_batch")
save_individual_sample_umaps(adata_all, FIGDIR, group_key="sample_batch")

if "cytotrace2_score" in adata_all.obs:
    add_ct2_umap_and_plots(adata_all, FIGDIR_CT2, title_suffix="(integrated)")
    # Also save CT2 on the main UMAP
    _save_umap(adata_all, "cytotrace2_score", "umap_cytotrace2_score.png",
               "CytoTRACE2 Score (scVI UMAP)", cmap="viridis")
elif RUN_CYTOTRACE2:
    print("  (warn) CytoTRACE2 scores missing despite RUN_CYTOTRACE2=True")
else:
    print("  Skipping CytoTRACE2 plots")

# ============================== SAVE =========================================
print("\n" + "=" * 84)
print("STEP 11 — Save outputs")
print("=" * 84)

try:
    model.save(OUTDIR / "scvi_model", overwrite=True)
    print(f"  ✓ scVI model → {OUTDIR / 'scvi_model'}")
except Exception as e:
    print(f"  (warn) Could not save scVI model: {e}")

adata_all.write(OUTDIR / "integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")
print(f"  ✓ AnnData → {OUTDIR / 'integrated_human_dlbcl_mouse_malignant_tonsil.h5ad'}")

adata_all.obs.to_csv(OUTDIR / "cell_metadata.csv")
print(f"  ✓ Metadata → {OUTDIR / 'cell_metadata.csv'}")

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTDIR}")
print(f"  Figures: {FIGDIR}")
print(f"  CT2 results: {CT2_OUTDIR}")
print(f"  Total cells: {adata_all.n_obs:,}")
print(f"  Total genes: {adata_all.n_vars:,}")

# Summary by disease state
print("\n  Cells by disease state:")
for ds, count in adata_all.obs["disease_state"].value_counts().items():
    print(f"    {ds}: {count:,}")

print("\nDONE.\n")


os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

import gc
import re
import shlex
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import urllib.request
import urllib.parse

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
OUTDIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration")
FIGDIR = OUTDIR / "figures"
FIGDIR_CT2 = OUTDIR / "figures_ct2"
CT2_WORKDIR = OUTDIR / "ct2_io"
CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"
CT2_OUTDIR = OUTDIR / "cytotrace2_results"

for p in (OUTDIR, FIGDIR, FIGDIR_CT2, CT2_WORKDIR, CT2_OUTDIR):
    p.mkdir(parents=True, exist_ok=True)

# Human DLBCL files
DLBCL_FILES = [
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL1_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL2_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/DLBCL3_raw.h5ad",
    "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/DLBCL/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad",
]

# Human Tonsil data files
TONSIL_GC_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_GCBC_RNA.h5ad"
TONSIL_MBC_PATH = "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad"
MAX_TONSIL_CELLS = 12500  # Maximum total GC + Memory B cells to keep

# Proliferation filtering for tonsil cells
FILTER_PROLIFERATING_TONSIL = False  # Changed to False to preserve Dark Zone (Centroblasts)
PROLIFERATION_THRESHOLD = 0.20  # Relaxed threshold if enabled (was 0.10)
EXCLUDE_PROLIFERATIVE_GC_TYPES = False  # Changed to False to preserve Dark Zone annotations

# Toggle CytoTRACE2 (Set to True later to run CT2 after scVI completes)
RUN_CYTOTRACE2 = False  # Set to False to skip CT2 and save memory

# Mouse malignant samples (CellBender filtered)
CELLBENDER_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/cellbender_filtered")
MOUSE_MALIGNANT_FILES = {
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
    "SIGAD5_Malignant_R2": {
        "path": CELLBENDER_DIR / "SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R2"
    },
    "SIGAH1_Malignant_R1": {
        "path": CELLBENDER_DIR / "SIGAH1_Malignant_R1_GEX_cellbender_filtered.h5",
        "condition": "Malignant",
        "replicate": "R1"
    },
}

# QC thresholds
MT_THRESHOLD = 10.0
MIN_GENES = 200
MIN_CELLS = 3

print("=" * 84)
print("HUMAN DLBCL + MALIGNANT MOUSE + TONSIL GC B INTEGRATION (scVI + CytoTRACE2)")
print("=" * 84)
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"Output dir    : {OUTDIR}")
print(f"MT threshold  : {MT_THRESHOLD}%")
print(f"HVGs for integration: 4000")
print(f"Tonsil B cell types: LZ/DZ + centroblasts + centrocytes + GC B cells + memory B cells")
print(f"Max tonsil cells: {MAX_TONSIL_CELLS:,}")
if FILTER_PROLIFERATING_TONSIL:
    print(f"Proliferation filtering: ON (S+G2M score threshold: {PROLIFERATION_THRESHOLD})")
    print(f"  -> WARNING: This may remove Dark Zone (Centroblast) populations.")
else:
    print(f"Proliferation filtering: OFF (Preserving Dark Zone/Centroblasts)")
if EXCLUDE_PROLIFERATIVE_GC_TYPES:
    print(f"Exclude proliferative GC types: ON")
else:
    print(f"Exclude proliferative GC types: OFF")
print()

# ============================== GENE SETS ====================================
OXPHOS_GENES = [
    'COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A2','COX7B','COX7C','COX8A',
    'CYC1','CYCS','NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9',
    'NDUFA10','NDUFA11','NDUFA12','NDUFA13','NDUFAB1','NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5',
    'NDUFB6','NDUFB7','NDUFB8','NDUFB9','NDUFB10','NDUFB11','NDUFC1','NDUFC2','NDUFS1','NDUFS2',
    'NDUFS3','NDUFS4','NDUFS5','NDUFS6','NDUFS7','NDUFS8','NDUFV1','NDUFV2','NDUFV3','SDHA','SDHB',
    'SDHC','SDHD','UQCR10','UQCR11','UQCRB','UQCRC1','UQCRC2','UQCRFS1','UQCRH','UQCRQ',
    'ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E','ATP5MC1','ATP5MC2','ATP5MC3','ATP5ME','ATP5MF',
    'ATP5MG','ATP5PB','ATP5PD','ATP5PF','ATP5PO'
]
BCR_GENES = [
    'CD79A','CD79B','CD19','CD22','CD72','CR2','FCRL1','FCRL2','FCRL3','FCRL4','FCRL5','MS4A1',
    'IGHM','IGHD','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHE','BTK','LYN','SYK','BLK',
    'BLNK','PIK3CD','PIK3AP1','PLCG2','PRKCB','NFKB1','NFKB2','REL','RELA','NFATC1','NFATC2',
    'BCL10','CARD11','MALT1','MAP3K7','IKBKB','IKBKG','CHUK','PTPN6','PTPRC','VAV1','VAV2',
    'VAV3','GRB2','SOS1','SOS2','HRAS','KRAS','NRAS','RAF1','MAP2K1','MAP2K2','MAPK1','MAPK3'
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


def enforce_raw_counts(adata: ad.AnnData, prefer_layers=("counts", "raw_counts"), tag=""):
    """Ensure layers['counts'] contains raw integer counts."""
    for lyr in prefer_layers:
        if lyr in adata.layers and _is_intlike(adata.layers[lyr]):
            adata.layers["counts"] = adata.layers[lyr]
            print(f"  [{tag}] Using '{lyr}' as counts.")
            return adata
    if adata.raw is not None and _is_intlike(adata.raw.X):
        adata.layers["counts"] = adata.raw.X
        print(f"  [{tag}] Copied raw.X as counts.")
        return adata
    if _is_intlike(adata.X):
        adata.layers["counts"] = adata.X.copy()
        print(f"  [{tag}] Using X as counts.")
        return adata
    raise ValueError(f"[{tag}] No integer-like counts found for scVI.")


def harmonize_preprocessing(adata, target_sum=1e4, tag=""):
    """Normalize and log-transform, preserving counts."""
    print(f"  [{tag}] normalize_total→log1p (counts preserved)")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["normalized"] = adata.X.copy()
    return adata


# ============================== ORTHOLOG MAPPING =============================
def fetch_biomart_m2h_one2one(host: str = "https://www.ensembl.org") -> Dict[str, str]:
    """BioMart mouse→human 1:1 via martservice XML."""
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="mmusculus_gene_ensembl" interface="default">
    <Filter name="with_hsapiens_homolog" excluded="0"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
    <Attribute name="hsapiens_homolog_orthology_type"/>
  </Dataset>
</Query>"""
    url = f"{host}/biomart/martservice?query=" + urllib.parse.quote(xml)
    print("  [ortholog] Querying Ensembl BioMart for 1:1 orthologs ...")
    txt = urllib.request.urlopen(url, timeout=120).read().decode("utf-8")
    m2h = {}
    for line in txt.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        mm, hs, typ = (p.strip() for p in parts)
        if typ == "ortholog_one2one" and mm and hs and mm not in m2h:
            m2h[mm] = hs.upper()
    if len(m2h) < 1000:
        raise RuntimeError("BioMart map unexpectedly small")
    print(f"  [ortholog] BioMart 1:1 pairs: {len(m2h):,}")
    return m2h


def load_one_to_one_orthologs_mgi(mgi_file: Optional[Path]) -> Dict[str, str]:
    """Load MGI ortholog mapping as fallback."""
    if mgi_file is None or not mgi_file.exists():
        return {}
    df = pd.read_csv(mgi_file, sep="\t", dtype=str)
    out: Dict[str, str] = {}
    for _, g in df.groupby("DB Class Key"):
        m = g[g["Common Organism Name"] == "mouse, laboratory"]["Symbol"].dropna().unique()
        h = g[g["Common Organism Name"] == "human"]["Symbol"].dropna().unique()
        if len(m) == 1 and len(h) == 1:
            out[m[0]] = h[0].upper()
    print(f"  [ortholog] MGI 1:1 pairs: {len(out):,}")
    return out


def get_m2h_map(prefer_biomart=True) -> Dict[str, str]:
    """Get mouse-to-human ortholog mapping."""
    m2h = {}
    if prefer_biomart:
        try:
            m2h.update(fetch_biomart_m2h_one2one())
        except Exception as e:
            print(f"  [ortholog] BioMart failed: {e}")
    rpt = OUTDIR / "HOM_MouseHumanSequence.rpt"
    if not rpt.exists():
        try:
            print("  [ortholog] Downloading MGI HOM report ...")
            urllib.request.urlretrieve(
                "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                str(rpt)
            )
        except Exception as e:
            print(f"  [ortholog] MGI download failed: {e}")
            return m2h
    mgi_map = load_one_to_one_orthologs_mgi(rpt)
    for k, v in mgi_map.items():
        m2h.setdefault(k, v)
    print(f"  [ortholog] total unique pairs: {len(m2h):,}")
    return m2h


def convert_mouse_to_human_genes(adata: ad.AnnData, m2h: Dict[str, str], tag=""):
    """Convert mouse gene symbols to human orthologs."""
    print(f"  [{tag}] mouse→human gene symbols (1:1 map → uppercase fallback)")
    mapped = [m2h.get(str(g), str(g).upper()) for g in adata.var_names]
    adata.var_names = pd.Index(mapped)
    adata.var_names_make_unique()
    return adata


def harmonize_human_genes(adata: ad.AnnData, tag=""):
    """Harmonize human gene symbols to uppercase."""
    adata.var_names = pd.Index([str(g).upper() for g in adata.var_names])
    adata.var_names_make_unique()
    return adata


# ============================== TONSIL FILTERING =============================
def find_label_column(adata: ad.AnnData) -> Optional[str]:
    """Find cell type annotation column in tonsil atlas."""
    candidates = [
        "annotation_20230508", "annotation_20220414", "cell_type", "CellType",
        "celltype", "label", "celltype_l2", "celltype_l1", "Azimuth.celltype",
        "celltype_final", "predicted.celltype", "predicted_celltype"
    ]
    for col in candidates:
        if col in adata.obs.columns:
            return col
    # Try regex search
    for col in adata.obs.columns:
        if re.search(r"(cell.?type|annotation|label)", col, re.I):
            return col
    return None


def filter_gc_b_cells(adata: ad.AnnData,
                      filter_proliferating: bool = False,
                      proliferation_threshold: float = 0.10,
                      exclude_proliferative_gc_types: bool = False) -> ad.AnnData:
    """Filter tonsil atlas for germinal center and memory B cell types, optionally filter proliferating cells.
    
    Includes: LZ/DZ B cells, centroblasts, centrocytes, GC B cells, memory B cells, and related populations.
    """
    # Convert to memory if backed (needed for filtering and copying)
    if adata.isbacked:
        print(f"  [tonsil] Converting backed AnnData to memory...")
        adata = adata.to_memory()
    
    label_col = find_label_column(adata)
    if label_col is None:
        print(f"  [tonsil] WARNING: No cell type column found!")
        print(f"  [tonsil] Available columns: {list(adata.obs.columns[:10])}")
        print(f"  [tonsil] Using all cells (no filtering)")
        return adata.copy()
    
    labels = adata.obs[label_col].astype(str).str.lower()
    
    # Comprehensive pattern matching for germinal center B cell types
    # LZ/DZ patterns
    is_gc = (
        labels.str.contains(r"\b(lz|light\s*zone)\b", case=False, na=False) |
        labels.str.contains(r"\b(dz|dark\s*zone)\b", case=False, na=False) |
        labels.str.contains(r"\bgcb.*(lz|dz)\b", case=False, na=False) |
        labels.str.contains(r"\b(lz|dz).*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"lz[-/]dz", case=False, na=False)
    )
    
    # Centroblast patterns
    is_gc |= (
        labels.str.contains(r"\bcentroblast", case=False, na=False) |
        labels.str.contains(r"\bcb\s*b\s*cells?\b", case=False, na=False)
    )
    
    # Centrocyte patterns
    is_gc |= (
        labels.str.contains(r"\bcentrocyte", case=False, na=False) |
        labels.str.contains(r"\bcc\s*b\s*cells?\b", case=False, na=False)
    )
    
    # General GC B cell patterns
    is_gc |= (
        labels.str.contains(r"\bgc\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bgerminal\s*center\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bgcb\b", case=False, na=False) |
        labels.str.contains(r"\bgc\s*b\s*lymphocyte", case=False, na=False)
    )
    
    # Memory B cell patterns
    is_gc |= (
        labels.str.contains(r"\bmemory\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bmbc\b", case=False, na=False) |
        labels.str.contains(r"\bmem\s*b\s*cells?\b", case=False, na=False) |
        labels.str.contains(r"\bb\s*memory", case=False, na=False)
    )
    
    # Activated B cells in GC context (but exclude plasma cells)
    is_gc |= (
        (labels.str.contains(r"\bactivated\s*b\s*cells?\b", case=False, na=False) |
         labels.str.contains(r"\bact\s*b\s*cells?\b", case=False, na=False)) &
        ~labels.str.contains(r"\bplasma\b", case=False, na=False)
    )
    
    n_matched = is_gc.sum()
    print(f"  [tonsil] Found {n_matched:,} GC + memory B cells (from {adata.n_obs:,} total)")
    
    if n_matched == 0:
        print(f"  [tonsil] WARNING: No GC/memory B cells found! Available cell types:")
        for ct, count in adata.obs[label_col].value_counts().head(15).items():
            print(f"    {ct}: {count:,} cells")
        print(f"  [tonsil] Using all cells instead")
        filtered = adata.copy()
    else:
        filtered = adata[is_gc].copy()
        # Report breakdown by cell type
        if label_col:
            print(f"  [tonsil] B cell type breakdown:")
            gc_types = filtered.obs[label_col].value_counts()
            for ct, count in gc_types.head(10).items():
                pct = 100 * count / len(filtered)
                print(f"    {ct}: {count:,} cells ({pct:.1f}%)")
        print(f"  [tonsil] Filtered to {filtered.n_obs:,} GC + memory B cells")
        
        # Print final breakdown
        if label_col:
            print(f"  [tonsil] Final composition:")
            final_types = filtered.obs[label_col].value_counts()
            for ct, count in final_types.head(10).items():
                 pct = 100 * count / len(filtered)
                 print(f"    - {ct}: {count:,} ({pct:.1f}%)")
    
    # Exclude proliferative GC cell types by annotation (e.g., "DZ late Sphase", "DZ early G2Mphase")
    if exclude_proliferative_gc_types and label_col and filtered.n_obs > 0:
        labels = filtered.obs[label_col].astype(str).str.lower()
        is_prolif_type = (
            labels.str.contains(r"\bproliferative\b", case=False, na=False) |
            labels.str.contains(r"\bs\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bg2m\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bg2\s*phase\b", case=False, na=False) |
            labels.str.contains(r"\bm\s*phase\b", case=False, na=False)
        )
        n_prolif_type = is_prolif_type.sum()
        if n_prolif_type > 0:
            filtered = filtered[~is_prolif_type].copy()
            print(f"  [tonsil] Excluded {n_prolif_type:,} cells with proliferative cell type annotations")
            print(f"  [tonsil] Retained {filtered.n_obs:,} non-proliferative B cells")
    
    # Filter proliferating cells if requested
    if filter_proliferating and filtered.n_obs > 0:
        print(f"  [tonsil] Computing cell cycle scores to filter proliferating cells...")
        # Preserve counts if not already in layers
        if "counts" not in filtered.layers:
            if _is_intlike(filtered.X):
                filtered.layers["counts"] = filtered.X.copy()
            elif filtered.raw is not None and _is_intlike(filtered.raw.X):
                filtered.layers["counts"] = filtered.raw.X.copy()
        # Ensure we have normalized data for cell cycle scoring
        if "normalized" not in filtered.layers:
            sc.pp.normalize_total(filtered, target_sum=1e4)
            sc.pp.log1p(filtered)
            filtered.layers["normalized"] = filtered.X.copy()
        
        # Cell cycle genes
        cc_genes = [x.strip() for x in """
MCM5,PCNA,TYMS,FEN1,MCM2,MCM4,RRM1,UNG,GINS2,MCM6,CDCA7,DTL,PRIM1,UHRF1,MLF1IP,HELLS,RFC2,RPA2,NASP,RAD51AP1,GMNN,WDR76,SLBP,CCNE2,UBR7,POLD3,MSH2,
ATAD2,RAD51,RRM2,CDC45,CDC6,EXO1,TIPIN,DSCC1,BLM,CASP8AP2,USP1,CLSPN,POLA1,CHAF1B,BRIP1,E2F8,
HMGB2,CDK1,NUSAP1,UBE2C,BIRC5,TPX2,TOP2A,NDC80,CKS2,NUF2,CKS1B,MKI67,TMPO,CENPF,TACC3,FAM64A,SMC4,CCNB2,CKAP2L,CKAP2,AURKB,BUB1,KIF11,ANP32E,TUBB4B,
GTSE1,KIF20B,HJURP,CDCA3,HN1,CDC20,TTK,CDC25C,KIF2C,RANGAP1,NCAPD2,DLGAP5,CDCA2,CDCA8,ECT2,KIF23,HMMR,AURKA,PSRC1,ANLN,LBR,CKAP5,CENPE,CTCF,NEK2,G2E3,GAS2L3,CBX5,CENPA
""".replace("\n",",").split(",") if x.strip()]
        s_genes = cc_genes[:43]
        g2m_genes = cc_genes[43:]
        
        # Score cell cycle
        sc.tl.score_genes_cell_cycle(filtered, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False)
        
        # Compute combined S+G2M score
        s_score = filtered.obs.get("S_score", pd.Series(0, index=filtered.obs.index))
        g2m_score = filtered.obs.get("G2M_score", pd.Series(0, index=filtered.obs.index))
        combined_cc_score = s_score + g2m_score
        filtered.obs["combined_cc_score"] = combined_cc_score
        
        # Filter out proliferating cells
        is_low_prolif = combined_cc_score <= proliferation_threshold
        n_before = filtered.n_obs
        filtered = filtered[is_low_prolif].copy()
        n_after = filtered.n_obs
        n_removed = n_before - n_after
        print(f"  [tonsil] Removed {n_removed:,} proliferating cells (S+G2M > {proliferation_threshold})")
        print(f"  [tonsil] Retained {n_after:,} low-proliferation B cells")
    
    return filtered


def load_and_filter_tonsil_combined(gc_path: str, mbc_path: str, max_cells: int,
                                    filter_proliferating: bool, proliferation_threshold: float,
                                    exclude_proliferative_gc_types: bool) -> ad.AnnData:
    """Load GC and MBC datasets using backed mode to save memory, filter, concatenate, and subsample."""
    parts = []
    
    # Helper to load subset from backed file
    def load_subset(path, tag):
        print(f"  [tonsil] Loading {tag} from: {Path(path).name} (backed mode)")
        try:
            # Load in backed mode - only reads metadata initially
            raw = sc.read_h5ad(path, backed="r")
            print(f"    Raw {tag}: {raw.n_obs:,} cells (on disk)")
            
            # Identify cells to keep based on obs (cell type)
            # We use filter_gc_b_cells but specialized to just return indices/mask first if possible?
            # Actually filter_gc_b_cells expects an AnnData. 
            # If we pass backed AnnData, it will read .obs which is fine.
            # But we need to ensure it doesn't try to read .X or copy the whole object.
            
            # Let's do the filtering manually here to ensure memory safety
            label_col = find_label_column(raw)
            if label_col is None:
                print(f"    (warn) No label column found, loading all...")
                subset = raw.to_memory()
            else:
                labels = raw.obs[label_col].astype(str).str.lower()
                
                # Re-use the regex logic (simplified/copied for safety or extract to helper if I could)
                # For now I will rely on the pattern matching being fast on just the obs series
                is_target = (
                    labels.str.contains(r"\b(lz|light\s*zone)\b", case=False, na=False) |
                    labels.str.contains(r"\b(dz|dark\s*zone)\b", case=False, na=False) |
                    labels.str.contains(r"\bgcb.*(lz|dz)\b", case=False, na=False) |
                    labels.str.contains(r"\b(lz|dz).*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"lz[-/]dz", case=False, na=False) |
                    labels.str.contains(r"\bcentroblast", case=False, na=False) |
                    labels.str.contains(r"\bcb\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bcentrocyte", case=False, na=False) |
                    labels.str.contains(r"\bcc\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgc\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgerminal\s*center\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bgcb\b", case=False, na=False) |
                    labels.str.contains(r"\bgc\s*b\s*lymphocyte", case=False, na=False) |
                    labels.str.contains(r"\bmemory\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bmbc\b", case=False, na=False) |
                    labels.str.contains(r"\bmem\s*b\s*cells?\b", case=False, na=False) |
                    labels.str.contains(r"\bb\s*memory", case=False, na=False) |
                    ((labels.str.contains(r"\bactivated\s*b\s*cells?\b", case=False, na=False) |
                      labels.str.contains(r"\bact\s*b\s*cells?\b", case=False, na=False)) &
                     ~labels.str.contains(r"\bplasma\b", case=False, na=False))
                )
                
                n_match = is_target.sum()
                print(f"    Found {n_match:,} target cells in {tag}")
                
                if n_match == 0:
                    print("    (warn) No target cells found, skipping...")
                    return None
                
                # Subset in backed mode (lazy)
                subset_backed = raw[is_target]
                
                # Now load ONLY the subset into memory
                print(f"    Loading {n_match:,} cells into memory...")
                subset = subset_backed.to_memory()
                
            # Close file handle
            if hasattr(raw, 'file') and raw.file is not None:
                try: raw.file.close()
                except: pass
                
            # Post-load processing
            subset = harmonize_human_genes(subset, tag=f"tonsil_{tag}")
            subset = enforce_raw_counts(subset, tag=f"tonsil_{tag}")
            
            # Apply fine-grained filtering (proliferation, etc.) in memory
            # We reuse filter_gc_b_cells but now it works on a much smaller object
            filt = filter_gc_b_cells(subset,
                                     filter_proliferating=filter_proliferating,
                                     proliferation_threshold=proliferation_threshold,
                                     exclude_proliferative_gc_types=exclude_proliferative_gc_types)
            return filt
            
        except Exception as e:
            print(f"    (error) Failed to load {tag}: {e}")
            return None
            
    # 1. Load GC B cells
    gc_filt = load_subset(gc_path, "gc")
    if gc_filt is not None:
        parts.append(gc_filt)
    gc.collect()
    
    # 2. Load Memory B cells
    mbc_filt = load_subset(mbc_path, "mbc")
    if mbc_filt is not None:
        parts.append(mbc_filt)
    gc.collect()
        
    if not parts:
        return None
        
    # 3. Concatenate
    print("  [tonsil] Concatenating datasets...")
    combined = ad.concat(parts, join="outer", index_unique=None, fill_value=0)
    # Restore counts
    combined.layers["counts"] = combined.layers.get("counts", combined.X)
    
    # Clean up parts to free memory
    del parts, gc_filt, mbc_filt
    gc.collect()
    
    print(f"  [tonsil] Combined GC + MBC: {combined.n_obs:,} cells")
    
    # 4. Subsample
    if combined.n_obs > max_cells:
        np.random.seed(0)
        idx = np.random.choice(combined.n_obs, size=max_cells, replace=False)
        combined = combined[idx].copy()
        print(f"  [tonsil] Subsampled to {max_cells:,} cells")
        
    return combined


# ============================== SCRUBLET =====================================
def run_scrublet_safe(adata: ad.AnnData, sample_name: str) -> pd.Series:
    """Run Scrublet for doublet detection with proper error handling."""
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
        scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            verbose=False
        )
        
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


# ============================== QC FUNCTIONS =================================
def compute_qc_metrics(adata: ad.AnnData, species: str = "mouse") -> ad.AnnData:
    """Compute QC metrics including mitochondrial gene percentage."""
    if species == "mouse":
        adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    else:  # human
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    
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
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= mt_threshold].copy()
    n_after = adata.n_obs
    print(f"    Filtered: {n_before:,} → {n_after:,} cells "
          f"({n_before - n_after:,} removed)")
    return adata


def load_cellbender_h5(h5_path: Path, sample_name: str) -> ad.AnnData:
    """Load CellBender filtered h5 file."""
    print(f"  Loading: {h5_path.name}")
    
    try:
        adata = sc.read_10x_h5(str(h5_path))
    except Exception as e1:
        print(f"    read_10x_h5 failed: {e1}, trying read_h5ad...")
        adata = sc.read_h5ad(str(h5_path))
    
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_name
    
    if not adata.obs_names.str.contains("-").any():
        adata.obs_names = [f"{bc}-{sample_name}" for bc in adata.obs_names]
    else:
        adata.obs_names = [f"{bc.split('-')[0]}-{sample_name}" for bc in adata.obs_names]
    
    print(f"    Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


# ============================== HVG FUNCTIONS ================================
def confounder_mask(varnames: pd.Index) -> pd.Series:
    """Create mask for confounder genes."""
    v = pd.Index([str(g) for g in varnames])
    is_mt = v.str.startswith("MT-")
    is_ribo = v.str.startswith(("RPS", "RPL"))
    is_ig = v.str.match(r"^(IGH|IGK|IGL)[VDJC].*", na=False)
    tcr_prefixes = ("TRAV", "TRBV", "TRGV", "TRDV", "TRAJ", "TRBJ", 
                    "TRGJ", "TRDJ", "TRAC", "TRBC", "TRGC", "TRDC")
    is_tcr = v.str.startswith(tcr_prefixes)
    return is_mt | is_ribo | is_ig | is_tcr


def build_sample_batch(adx, tag=""):
    """Build sample_batch column for HVG detection."""
    for cand in ["sample_id", "donor", "library_id", "orig.ident", "sample", 
                 "Sample", "Library", "Donor", "batch"]:
        if cand in adx.obs:
            adx.obs["sample_batch"] = adx.obs[cand].astype(str)
            print(f"  [{tag}] sample_batch ← '{cand}'")
            return adx
    adx.obs["sample_batch"] = (adx.obs.get("study", "study").astype(str) + "_" +
                               adx.obs.get("sample_type", "sample").astype(str))
    print(f"  [{tag}] sample_batch ← fallback")
    return adx


def compute_shared_hvgs(adata_list, n_top=3500, min_datasets=2):
    """Compute HVGs shared across datasets."""
    hvg_counts = {}
    for i, adx in enumerate(adata_list):
        print(f"    HVGs for dataset {i+1}/{len(adata_list)} (seurat_v3 on counts)")
        bk = "sample_batch" if "sample_batch" in adx.obs.columns else None
        sc.pp.highly_variable_genes(
            adx, layer="counts", n_top_genes=n_top, flavor="seurat_v3",
            subset=False, batch_key=bk
        )
        for g in adx.var_names[adx.var["highly_variable"]]:
            hvg_counts[g] = hvg_counts.get(g, 0) + 1
    min_count = max(min_datasets, len(adata_list) // 2)
    shared = [g for g, c in hvg_counts.items() if c >= min_count]
    print(f"  Genes HVG in ≥{min_count} datasets: {len(shared):,}")
    return set(shared)


# ============================== CT2 FUNCTIONS ================================
def stream_counts_to_ct2_txt(adata: ad.AnnData, dest: Path) -> None:
    """Write CytoTRACE2-compatible matrix."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    
    if "counts" not in adata.layers:
        raise RuntimeError("[ct2-io] layers['counts'] missing.")
    
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
            g = str(adata.var_names[j]).strip().upper()
            if g in na_like:
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
        
        id_col = next((c for c in ("cell", "cell_id", "barcode") if c in df.columns), df.columns[0])
        score_cols = [c for c in df.columns if ("cytotrace2" in c and "score" in c)]
        
        if not score_cols:
            continue
        
        tmp = pd.DataFrame(index=df[id_col].astype(str).values)
        tmp["cytotrace2_score"] = _coerce_float_series(df[score_cols[0]]).values
        
        pot_cols = [c for c in df.columns if "potency" in c]
        if pot_cols:
            tmp["cytotrace2_potency"] = df[pot_cols[0]].astype(str).values
        
        n_match = len(set(tmp.index) & set(map(str, adata.obs_names)))
        if best is None or n_match > best[0]:
            best = (n_match, f, tmp)
    
    if best is None:
        raise FileNotFoundError(f"[ct2] No results found under {outdir}")
    
    print(f"[ct2] Using: {best[1]} (matched {best[0]:,} cells)")
    df_best = best[2]
    df_best = df_best[~df_best.index.duplicated(keep="first")]
    return df_best.reindex(adata.obs_names.astype(str))


# ============================== PLOTTING FUNCTIONS ===========================
def _save_umap(adata, color, fname, title=None, palette=None, cmap=None):
    """Save UMAP plot in PNG, PDF, and SVG formats."""
    fig, ax = plt.subplots(figsize=(10, 9))
    sc.pl.umap(adata, color=color, title=(title or color), palette=palette, cmap=cmap,
               frameon=False, legend_loc="right margin", ax=ax, show=False, s=25)
    plt.tight_layout()
    
    # Save in multiple formats
    base_name = fname.replace('.png', '')
    plt.savefig(FIGDIR / f"{base_name}.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIGDIR / f"{base_name}.pdf", bbox_inches="tight")
    plt.savefig(FIGDIR / f"{base_name}.svg", bbox_inches="tight")
    plt.close()
    print(f"  ✓ {base_name} (.png, .pdf, .svg)")


def save_umaps_by_sample(adata: ad.AnnData, figdir: Path, group_key: str = "sample_batch"):
    """Save individual UMAPs highlighting each sample in PNG, PDF, and SVG formats."""
    if group_key not in adata.obs:
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
            adata, color="_highlight", palette=palette, frameon=False,
            legend_loc=None, show=False, s=25, ax=ax, title=f"{group_key} → {grp}"
        )
        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.pdf", bbox_inches="tight")
        fig.savefig(outdir / f"umap_{group_key}_{safe_grp}.svg", bbox_inches="tight")
        plt.close(fig)
        adata.obs.drop(columns="_highlight", inplace=True)
    
    print(f"  ✓ Per-group UMAPs saved to {outdir} (.png, .pdf, .svg)")


def save_individual_sample_umaps(adata: ad.AnnData, figdir: Path, 
                                  group_key: str = "sample_batch"):
    """Save UMAPs showing ONLY cells from each sample in PNG, PDF, and SVG formats."""
    if group_key not in adata.obs:
        return
    
    outdir = figdir / f"umap_individual_{group_key}"
    outdir.mkdir(parents=True, exist_ok=True)
    
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    
    # Get global UMAP limits
    umap_coords = adata.obsm["X_umap"]
    x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
    y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    
    for grp in unique_groups:
        mask = labels == grp
        adata_sub = adata[mask].copy()
        n_cells = adata_sub.n_obs
        
        safe_grp = re.sub(r"[^A-Za-z0-9._-]+", "_", grp)
        
        fig, ax = plt.subplots(figsize=(9, 8))
        ax.scatter(
            adata_sub.obsm["X_umap"][:, 0],
            adata_sub.obsm["X_umap"][:, 1],
            c="#d62728", s=15, alpha=0.7, rasterized=True
        )
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_title(f"{grp} (n={n_cells:,})", fontsize=14)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.pdf", bbox_inches="tight")
        fig.savefig(outdir / f"umap_{safe_grp}_cells_only.svg", bbox_inches="tight")
        plt.close(fig)
    
    print(f"  ✓ Individual sample UMAPs saved to {outdir} (.png, .pdf, .svg)")


def add_ct2_umap_and_plots(adata: ad.AnnData, figdir: Path, title_suffix=""):
    """Add CytoTRACE2 UMAP and histogram plots in PNG, PDF, and SVG formats."""
    figdir.mkdir(parents=True, exist_ok=True)
    
    if "cytotrace2_score" in adata.obs:
        fig, ax = plt.subplots(figsize=(10, 9))
        sc.pl.umap(adata, color="cytotrace2_score", ax=ax, show=False, s=30, frameon=False,
                   cmap="viridis", title=f"CytoTRACE2 score {title_suffix}")
        fig.savefig(figdir / "umap_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        fig.savefig(figdir / "umap_cytotrace2_score.pdf", bbox_inches="tight")
        fig.savefig(figdir / "umap_cytotrace2_score.svg", bbox_inches="tight")
        plt.close(fig)
        
        fig, ax = plt.subplots(figsize=(7, 5))
        adata.obs["cytotrace2_score"].astype(float).plot(kind="hist", bins=50, ax=ax)
        ax.set_xlabel("CytoTRACE2 score")
        ax.set_ylabel("Cell count")
        ax.set_title("Distribution of CT2 scores")
        fig.savefig(figdir / "hist_cytotrace2_score.png", dpi=300, bbox_inches="tight")
        fig.savefig(figdir / "hist_cytotrace2_score.pdf", bbox_inches="tight")
        fig.savefig(figdir / "hist_cytotrace2_score.svg", bbox_inches="tight")
        plt.close(fig)


# ============================== MAIN PIPELINE ================================
print("=" * 84)
print("STEP 0 — Ortholog mapping (BioMart preferred)")
print("=" * 84)
M2H_1TO1 = get_m2h_map(prefer_biomart=True)

# ============================== LOAD MOUSE MALIGNANT =========================
print("\n" + "=" * 84)
print("STEP 1 — Load and QC mouse malignant samples (CellBender filtered)")
print("=" * 84)

mouse_list = []
qc_stats = []

for sample_name, sample_info in MOUSE_MALIGNANT_FILES.items():
    h5_path = sample_info["path"]
    
    if not h5_path.exists():
        print(f"  WARNING: {h5_path} not found, skipping...")
        continue
    
    # Load sample
    adata = load_cellbender_h5(h5_path, sample_name)
    
    # Add metadata
    adata.obs["condition"] = sample_info["condition"]
    adata.obs["replicate"] = sample_info["replicate"]
    adata.obs["species"] = "mouse"
    adata.obs["study"] = "mouse_malignant"
    adata.obs["disease_state"] = f"Mouse_{sample_info['condition']}"
    
    # Compute QC metrics
    adata = compute_qc_metrics(adata, species="mouse")
    n_pre = adata.n_obs
    
    # Run Scrublet
    doublet_mask = run_scrublet_safe(adata, sample_name)
    adata.obs["predicted_doublet"] = doublet_mask.values
    
    # Filter cells
    adata = filter_cells(adata, mt_threshold=MT_THRESHOLD, min_genes=MIN_GENES)
    
    # Remove doublets
    n_pre_doublet = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_doublets_removed = n_pre_doublet - adata.n_obs
    print(f"    Removed {n_doublets_removed:,} doublets → {adata.n_obs:,} cells")
    
    # Store counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    
    qc_stats.append({
        "sample": sample_name,
        "condition": sample_info["condition"],
        "species": "mouse",
        "cells_raw": n_pre,
        "cells_after_qc": adata.n_obs
    })
    
    mouse_list.append(adata)
    print(f"    Final: {adata.n_obs:,} cells\n")

# Concatenate mouse samples
mouse_adata = ad.concat(mouse_list, join="outer", index_unique=None)
print(f"  Mouse combined: {mouse_adata.n_obs:,} cells × {mouse_adata.n_vars:,} genes")

# Ensure counts layer
if "counts" not in mouse_adata.layers:
    mouse_adata.layers["counts"] = mouse_adata.X.copy()

# Convert mouse genes to human
mouse_adata = convert_mouse_to_human_genes(mouse_adata, M2H_1TO1, tag="mouse")
mouse_adata = enforce_raw_counts(mouse_adata, tag="mouse")

# ============================== LOAD HUMAN DLBCL =============================
print("\n" + "=" * 84)
print("STEP 2 — Load human DLBCL samples")
print("=" * 84)

dlbcl_parts = []
for p in DLBCL_FILES:
    if not Path(p).exists():
        print(f"  WARNING: {p} not found, skipping...")
        continue
    
    a = sc.read_h5ad(p)
    nm = Path(p).stem
    a.obs["species"] = "human"
    a.obs["study"] = ("DLBCL_Alizadeh" if "Alizadeh" in p else "DLBCL_Roider")
    a.obs["sample_type"] = "dlbcl"
    a.obs["disease_state"] = "DLBCL"
    a = harmonize_human_genes(a, tag=nm)
    a = enforce_raw_counts(a, tag=nm)
    if "sample" not in a.obs:
        a.obs["sample"] = nm
    dlbcl_parts.append(a)
    print(f"  + {Path(p).name}: {a.n_obs:,} cells, {a.n_vars:,} genes")
    
    qc_stats.append({
        "sample": nm,
        "condition": "DLBCL",
        "species": "human",
        "cells_raw": a.n_obs,
        "cells_after_qc": a.n_obs
    })

dlbcl = ad.concat(dlbcl_parts, join="outer", index_unique=None)
dlbcl.layers["counts"] = dlbcl.layers.get("counts", dlbcl.X)
print(f"  DLBCL combined: {dlbcl.n_obs:,} × {dlbcl.n_vars:,}")

# ============================== LOAD HUMAN TONSIL ============================
print("\n" + "=" * 84)
print("STEP 2.5 — Load Human Tonsil GC B cells")
print("=" * 84)

# Load and filter Tonsil data (GC + MBC combined)
tonsil = load_and_filter_tonsil_combined(
    gc_path=TONSIL_GC_PATH,
    mbc_path=TONSIL_MBC_PATH,
    max_cells=MAX_TONSIL_CELLS,
    filter_proliferating=FILTER_PROLIFERATING_TONSIL,
    proliferation_threshold=PROLIFERATION_THRESHOLD,
    exclude_proliferative_gc_types=EXCLUDE_PROLIFERATIVE_GC_TYPES
)

if tonsil is None:
    raise RuntimeError("Failed to load tonsil data.")

tonsil.obs["species"] = "human"
tonsil.obs["study"] = "Tonsil_GCBC_MBC"
tonsil.obs["sample_type"] = "tonsil_b_cells"
tonsil.obs["disease_state"] = "Tonsil_Normal"
print(f"  Tonsil (GC+MBC) final: {tonsil.n_obs:,} × {tonsil.n_vars:,}")

qc_stats.append({
    "sample": "Tonsil_GC_MBC",
    "condition": "Tonsil_Normal",
    "species": "human",
    "cells_raw": "N/A (combined)",
    "cells_after_qc": tonsil.n_obs
})

# Save QC stats
qc_df = pd.DataFrame(qc_stats)
qc_df.to_csv(OUTDIR / "qc_statistics.csv", index=False)
print(f"\n  QC stats saved to: {OUTDIR / 'qc_statistics.csv'}")
print(qc_df.to_string())

# ============================== NORMALIZE ====================================
print("\n" + "=" * 84)
print("STEP 3 — Normalize datasets")
print("=" * 84)

mouse_adata = harmonize_preprocessing(mouse_adata, tag="mouse")
dlbcl = harmonize_preprocessing(dlbcl, tag="DLBCL")
tonsil = harmonize_preprocessing(tonsil, tag="tonsil")

# ============================== COMMON GENES + HVGs ==========================
print("\n" + "=" * 84)
print("STEP 4 — Common genes + HVGs")
print("=" * 84)

datasets = [mouse_adata, dlbcl, tonsil]
common_genes = sorted(list(set(datasets[0].var_names).intersection(
    *[set(x.var_names) for x in datasets[1:]])))
print(f"  Common genes across datasets: {len(common_genes):,}")

# Build sample_batch
mouse_adata = build_sample_batch(mouse_adata, "mouse")
dlbcl = build_sample_batch(dlbcl, "DLBCL")
tonsil = build_sample_batch(tonsil, "tonsil")

adata_for_hvg = [x[:, common_genes].copy() for x in datasets]
hvgs = compute_shared_hvgs(adata_for_hvg, n_top=4000, min_datasets=2)
hvgi = pd.Index(sorted(hvgs))
hvgs_clean = sorted(list(hvgi[~confounder_mask(hvgi)]))
print(f"  HVGs after removing confounders: {len(hvgs_clean):,}")

# ============================== CONCATENATE ==================================
print("\n" + "=" * 84)
print("STEP 5 — Concatenate datasets")
print("=" * 84)

# For scVI: use common genes (inner join)
adata_full = ad.concat([x[:, common_genes] for x in datasets], 
                        merge="same", join="inner", index_unique=None)
if "counts" not in adata_full.layers:
    raise RuntimeError("layers['counts'] missing in adata_full.")
print(f"  Combined FULL (common genes): {adata_full.n_obs:,} × {adata_full.n_vars:,}")

# For CT2: use UNION of all genes (outer join) - CT2 explicitly asks NOT to pre-filter genes
if RUN_CYTOTRACE2:
    print("  Creating CT2 object with UNION of all genes (outer join, fill missing with zeros)...")
    all_genes_union = sorted(list(set().union(*[set(x.var_names) for x in datasets])))
    print(f"  Union of all genes: {len(all_genes_union):,} (vs {len(common_genes):,} common genes)")

    # Verify all datasets have counts layer before concatenation
    for i, ds in enumerate(datasets):
        if "counts" not in ds.layers:
            raise RuntimeError(f"[ct2] Dataset {i} missing layers['counts'].")

    # Concatenate with outer join - fills missing genes with zeros
    adata_ct2 = ad.concat(datasets, merge="same", join="outer", index_unique=None, fill_value=0)

    # Verify/reconstruct counts layer
    if "counts" not in adata_ct2.layers:
        if _is_intlike(adata_ct2.X):
            adata_ct2.layers["counts"] = adata_ct2.X.copy()
            print("  [ct2] Reconstructed counts layer from X")
        else:
            raise RuntimeError("[ct2] Could not find or reconstruct raw counts for CT2.")
    print(f"  CT2 object: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes (UNION, no pre-filtering)")
else:
    print("  Skipping CT2 object creation (RUN_CYTOTRACE2=False)")
    adata_ct2 = None

# Cross-species HVG filtering
if "species" in adata_full.obs and adata_full.obs["species"].nunique() > 1:
    print("  Ensuring HVGs detected across species (≥5% of cells each)")
    hvg_view = adata_full[:, hvgs_clean]
    counts_layer = hvg_view.layers["counts"] if "counts" in hvg_view.layers else hvg_view.X
    is_sparse = sp.issparse(counts_layer)
    if is_sparse:
        counts_layer = counts_layer.tocsr()
    species_labels = hvg_view.obs["species"].astype(str).values
    keep_mask = np.ones(len(hvgs_clean), dtype=bool)
    for sp_name in np.unique(species_labels):
        row_idx = np.where(species_labels == sp_name)[0]
        if row_idx.size == 0:
            continue
        if is_sparse:
            sub = counts_layer[row_idx, :]
            detected = np.asarray((sub > 0).mean(axis=0)).ravel()
        else:
            sub = counts_layer[row_idx, :]
            detected = (sub > 0).astype(float).mean(axis=0)
        keep_mask &= detected >= 0.05
    hvgs_filtered = [g for g, keep in zip(hvgs_clean, keep_mask) if keep]
    if hvgs_filtered:
        hvgs_clean = hvgs_filtered
        print(f"  HVGs after cross-species detection ≥5%: {len(hvgs_clean):,}")

adata_all = adata_full[:, hvgs_clean].copy()
print(f"  Combined HVG: {adata_all.n_obs:,} × {adata_all.n_vars:,}")

# ============================== scVI =========================================
print("\n" + "=" * 84)
print("STEP 6 — scVI integration (NO cell cycle regression)")
print("=" * 84)

for cat in ("sample_batch", "species"):
    if cat in adata_all.obs:
        adata_all.obs[cat] = adata_all.obs[cat].astype("category")

scvi.model.SCVI.setup_anndata(
    adata_all, layer="counts", batch_key="sample_batch",
    categorical_covariate_keys=["species"]
)

model = scvi.model.SCVI(
    adata_all,
    n_latent=96,
    n_layers=2,
    dropout_rate=0.15,
    gene_likelihood="nb",
    dispersion="gene-batch",
    use_layer_norm="both",
    use_batch_norm="none"
)

max_epochs = 150 if torch.cuda.is_available() else 600
print(f"  Training scVI for {max_epochs} epochs (early stopping patience=60)...")
model.train(
    max_epochs=max_epochs,
    plan_kwargs={"lr": 3e-4, "n_epochs_kl_warmup": 200, "reduce_lr_on_plateau": True},
    check_val_every_n_epoch=20,
    early_stopping=True,
    early_stopping_patience=60
)

adata_all.obsm["X_scvi"] = model.get_latent_representation()

# ============================== UMAP + CLUSTERING ============================
print("\n" + "=" * 84)
print("STEP 7 — UMAP and clustering")
print("=" * 84)

sc.pp.neighbors(adata_all, use_rep="X_scvi", n_neighbors=30)
sc.tl.umap(adata_all, min_dist=0.2, spread=1.5)

for res in [0.5, 1.0, 1.5]:
    sc.tl.leiden(adata_all, resolution=res, key_added=f"leiden_{res}")
    print(f"  Leiden res={res}: {adata_all.obs[f'leiden_{res}'].nunique()} clusters")

# ============================== CELL CYCLE + GENE SET SCORES =================
print("\n" + "=" * 84)
print("STEP 8 — Cell cycle scoring + gene-set scores")
print("=" * 84)

# Cell cycle scoring (for visualization only, NOT regressed out)
cc_genes = [x.strip() for x in """
MCM5,PCNA,TYMS,FEN1,MCM2,MCM4,RRM1,UNG,GINS2,MCM6,CDCA7,DTL,PRIM1,UHRF1,MLF1IP,HELLS,RFC2,RPA2,NASP,RAD51AP1,GMNN,WDR76,SLBP,CCNE2,UBR7,POLD3,MSH2,
ATAD2,RAD51,RRM2,CDC45,CDC6,EXO1,TIPIN,DSCC1,BLM,CASP8AP2,USP1,CLSPN,POLA1,CHAF1B,BRIP1,E2F8,
HMGB2,CDK1,NUSAP1,UBE2C,BIRC5,TPX2,TOP2A,NDC80,CKS2,NUF2,CKS1B,MKI67,TMPO,CENPF,TACC3,FAM64A,SMC4,CCNB2,CKAP2L,CKAP2,AURKB,BUB1,KIF11,ANP32E,TUBB4B,
GTSE1,KIF20B,HJURP,CDCA3,HN1,CDC20,TTK,CDC25C,KIF2C,RANGAP1,NCAPD2,DLGAP5,CDCA2,CDCA8,ECT2,KIF23,HMMR,AURKA,PSRC1,ANLN,LBR,CKAP5,CENPE,CTCF,NEK2,G2E3,GAS2L3,CBX5,CENPA
""".replace("\n", ",").split(",") if x.strip()]
s_genes = cc_genes[:43]
g2m_genes = cc_genes[43:]

# Filter to genes present
s_genes_present = [g for g in s_genes if g in adata_all.var_names]
g2m_genes_present = [g for g in g2m_genes if g in adata_all.var_names]

if s_genes_present and g2m_genes_present:
    sc.tl.score_genes_cell_cycle(adata_all, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
    print(f"  Cell cycle scoring: {len(s_genes_present)} S genes, {len(g2m_genes_present)} G2M genes")
    print(f"  Phase distribution:")
    print(adata_all.obs['phase'].value_counts().to_string())
else:
    print(f"  WARNING: Not enough cell cycle genes found for scoring")

# Gene set scores
avail = set(adata_all.var_names)
ox = [g for g in OXPHOS_GENES if g in avail]
bcr = [g for g in BCR_GENES if g in avail]
if ox:
    sc.tl.score_genes(adata_all, gene_list=ox, score_name='oxphos_score')
    print(f"  OXPHOS score: {len(ox)} genes")
if bcr:
    sc.tl.score_genes(adata_all, gene_list=bcr, score_name='bcr_score')
    print(f"  BCR score: {len(bcr)} genes")

# ============================== CytoTRACE2 ===================================
print("\n" + "=" * 84)
print("STEP 9 — CytoTRACE2")
print("=" * 84)

if RUN_CYTOTRACE2:
    # Use adata_ct2 (union of all genes) for CT2, as recommended by CT2 authors
    print("  Using CT2 object with UNION of all genes (no pre-filtering, as recommended by CT2)")
    print(f"  CT2 input size: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes")

    # Filter to genes expressed in at least some cells to reduce file size/memory
    if adata_ct2.n_vars > 30000:
        print(f"  [ct2] Filtering to genes expressed in ≥0.1% of cells to reduce memory...")
        X_counts = adata_ct2.layers["counts"] if "counts" in adata_ct2.layers else adata_ct2.X
        if sp.issparse(X_counts):
            n_expressed = np.array((X_counts > 0).sum(axis=0)).ravel()
        else:
            n_expressed = (X_counts > 0).sum(axis=0)
        min_cells = max(1, int(adata_ct2.n_obs * 0.001))  # At least 0.1% of cells
        keep_genes = n_expressed >= min_cells
        n_keep = keep_genes.sum()
        print(f"  [ct2] Keeping {n_keep:,} genes (expressed in ≥{min_cells} cells) out of {adata_ct2.n_vars:,}")
        adata_ct2 = adata_ct2[:, keep_genes].copy()
        print(f"  [ct2] Filtered CT2 input: {adata_ct2.n_obs:,} cells × {adata_ct2.n_vars:,} genes")

    stream_counts_to_ct2_txt(adata_ct2, CT2_INPUT_TXT)

    ct2_obj = run_ct2_python(CT2_INPUT_TXT, species="human", outdir=CT2_OUTDIR)

    if ct2_obj is None:
        cmd = f"cytotrace2 -f {shlex.quote(str(CT2_INPUT_TXT))} -sp human --output-dir {shlex.quote(str(CT2_OUTDIR))} --disable-plotting"
        print("[ct2] CLI:", cmd)
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode != 0:
            print(f"[ct2] WARNING: CytoTRACE2 CLI failed")

    # Parse CT2 scores
    df_scores = None
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
            idx = adata_ct2.obs_names.astype(str)
            
            df_scores = pd.DataFrame(index=idx)
            sname = next((c for c in obs.columns if ("cytotrace2" in c and "score" in c)), None)
            if sname:
                df_scores["cytotrace2_score"] = _coerce_float_series(obs.reindex(idx)[sname]).values
            
            p = next((c for c in obs.columns if "potency" in c), None)
            if p:
                df_scores["cytotrace2_potency"] = obs.reindex(idx)[p].astype(str).values
        except Exception as e:
            print(f"[ct2] Could not parse object: {e}")

    if df_scores is None:
        try:
            df_scores = parse_ct2_scores(CT2_OUTDIR, adata_ct2)
        except FileNotFoundError:
            print("[ct2] WARNING: CT2 results not found")

    if df_scores is not None and not df_scores.empty:
        for col in df_scores.columns:
            ser_full = df_scores[col].reindex(adata_full.obs_names.astype(str))
            ser_all = df_scores[col].reindex(adata_all.obs_names.astype(str))
            adata_full.obs[col] = ser_full
            adata_all.obs[col] = ser_all
        print("✓ CytoTRACE2 scores attached")
else:
    print("  Skipping CytoTRACE2 execution (RUN_CYTOTRACE2=False)")


# ============================== FIGURES ======================================
print("\n" + "=" * 84)
print("STEP 10 — Generate figures")
print("=" * 84)

disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000",
    "Tonsil_GC_B": "#4169E1",
    "Tonsil_Normal": "#4169E1"  # Same color as Tonsil_GC_B (royal blue)
}
species_colors = {"mouse": "#98FB98", "human": "#6495ED"}
phase_colors = {"G1": "#3498db", "S": "#e74c3c", "G2M": "#2ecc71"}

_save_umap(adata_all, "disease_state", "umap_disease_state.png", "Disease State", disease_colors)
_save_umap(adata_all, "species", "umap_species.png", "Species", species_colors)
_save_umap(adata_all, "sample_batch", "umap_sample_batch.png", "Sample Batch")
_save_umap(adata_all, "leiden_1.0", "umap_leiden_1.0.png", "Leiden (res=1.0)")
_save_umap(adata_all, "phase", "umap_cell_cycle.png", "Cell Cycle Phase", phase_colors)

if "oxphos_score" in adata_all.obs:
    _save_umap(adata_all, "oxphos_score", "umap_oxphos_score.png", "OXPHOS Score", cmap="RdYlBu_r")
if "bcr_score" in adata_all.obs:
    _save_umap(adata_all, "bcr_score", "umap_bcr_score.png", "BCR Score", cmap="RdYlBu_r")

save_umaps_by_sample(adata_all, FIGDIR, group_key="sample_batch")
save_individual_sample_umaps(adata_all, FIGDIR, group_key="sample_batch")

if "cytotrace2_score" in adata_all.obs:
    add_ct2_umap_and_plots(adata_all, FIGDIR_CT2, title_suffix="(integrated)")
    # Also save CT2 on the main UMAP
    _save_umap(adata_all, "cytotrace2_score", "umap_cytotrace2_score.png", 
               "CytoTRACE2 Score (scVI UMAP)", cmap="viridis")
elif RUN_CYTOTRACE2:
    print("  (warn) CytoTRACE2 scores missing despite RUN_CYTOTRACE2=True")
else:
    print("  Skipping CytoTRACE2 plots")

# ============================== SAVE =========================================
print("\n" + "=" * 84)
print("STEP 11 — Save outputs")
print("=" * 84)

try:
    model.save(OUTDIR / "scvi_model", overwrite=True)
    print(f"  ✓ scVI model → {OUTDIR / 'scvi_model'}")
except Exception as e:
    print(f"  (warn) Could not save scVI model: {e}")

adata_all.write(OUTDIR / "integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")
print(f"  ✓ AnnData → {OUTDIR / 'integrated_human_dlbcl_mouse_malignant_tonsil.h5ad'}")

adata_all.obs.to_csv(OUTDIR / "cell_metadata.csv")
print(f"  ✓ Metadata → {OUTDIR / 'cell_metadata.csv'}")

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  Output directory: {OUTDIR}")
print(f"  Figures: {FIGDIR}")
print(f"  CT2 results: {CT2_OUTDIR}")
print(f"  Total cells: {adata_all.n_obs:,}")
print(f"  Total genes: {adata_all.n_vars:,}")

# Summary by disease state
print("\n  Cells by disease state:")
for ds, count in adata_all.obs["disease_state"].value_counts().items():
    print(f"    {ds}: {count:,}")

print("\nDONE.\n")



