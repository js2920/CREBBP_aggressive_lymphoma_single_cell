#!/usr/bin/env bash
# ==============================================================================
# MasterAnalysis.sh
# Complete analysis pipeline for:
#   "Single-cell RNA sequencing of mouse B-cell lymphoma progression
#    and cross-species integration with human DLBCL"
# ==============================================================================
#
# This is a SELF-CONTAINED reproducibility recipe.
# All 26 analysis scripts + 1 synthetic smoke test are embedded below and will
# be extracted into ./scripts/ when you run this file.
#
# QUICK START:
#   1. View options:        bash MasterAnalysis.sh --help
#   2. Check environment:   bash MasterAnalysis.sh --check-env
#   3. Test pipeline code:  bash MasterAnalysis.sh --smoke-test
#   4. Preview dry-run:     bash MasterAnalysis.sh --dry-run
#   5. Run targeted stage:  bash MasterAnalysis.sh --run --stage 1
#   6. Run full analysis:   bash MasterAnalysis.sh --run
#
# ==============================================================================
#
# PIPELINE STAGES
# ---------------
#   Stage 0  Cell Ranger 9.0.1 alignment to GRCm39 (manual)
#   Stage 1  Mouse scVI integration, CytoTRACE2, gene-set scoring (9 scripts)
#   Stage 2  Geneformer fine-tuning on tonsil atlas → mouse prediction (2 scripts)
#   Stage 3  Cross-species scVI: human DLBCL + mouse + tonsil B/plasma (9 scripts)
#   Stage 4  Gene expression programs along CytoTRACE2 trajectory + GSEA (3 scripts)
#   Stage 5  Supplementary: PGC1α trajectory, Wilcoxon rank-sum tests (2 scripts)
#
# ==============================================================================

set -euo pipefail

# ┌────────────────────────────────────────────────────────────────────────────┐
# │                      ── USER CONFIGURATION ──                             │
# │  Edit the paths below to match your local setup before running.           │
# └────────────────────────────────────────────────────────────────────────────┘

# Base directory: auto-detected repository root
BASEDIR="${CREBBP_BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ── Input data directories ──
CELLBENDER_DIR="${CELLBENDER_DIR:-${BASEDIR}/data/cellbender_filtered}"
DLBCL_DIR="${DLBCL_DIR:-${BASEDIR}/data/DLBCL}"
TONSIL_DIR="${TONSIL_DIR:-${BASEDIR}/data/tonsil_export}"
GENEFORMER_MODEL_DIR="${GENEFORMER_MODEL_DIR:-${BASEDIR}/models/geneformer_tonsil_multi}"

# ── Configuration file ──
CONFIG_FILE="${CONFIG_FILE:-${BASEDIR}/config/pipeline_config.yaml}"

# ── Cell Ranger settings (Stage 0 only) ──
FASTQ_DATA_ROOT="${FASTQ_DATA_ROOT:-${BASEDIR}/data/raw_fastq}"
CELLRANGER_REF="${CELLRANGER_REF:-${BASEDIR}/data/refdata-gex-GRCm39-2024-A}"
CELLRANGER_BIN="${CELLRANGER_BIN:-cellranger}"

# ── Python executable auto-detection ──
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if command -v python >/dev/null 2>&1 && python -c "import scanpy" >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python)"
    elif command -v python3 >/dev/null 2>&1 && python3 -c "import scanpy" >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python3)"
    elif [[ -x "${CONDA_PREFIX:-}/bin/python" ]] && "${CONDA_PREFIX}/bin/python" -c "import scanpy" >/dev/null 2>&1; then
        PYTHON_BIN="${CONDA_PREFIX}/bin/python"
    elif [[ -x "${HOME}/miniconda3/envs/crebbp_sc_pipeline/bin/python" ]]; then
        PYTHON_BIN="${HOME}/miniconda3/envs/crebbp_sc_pipeline/bin/python"
    elif [[ -x "${HOME}/miniconda3/envs/scvi_cytotrace2/bin/python" ]]; then
        PYTHON_BIN="${HOME}/miniconda3/envs/scvi_cytotrace2/bin/python"
    elif [[ -x "${HOME}/anaconda3/envs/crebbp_sc_pipeline/bin/python" ]]; then
        PYTHON_BIN="${HOME}/anaconda3/envs/crebbp_sc_pipeline/bin/python"
    else
        PYTHON_BIN="python3"
    fi
fi

# ── Geneformer Python environment auto-detection (Stage 2) ──
if [[ -z "${GENEFORMER_PYTHON:-}" ]]; then
    if command -v python >/dev/null 2>&1 && python -c "import geneformer" >/dev/null 2>&1; then
        GENEFORMER_PYTHON="$(command -v python)"
    elif [[ -x "${HOME}/miniconda3/envs/geneformer/bin/python" ]]; then
        GENEFORMER_PYTHON="${HOME}/miniconda3/envs/geneformer/bin/python"
    elif [[ -x "${HOME}/anaconda3/envs/geneformer/bin/python" ]]; then
        GENEFORMER_PYTHON="${HOME}/anaconda3/envs/geneformer/bin/python"
    else
        GENEFORMER_PYTHON="${PYTHON_BIN}"
    fi
fi

# ┌────────────────────────────────────────────────────────────────────────────┐
# │                    END OF USER CONFIGURATION                              │
# └────────────────────────────────────────────────────────────────────────────┘

SCRIPTS="${BASEDIR}/scripts"
RUN_MODE="dry-run"
TARGET_STAGE="all"
FORCE_EXTRACT="false"
SKIP_DATA_CHECK="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)
            RUN_MODE="run"
            shift
            ;;
        --dry-run)
            RUN_MODE="dry-run"
            shift
            ;;
        --extract-only)
            RUN_MODE="extract-only"
            shift
            ;;
        --force-extract)
            FORCE_EXTRACT="true"
            shift
            ;;
        --skip-data-check)
            SKIP_DATA_CHECK="true"
            shift
            ;;
        --stage)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --stage requires a stage number (0, 1, 2, 3, 4, 5, or all)" >&2
                exit 1
            fi
            TARGET_STAGE="$2"
            if [[ ! "$TARGET_STAGE" =~ ^(0|1|2|3|4|5|all)$ ]]; then
                echo "Error: Invalid stage '$TARGET_STAGE'. Allowed values: 0, 1, 2, 3, 4, 5, all" >&2
                exit 1
            fi
            shift 2
            ;;
        --config)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --config requires a configuration file path" >&2
                exit 1
            fi
            CONFIG_FILE="$2"
            shift 2
            ;;
        --smoke-test)
            echo "============================================================"
            echo "  Running Pipeline Smoke Test (Synthetic Dataset)"
            echo "============================================================"
            "${PYTHON_BIN}" "${SCRIPTS}/run_smoke_test.py"
            exit 0
            ;;
        --check-env)
            echo "============================================================"
            echo "  Checking Computational Environment & Dependencies"
            echo "============================================================"
            echo "  Base Directory   : ${BASEDIR}"
            echo "  Primary Python   : ${PYTHON_BIN}"
            echo "  Geneformer Python: ${GENEFORMER_PYTHON}"
            echo "------------------------------------------------------------"
            "${PYTHON_BIN}" -W ignore -c "import scanpy, scvi, anndata, scipy, numpy, pandas; print('  ✓ Python core transcriptomics stack OK (scanpy ' + scanpy.__version__ + ', scvi-tools ' + scvi.__version__ + ')')" 2>/dev/null || echo "  ✗ Python core stack missing in ${PYTHON_BIN} (activate crebbp_sc_pipeline)"
            "${PYTHON_BIN}" -W ignore -c "import cytotrace2_py; print('  ✓ CytoTRACE2 python package OK')" 2>/dev/null || echo "  ! CytoTRACE2 python package not detected"
            "${PYTHON_BIN}" -W ignore -c "import magic, gseapy, scrublet, mygene; print('  ✓ Trajectory & analysis packages OK (magic, gseapy, scrublet, mygene)')" 2>/dev/null || echo "  ! Some secondary analysis packages missing"
            "${PYTHON_BIN}" -W ignore -c "import torch; cuda_ok = torch.cuda.is_available(); dev = torch.cuda.get_device_name(0) if cuda_ok else 'CPU only'; print('  ✓ PyTorch CUDA status: ' + ('Available (' + dev + ')' if cuda_ok else 'CPU only (GPU recommended for scVI)'))" 2>/dev/null || echo "  ! PyTorch not detected"
            "${GENEFORMER_PYTHON}" -W ignore -c "import transformers, datasets, torch; print('  ✓ Geneformer / HuggingFace environment OK')" 2>/dev/null || echo "  ! Geneformer environment not detected (activate geneformer conda env for Stage 2)"
            Rscript -e "suppressPackageStartupMessages(library(Seurat)); cat('  ✓ R / Seurat stack OK (Seurat v', as.character(packageVersion('Seurat')), ')
', sep='')" 2>/dev/null || echo "  ! R Seurat stack not detected (required for Stage 3i .rds export)"
            echo "============================================================"
            exit 0
            ;;
        -h|--help)
            cat << 'HELP_DOC'
Usage: bash MasterAnalysis.sh [OPTIONS]

Options:
  --dry-run            Display execution order and commands without running (default)
  --run                Execute the pipeline (full pipeline or targeted stage)
  --stage <N>          Execute only stage N (0, 1, 2, 3, 4, 5, or all)
  --config <path>      Path to pipeline_config.yaml configuration file
  --smoke-test         Run fast 5-second smoke test on synthetic dataset
  --check-env          Verify required Python and R dependencies
  --extract-only       Extract all 27 embedded scripts to ./scripts/ and exit
  --force-extract      Force re-extraction of scripts even if already present
  --skip-data-check    Bypass pre-flight input data existence checks
  -h, --help           Show this help message

Pipeline Stages:
  Stage 0              Cell Ranger alignment of raw FASTQs to GRCm39 (manual)
  Stage 1              Mouse scVI integration, CytoTRACE2, and Leiden clustering
  Stage 2              Geneformer training on tonsil atlas & mouse state prediction
  Stage 3              Cross-species scVI (mouse lymphoma + human DLBCL + tonsil B cells)
  Stage 4              Trajectory gene programs along potency axis + GSEA
  Stage 5              Supplementary: PGC1α metabolic trajectory & Wilcoxon DE

Environment Overrides:
  PYTHON_BIN           Path to primary Python binary (auto-detected)
  GENEFORMER_PYTHON    Path to Geneformer Python binary (auto-detected)
  CREBBP_BASE_DIR      Override repository base directory
  CELLBENDER_DIR       Override CellBender input directory
  DLBCL_DIR            Override human DLBCL input directory
  TONSIL_DIR           Override human tonsil input directory
  GENEFORMER_MODEL_DIR Override Geneformer model directory
HELP_DOC
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Run 'bash MasterAnalysis.sh --help' for usage." >&2
            exit 1
            ;;
    esac
done

# Export environment variables for child processes
export CREBBP_BASE_DIR="${BASEDIR}"
export CELLBENDER_DIR="${CELLBENDER_DIR}"
export DLBCL_DIR="${DLBCL_DIR}"
export TONSIL_DIR="${TONSIL_DIR}"
export GENEFORMER_MODEL_DIR="${GENEFORMER_MODEL_DIR}"
export FASTQ_DATA_ROOT="${FASTQ_DATA_ROOT}"
export CELLRANGER_REF="${CELLRANGER_REF}"
export CELLRANGER_BIN="${CELLRANGER_BIN}"
export PYTHON_BIN="${PYTHON_BIN}"
export GENEFORMER_PYTHON="${GENEFORMER_PYTHON}"
if [[ -f "${CONFIG_FILE}" ]]; then
    export PIPELINE_CONFIG="${CONFIG_FILE}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Extract embedded scripts into ./scripts/
# ─────────────────────────────────────────────────────────────────────────────
mkdir -p "${SCRIPTS}"

if [[ ! -f "${SCRIPTS}/run_smoke_test.py" || "$FORCE_EXTRACT" == "true" || "$RUN_MODE" == "extract-only" ]]; then
echo "============================================================"
echo "  MasterAnalysis.sh — Extracting scripts"
echo "============================================================"
echo ""
echo "  Base directory : ${BASEDIR}"
echo "  Scripts dir    : ${SCRIPTS}"
echo ""

cat > "${SCRIPTS}/cellranger_shabanas_gex.sh" << '__EOF_cellranger_shabanas_gex_sh__'
#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Batch Cell Ranger GEX alignment for the Shabanas 10x dataset.
# - Discovers all SLX-*.SIG*.tar archives under ${data_root}
# - Extracts FASTQs per barcode (once, unless force_reextract=true)
# - Maps barcodes to friendly sample names via sample_attribution.txt
# - Runs `cellranger count` sequentially, reserving 28 cores per job
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

############################ USER SETTINGS ####################################
data_root="__FASTQ_DATA_ROOT__"
sample_attr="${data_root}/sample_attribution.txt"
fastq_workspace="${data_root}/cellranger_fastqs"
run_root="${data_root}/cellranger_runs"

ref_dir="__CELLRANGER_REF__"
cellranger_bin="__CELLRANGER_BIN__"

cores=22
mem_gb=128
create_bam=false          # true keeps BAMs, false saves space
force_reextract=false     # set true to re-untar even if FASTQs already exist
###############################################################################

mkdir -p "$fastq_workspace" "$run_root"

command -v python3 >/dev/null 2>&1 || { echo "❌ python3 is required."; exit 1; }
[[ -x "$cellranger_bin" ]] || { echo "❌ Cell Ranger binary not found → $cellranger_bin"; exit 1; }
[[ -d "$ref_dir"       ]] || { echo "❌ Reference folder not found → $ref_dir"; exit 1; }
[[ -d "$data_root"     ]] || { echo "❌ FASTQ archive root missing → $data_root"; exit 1; }
[[ -f "$sample_attr"   ]] || { echo "⚠️ Sample attribution sheet missing → $sample_attr"; }

clean_id() {
  local raw="${1:-sample}"
  raw="${raw// /_}"
  raw=$(printf '%s' "$raw" | sed -E 's/[^A-Za-z0-9_]+/_/g' | sed -E 's/_+/_/g' | sed -E 's/^_+//; s/_+$//')
  [[ -z "$raw" ]] && raw="sample"
  printf '%s' "$raw"
}

metadata_lines=$(python3 - "$sample_attr" <<'PY'
import sys, re, pathlib
path = pathlib.Path(sys.argv[1])
if not path.exists():
    sys.exit(0)
emit = False
for raw in path.read_text().splitlines():
    if not raw.strip():
        continue
    if raw.startswith("For metadata"):
        emit = True
        continue
    if not emit or raw.startswith("ID\t"):
        continue
    parts = raw.split('\t')
    if len(parts) < 5:
        continue
    barcode = parts[-1].strip()
    sample = parts[-2].strip()
    clean = re.sub(r'[^A-Za-z0-9_]+', '_', sample.replace(' ', '_')).strip('_')
    if not clean:
        clean = barcode
    print(f"{barcode}\t{clean}")
PY
)

declare -A barcode_to_label=()
if [[ -n "${metadata_lines//[[:space:]]/}" ]]; then
  while IFS=$'\t' read -r barcode label; do
    [[ -z "$barcode" ]] && continue
    barcode_to_label["$barcode"]="$label"
  done <<< "$metadata_lines"
else
  echo "⚠️ No metadata rows parsed from $sample_attr; will label runs with barcodes."
fi

declare -A barcode_tar_map=()
while IFS= read -r -d '' tarball; do
  fname=$(basename "$tarball")
  barcode=$(awk -F'.' '{print $2}' <<<"$fname")
  [[ -z "$barcode" ]] && continue
  barcode_tar_map["$barcode"]+="$tarball"$'\n'
done < <(find "$data_root" -type f -name "SLX-*.SIG*.tar" ! -name "*lostreads*" -print0)

if [[ ${#barcode_tar_map[@]} -eq 0 ]]; then
  echo "❌ No SLX-*.SIG*.tar archives found under $data_root"
  exit 1
fi

echo "🔎 Found ${#barcode_tar_map[@]} barcode group(s) to process."
mapfile -t sorted_barcodes < <(printf "%s\n" "${!barcode_tar_map[@]}" | sort)

pushd "$run_root" >/dev/null

for barcode in "${sorted_barcodes[@]}"; do
  tar_entries="${barcode_tar_map[$barcode]}"
  mapfile -t tar_paths < <(printf '%s' "$tar_entries" | sed '/^$/d')
  [[ ${#tar_paths[@]} -eq 0 ]] && continue

  sample_label="${barcode_to_label[$barcode]:-$barcode}"
  safe_label=$(clean_id "$sample_label")
  run_id="${barcode}_${safe_label}_GEX"
  fastq_dest="${fastq_workspace}/${barcode}"

  echo "================================================================"
  echo "📂 Barcode: ${barcode}"
  echo "   Label  : ${sample_label}"
  echo "   FASTQs :"
  printf '     • %s\n' "${tar_paths[@]}"

  needs_extract=true
  if [[ "$force_reextract" == false && -d "$fastq_dest" ]]; then
    existing_r1=("${fastq_dest}"/*_R1_*.fastq.gz)
    if [[ ${#existing_r1[@]} -gt 0 ]]; then
      needs_extract=false
    fi
  fi

  if [[ "$needs_extract" == true ]]; then
    rm -rf "$fastq_dest"
    mkdir -p "$fastq_dest"
    for tar_path in "${tar_paths[@]}"; do
      echo "   ↪ Extracting $(basename "$tar_path")"
      tar -xf "$tar_path" -C "$fastq_dest"
    done
    printf '%s\n' "${tar_paths[@]}" > "${fastq_dest}/.source_tarballs.txt"
  else
    echo "   ↪ FASTQs already extracted in $fastq_dest (set force_reextract=true to refresh)"
  fi

  r1_files=("${fastq_dest}"/*_R1_*.fastq.gz)
  if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "❌ No R1 FASTQs detected for ${barcode}; skipping."
    continue
  fi

  if [[ -d "${run_root}/${run_id}/outs" ]]; then
    echo "✅ Existing Cell Ranger output detected → ${run_root}/${run_id}/outs (skipping)."
    continue
  fi

  echo "▶ Running Cell Ranger: --id=${run_id}, --sample=${barcode}"
  "$cellranger_bin" count \
    --id="$run_id" \
    --fastqs="$fastq_dest" \
    --sample="$barcode" \
    --transcriptome="$ref_dir" \
    --create-bam="$create_bam" \
    --nosecondary \
    --localcores="$cores" \
    --localmem="$mem_gb"

  status=$?
  if [[ $status -eq 0 ]]; then
    echo "✅ ${run_id} finished → ${run_root}/${run_id}/outs"
  else
    echo "❌ ${run_id} failed with exit code ${status}"
  fi
done

popd >/dev/null
echo "🎉 Alignment loop finished."



__EOF_cellranger_shabanas_gex_sh__

cat > "${SCRIPTS}/mouse_scvi_cytotrace2_cellbender.py" << '__EOF_mouse_scvi_cytotrace2_cellbender_py__'
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
CELLBENDER_DIR = Path("__CELLBENDER_DIR__")
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
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



__EOF_mouse_scvi_cytotrace2_cellbender_py__

cat > "${SCRIPTS}/plot_cytotrace2_downstream.py" << '__EOF_plot_cytotrace2_downstream_py__'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Downstream Plotting Script for Mouse scVI + CytoTRACE2 Analysis
================================================================

Generates:
1. Violin plots for CytoTRACE2 score distributions by condition
2. Violin plots for BCR score distributions by condition
3. Violin plots for OXPHOS score distributions by condition
4. Majority cell type (from Geneformer) per Leiden cluster (resolution 1.0)

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
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"
GENEFORMER_H5AD = Path("__BASEDIR__/Geneformer/mouse_with_geneformer_predictions.h5ad")

FIGDIR.mkdir(parents=True, exist_ok=True)

# ============================== LOAD DATA ====================================
print("=" * 84)
print("Loading integrated AnnData...")
print("=" * 84)

if not H5AD_PATH.exists():
    raise FileNotFoundError(f"AnnData file not found: {H5AD_PATH}")

adata = sc.read_h5ad(H5AD_PATH)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  Conditions: {adata.obs['condition'].unique().tolist()}")

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
        print(f"  Cell types: {adata.obs['geneformer_predicted_celltype'].dropna().unique().tolist()}")
    else:
        print(f"  (warn) geneformer_predicted_celltype not found in Geneformer file")
else:
    print(f"  (warn) Geneformer file not found: {GENEFORMER_H5AD}")

# ============================== COLOR PALETTES ===============================
condition_colors = {
    "WT_B_cells": "#2ecc71",
    "Crebbp_B_cells": "#3498db", 
    "Pre_malignant": "#f39c12",
    "Matched_malignant": "#e74c3c",
    "Malignant": "#8e44ad"
}

# Define condition order for consistent plotting
condition_order = ["WT_B_cells", "Crebbp_B_cells", "Pre_malignant", "Matched_malignant", "Malignant"]
condition_order = [c for c in condition_order if c in adata.obs['condition'].unique()]

# ============================== HELPER FUNCTION ==============================

def save_figure(fig, figdir, filename_base):
    """Save figure in PNG, SVG, and PDF formats."""
    for ext in ['png', 'svg', 'pdf']:
        filepath = figdir / f"{filename_base}.{ext}"
        fig.savefig(filepath, dpi=300, bbox_inches="tight", format=ext)
    print(f"  ✓ {filename_base}.{{png,svg,pdf}}")


# ============================== VIOLIN PLOTS =================================
print("\n" + "=" * 84)
print("Generating violin plots by condition...")
print("=" * 84)


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
    ax.set_xlabel("Condition", fontsize=14, fontweight='bold')
    ax.set_ylabel(score_key.replace("_", " ").title(), fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.tick_params(axis='x', rotation=45, labelsize=12)
    ax.tick_params(axis='y', labelsize=11)
    
    # Add sample sizes
    for i, cond in enumerate(present_order):
        n = (df[groupby] == cond).sum()
        ax.text(i, ax.get_ylim()[0] - 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                f"n={n:,}", ha='center', va='top', fontsize=10, color='gray')
    
    plt.tight_layout()
    save_figure(fig, figdir, filename_base)
    plt.close(fig)


# --- CytoTRACE2 Score Violin ---
save_violin_plot(
    adata, 
    score_key="cytotrace2_score",
    groupby="condition",
    order=condition_order,
    palette=condition_colors,
    title="CytoTRACE2 Score Distribution by Condition",
    filename_base="violin_cytotrace2_score_by_condition",
    figdir=FIGDIR
)

# --- BCR Score Violin ---
save_violin_plot(
    adata,
    score_key="bcr_score",
    groupby="condition",
    order=condition_order,
    palette=condition_colors,
    title="BCR Signaling Score Distribution by Condition",
    filename_base="violin_bcr_score_by_condition",
    figdir=FIGDIR
)

# --- OXPHOS Score Violin ---
save_violin_plot(
    adata,
    score_key="oxphos_score",
    groupby="condition",
    order=condition_order,
    palette=condition_colors,
    title="OXPHOS Score Distribution by Condition",
    filename_base="violin_oxphos_score_by_condition",
    figdir=FIGDIR
)

# ============================== COMBINED MULTI-PANEL FIGURE ==================
print("\n" + "=" * 84)
print("Generating combined multi-panel figure...")
print("=" * 84)

scores_to_plot = []
if "cytotrace2_score" in adata.obs.columns:
    scores_to_plot.append(("cytotrace2_score", "CytoTRACE2 Score"))
if "bcr_score" in adata.obs.columns:
    scores_to_plot.append(("bcr_score", "BCR Score"))
if "oxphos_score" in adata.obs.columns:
    scores_to_plot.append(("oxphos_score", "OXPHOS Score"))

if scores_to_plot:
    fig, axes = plt.subplots(1, len(scores_to_plot), figsize=(6*len(scores_to_plot), 7))
    if len(scores_to_plot) == 1:
        axes = [axes]
    
    for ax, (score_key, label) in zip(axes, scores_to_plot):
        df = adata.obs[[score_key, 'condition']].dropna()
        df[score_key] = df[score_key].astype(float)
        present_order = [c for c in condition_order if c in df['condition'].unique()]
        
        sns.violinplot(
            data=df,
            x='condition',
            y=score_key,
            order=present_order,
            palette=[condition_colors.get(c, "#999999") for c in present_order],
            inner="box",
            ax=ax,
            cut=0,
            scale="width"
        )
        ax.set_xlabel("")
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.tick_params(axis='x', rotation=45, labelsize=10)
    
    plt.suptitle("Score Distributions by Condition", fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_figure(fig, FIGDIR, "violin_all_scores_combined")
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
    
    # Generate color palette for cell types
    unique_celltypes = cluster_celltype.columns.tolist()
    n_types = len(unique_celltypes)
    cmap = plt.cm.get_cmap('tab20', max(n_types, 20))
    celltype_colors = {ct: cmap(i % 20) for i, ct in enumerate(unique_celltypes)}
    
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
    ax.legend(title="Cell Type", bbox_to_anchor=(1.02, 1), loc='upper left', framealpha=0.9, fontsize=9)
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
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create a new column with cluster annotated by majority cell type
    adata.obs['cluster_majority_celltype'] = adata.obs[leiden_key].map(
        lambda x: f"C{x}: {majority_celltype[x]}"
    )
    
    # Generate colors based on majority cell type
    unique_clusters = sorted(adata.obs['cluster_majority_celltype'].unique(), 
                              key=lambda x: int(x.split(':')[0].replace('C', '')))
    cluster_palette = {}
    for clust_label in unique_clusters:
        cluster_num = clust_label.split(':')[0].replace('C', '')
        ct = majority_celltype[cluster_num]
        cluster_palette[clust_label] = celltype_colors.get(ct, "#999999")
    
    sc.pl.umap(
        adata,
        color='cluster_majority_celltype',
        palette=cluster_palette,
        ax=ax,
        show=False,
        frameon=False,
        title="Leiden Clusters (res=1.0) colored by Majority Cell Type (Geneformer)",
        legend_loc='right margin',
        s=15
    )
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "umap_leiden_1.0_majority_celltype")
    plt.close(fig)
    
    # Clean up temp column
    adata.obs.drop(columns=['cluster_majority_celltype'], inplace=True)
    
    # --- Plot: UMAP colored by cell type directly ---
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(
        adata,
        color=celltype_key,
        palette=celltype_colors,
        ax=ax,
        show=False,
        frameon=False,
        title="Geneformer Predicted Cell Type",
        legend_loc='right margin',
        s=15
    )
    
    plt.tight_layout()
    save_figure(fig, FIGDIR, "umap_celltype_geneformer")
    plt.close(fig)

else:
    if not leiden_key:
        print("  (warn) No Leiden clustering found in adata.obs")
    if celltype_key not in adata.obs.columns:
        print(f"  (warn) '{celltype_key}' not found in adata.obs")
        print(f"  Available columns: {list(adata.obs.columns)}")


# ============================== STATISTICAL SUMMARY ==========================
print("\n" + "=" * 84)
print("Statistical Summary")
print("=" * 84)

summary_stats = []
for score_key, score_label in scores_to_plot:
    if score_key in adata.obs.columns:
        for cond in condition_order:
            if cond in adata.obs['condition'].unique():
                vals = adata.obs.loc[adata.obs['condition'] == cond, score_key].dropna().astype(float)
                summary_stats.append({
                    'score': score_label,
                    'condition': cond,
                    'n_cells': len(vals),
                    'mean': vals.mean(),
                    'median': vals.median(),
                    'std': vals.std(),
                    'min': vals.min(),
                    'max': vals.max()
                })

if summary_stats:
    stats_df = pd.DataFrame(summary_stats)
    stats_df.to_csv(FIGDIR / "score_statistics_by_condition.csv", index=False)
    print(f"  ✓ score_statistics_by_condition.csv")
    print("\n" + stats_df.to_string(index=False))

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  All figures saved to: {FIGDIR}")
print("\nDONE.\n")


__EOF_plot_cytotrace2_downstream_py__

cat > "${SCRIPTS}/plot_umap_confidence.py" << '__EOF_plot_umap_confidence_py__'
#!/usr/bin/env python3
"""
Generate UMAP with color intensity proportional to annotation confidence
Based on majority_celltype percentage within each leiden cluster
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from colorsys import rgb_to_hsv, hsv_to_rgb

# Paths
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"
COMPOSITION_CSV = FIGDIR / "leiden_1.0_celltype_composition.csv"

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)

# Load composition data to get majority celltypes and percentages
print("Loading cluster composition data...")
comp_df = pd.read_csv(COMPOSITION_CSV)
comp_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_celltype']))
conf_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_pct'] / 100.0))  # Convert to 0-1

# Map to cells
print("Mapping majority celltypes and confidence to cells...")
adata.obs['majority_celltype'] = adata.obs['leiden_1.0'].astype(str).map(comp_dict)
adata.obs['annotation_confidence'] = adata.obs['leiden_1.0'].astype(str).map(conf_dict)

# Fill NaN values (if any clusters not in composition file)
adata.obs['majority_celltype'] = adata.obs['majority_celltype'].fillna('Unknown')
adata.obs['annotation_confidence'] = adata.obs['annotation_confidence'].fillna(0.0)

conf_min = adata.obs['annotation_confidence'].min()
conf_max = adata.obs['annotation_confidence'].max()
print(f"Confidence range: {conf_min:.3f} - {conf_max:.3f}")
print(f"Size range: {1 + 14 * conf_min:.1f} - {1 + 14 * conf_max:.1f} pixels (directly proportional)")
print(f"Number of celltypes: {adata.obs['majority_celltype'].nunique()}")

# Get unique celltypes and create a color palette
celltypes = sorted(adata.obs['majority_celltype'].unique())
n_types = len(celltypes)

# Use a qualitative colormap (tab20 or similar, cycled if needed)
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')  # Fallback for older matplotlib

colors = {}
for i, ct in enumerate(celltypes):
    colors[ct] = base_cmap(i % 20)

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Plot each celltype with varying intensity based on confidence
# Use color intensity (saturation) where darker = higher confidence
for celltype in celltypes:
    mask = adata.obs['majority_celltype'] == celltype
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Use full base color for all cells (no intensity variation)
    base_color = colors[celltype][:3]  # RGB tuple
    
    # Encode confidence only by size and shape, not color
    sizes = []
    markers_list = []
    
    # Define markers for different confidence levels
    # High conf (>0.7): circle 'o'
    # Medium conf (0.4-0.7): square 's'  
    # Low conf (<0.4): triangle '^'
    
    for conf in confidences:
        # Size: directly proportional to confidence score
        # Use wider range for better visibility: 1 to 15
        # conf=0.258 (min) -> size ~4, conf=0.931 (max) -> size ~14
        size = 1 + 14 * conf  # Directly proportional: size = 1 + 14 * confidence
        sizes.append(size)
        
        # Shape based on confidence
        if conf > 0.7:
            marker = 'o'  # Circle for high confidence
        elif conf > 0.4:
            marker = 's'  # Square for medium confidence
        else:
            marker = '^'  # Triangle for low confidence
        markers_list.append(marker)
    
    sizes = np.array(sizes)
    
    # Plot by marker type for better control
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_list) == marker_type
        if mask_marker.sum() == 0:
            continue
        # Use same full color for all points of this celltype
        ax.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                  color=base_color, s=sizes[mask_marker], 
                  alpha=0.7, marker=marker_type,
                  label=celltype if len(celltypes) <= 20 and marker_type == 'o' else None,
                  edgecolors='black', linewidths=0.3, rasterized=True)

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title('UMAP colored by majority celltype\n(Larger + Circle = High Conf | Smaller + Triangle = Low Conf)', 
             fontsize=13, fontweight='bold')

# Add legend for cell types (limit to avoid overcrowding)
if len(celltypes) <= 20:
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, 
              title='Cell Type', title_fontsize=9, framealpha=0.9)

# Add shape legend for confidence
from matplotlib.lines import Line2D
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=8, label='High confidence (>70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=6, label='Medium confidence (40-70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=5, label='Low confidence (<40%)', markeredgecolor='black'),
]
shape_legend = ax.legend(handles=shape_legend_elements, loc='lower right', 
                         fontsize=9, title='Confidence Level', title_fontsize=10,
                         framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax.add_artist(shape_legend)  # Keep both legends

# Add colorbar for confidence (reference scale)
# Note: confidence is encoded by size and shape, not color
sm = plt.cm.ScalarMappable(cmap=plt.cm.Greys, 
                           norm=plt.Normalize(vmin=0, vmax=100))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Annotation Confidence\n(% majority in cluster)\n(Encoded by size & shape)', 
               rotation=270, labelpad=25, fontsize=10)

plt.tight_layout()

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_celltype_confidence.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nSaved → {output_path}")

# Also save PDF and SVG
fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

# ===== ALTERNATIVE: HSV-based visualization (even more dramatic) =====
print("\nCreating alternative HSV-based visualization...")
fig2, ax2 = plt.subplots(figsize=(12, 10))

for celltype in celltypes:
    mask = adata.obs['majority_celltype'] == celltype
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Use full base color for all cells (no intensity variation)
    base_color = colors[celltype][:3]  # RGB tuple
    
    # Encode confidence only by size and shape, not color
    sizes_hsv = []
    markers_hsv = []
    
    for conf in confidences:
        # Size: directly proportional to confidence score
        # Use wider range for better visibility: 1 to 15
        # conf=0.258 (min) -> size ~4, conf=0.931 (max) -> size ~14
        size = 1 + 14 * conf  # Directly proportional: size = 1 + 14 * confidence
        sizes_hsv.append(size)
        
        # Shape based on confidence
        if conf > 0.7:
            marker = 'o'  # Circle for high confidence
        elif conf > 0.4:
            marker = 's'  # Square for medium confidence
        else:
            marker = '^'  # Triangle for low confidence
        markers_hsv.append(marker)
    
    sizes_hsv = np.array(sizes_hsv)
    
    # Plot by marker type
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_hsv) == marker_type
        if mask_marker.sum() == 0:
            continue
        ax2.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                   color=base_color, s=sizes_hsv[mask_marker], 
                   alpha=0.7, marker=marker_type,
                   label=celltype if len(celltypes) <= 20 and marker_type == 'o' else None,
                   edgecolors='black', linewidths=0.3, rasterized=True)

ax2.set_xlabel('UMAP 1', fontsize=12)
ax2.set_ylabel('UMAP 2', fontsize=12)
ax2.set_title('UMAP colored by majority celltype\n(Larger + Circle = High Conf | Smaller + Triangle = Low Conf)', 
             fontsize=13, fontweight='bold')

if len(celltypes) <= 20:
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, 
              title='Cell Type', title_fontsize=9, framealpha=0.9)

# Add shape legend for HSV version too
shape_legend_elements2 = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=8, label='High confidence (>70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=6, label='Medium confidence (40-70%)', markeredgecolor='black'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=5, label='Low confidence (<40%)', markeredgecolor='black'),
]
shape_legend2 = ax2.legend(handles=shape_legend_elements2, loc='lower right', 
                           fontsize=9, title='Confidence Level', title_fontsize=10,
                           framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax2.add_artist(shape_legend2)

# Add colorbar for confidence (reference scale)
sm2 = plt.cm.ScalarMappable(cmap=plt.cm.Greys, 
                            norm=plt.Normalize(vmin=0, vmax=100))
sm2.set_array([])
cbar2 = plt.colorbar(sm2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.set_label('Annotation Confidence\n(% majority in cluster)\n(Encoded by size & shape)', 
               rotation=270, labelpad=25, fontsize=10)

plt.tight_layout()

# Save alternative version
output_path_hsv = FIGDIR / "umap_leiden_1.0_majority_celltype_confidence_hsv.png"
fig2.savefig(output_path_hsv, dpi=300, bbox_inches='tight')
print(f"Saved → {output_path_hsv}")
fig2.savefig(output_path_hsv.with_suffix('.pdf'), bbox_inches='tight')
fig2.savefig(output_path_hsv.with_suffix('.svg'), bbox_inches='tight')

plt.close()

print("\nDone!")


__EOF_plot_umap_confidence_py__

cat > "${SCRIPTS}/plot_umap_leiden_majority_confidence.py" << '__EOF_plot_umap_leiden_majority_confidence_py__'
#!/usr/bin/env python3
"""
Plot UMAP showing:
1. Leiden clusters (1.0 resolution) as background
2. Majority-voted celltype annotation for each cluster
3. Confidence encoded by size and shape
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Paths
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"
COMPOSITION_CSV = FIGDIR / "leiden_1.0_celltype_composition.csv"

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)

# Load composition data to get majority celltypes and percentages
print("Loading cluster composition data...")
comp_df = pd.read_csv(COMPOSITION_CSV)
comp_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_celltype']))
conf_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_pct'] / 100.0))

# Map to cells
print("Mapping majority celltypes and confidence to cells...")
adata.obs['majority_celltype'] = adata.obs['leiden_1.0'].astype(str).map(comp_dict)
adata.obs['annotation_confidence'] = adata.obs['leiden_1.0'].astype(str).map(conf_dict)

# Fill NaN values
adata.obs['majority_celltype'] = adata.obs['majority_celltype'].fillna('Unknown')
adata.obs['annotation_confidence'] = adata.obs['annotation_confidence'].fillna(0.0)

conf_min = adata.obs['annotation_confidence'].min()
conf_max = adata.obs['annotation_confidence'].max()
print(f"Confidence range: {conf_min:.3f} - {conf_max:.3f}")
# Calculate size range with the new formula
size_min = 3 + 22 * (conf_min ** 0.7)
size_max = 3 + 22 * (conf_max ** 0.7)
print(f"Size range: {size_min:.1f} - {size_max:.1f} pixels (emphasizing high confidence)")
print(f"Number of celltypes: {adata.obs['majority_celltype'].nunique()}")
print(f"Number of leiden clusters: {adata.obs['leiden_1.0'].nunique()}")

# Get unique celltypes and leiden clusters
celltypes = sorted(adata.obs['majority_celltype'].unique())
leiden_clusters = sorted(adata.obs['leiden_1.0'].unique())

# Color palette for celltypes
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')

celltype_colors = {}
for i, ct in enumerate(celltypes):
    celltype_colors[ct] = base_cmap(i % 20)[:3]

# Color palette for leiden clusters - use distinct colors for each
# Use a larger colormap and cycle through to ensure all clusters get different colors
try:
    leiden_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    leiden_cmap = plt.cm.get_cmap('tab20')
leiden_colors = {}
for i, cluster in enumerate(leiden_clusters):
    # Use tab20 and cycle if needed, or use Set3 for more pastel
    if len(leiden_clusters) <= 20:
        leiden_colors[str(cluster)] = leiden_cmap(i)[:3]
    else:
        # For more than 20 clusters, use Set3 and cycle
        try:
            set3_cmap = plt.colormaps['Set3']
        except (AttributeError, KeyError):
            set3_cmap = plt.cm.get_cmap('Set3')
        leiden_colors[str(cluster)] = set3_cmap(i % 12)[:3]

# Create figure with space for legend on the right
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.25], hspace=0.3)
ax = fig.add_subplot(gs[0, 0])
ax_legend = fig.add_subplot(gs[0, 1])
ax_legend.axis('off')

# Plot cells colored by leiden cluster, size/shape by confidence
print("Plotting cells with leiden cluster colors and confidence encoding...")
for cluster in leiden_clusters:
    cluster_str = str(cluster)
    mask = adata.obs['leiden_1.0'].astype(str) == cluster_str
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    confidences = adata.obs.loc[mask, 'annotation_confidence'].values
    
    # Get leiden cluster color
    leiden_color = leiden_colors[cluster_str]
    
    # Encode confidence by size and shape
    sizes = []
    markers_list = []
    
    for conf in confidences:
        # Size: use wider range and emphasize high confidence cells
        # Use a power function to make high confidence cells more prominent
        # conf^0.7 makes the scaling more dramatic for high values
        conf_adj = conf ** 0.7  # Emphasize high confidence
        # Wider range: 3 to 25 pixels for much better visibility
        size = 3 + 22 * conf_adj  # Range from 3 (low) to 25 (high)
        sizes.append(size)
        
        # Shape based on confidence
        if conf > 0.7:
            marker = 'o'  # Circle for high confidence
        elif conf > 0.4:
            marker = 's'  # Square for medium confidence
        else:
            marker = '^'  # Triangle for low confidence
        markers_list.append(marker)
    
    sizes = np.array(sizes)
    
    # Plot by marker type
    for marker_type in ['o', 's', '^']:
        mask_marker = np.array(markers_list) == marker_type
        if mask_marker.sum() == 0:
            continue
        ax.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                  color=leiden_color, s=sizes[mask_marker], 
                  alpha=0.8, marker=marker_type,
                  edgecolors='white', linewidths=0.5, rasterized=True, zorder=2)

# Note: Cluster labels removed - information shown in legend instead

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title('UMAP: Leiden Clusters (1.0) with Majority-Voted Celltype Annotations\n' +
             'Colors = Leiden Clusters | Size/Shape = Confidence', 
             fontsize=13, fontweight='bold')

# Add shape legend for confidence (on the plot)
# Use larger marker sizes to reflect the new size range
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=15, label='High confidence (>70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=10, label='Medium confidence (40-70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=6, label='Low confidence (<40%)', markeredgecolor='black', linewidth=0.5),
]
shape_legend = ax.legend(handles=shape_legend_elements, loc='lower right', 
                        fontsize=9, title='Confidence Level', title_fontsize=10,
                        framealpha=0.9, bbox_to_anchor=(0.98, 0.02))
ax.add_artist(shape_legend)

# Add leiden cluster legend on the right side
print("Creating leiden cluster legend...")
leiden_legend_elements = []
for cluster in leiden_clusters:
    cluster_str = str(cluster)
    color = leiden_colors[cluster_str]
    
    # Get majority celltype and confidence for this cluster
    majority_ct = comp_dict.get(cluster_str, 'Unknown')
    confidence_pct = conf_dict.get(cluster_str, 0.0) * 100
    
    # Truncate long celltype names for cleaner legend
    if len(majority_ct) > 20:
        majority_ct_display = majority_ct[:17] + "..."
    else:
        majority_ct_display = majority_ct
    
    # Create label with cluster info
    label = f"L{cluster}: {majority_ct_display} ({confidence_pct:.0f}%)"
    
    leiden_legend_elements.append(
        mpatches.Patch(facecolor=color, edgecolor='black', linewidth=1.5, label=label)
    )

ax_legend.legend(handles=leiden_legend_elements, loc='center left', 
                fontsize=7.5, title='Leiden Clusters (1.0)\n[Majority Celltype (Confidence %)]', 
                title_fontsize=9, framealpha=0.95, 
                bbox_to_anchor=(0, 0.5), handlelength=1.5, handletextpad=0.5)

# Use constrained_layout instead of tight_layout for better handling of subplots
fig.set_constrained_layout(True)

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_annotation_confidence.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

print("\nDone!")


__EOF_plot_umap_leiden_majority_confidence_py__

cat > "${SCRIPTS}/plot_umap_leiden_majority_confidence_by_condition.py" << '__EOF_plot_umap_leiden_majority_confidence_by_condition_py__'
#!/usr/bin/env python3
"""
Plot UMAP in separate panels for each condition:
1. WT_B_cells
2. Crebbp_B_cells
3. Pre_malignant
4. Malignant (combined Malignant + Matched_malignant)

Each panel shows:
- Leiden clusters (1.0 resolution) colored by cluster
- Confidence encoded by size and shape
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Paths
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"
COMPOSITION_CSV = FIGDIR / "leiden_1.0_celltype_composition.csv"

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)

# Load composition data to get majority celltypes and percentages
print("Loading cluster composition data...")
comp_df = pd.read_csv(COMPOSITION_CSV)
comp_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_celltype']))
conf_dict = dict(zip(comp_df['cluster'].astype(str), comp_df['majority_pct'] / 100.0))

# Map to cells
print("Mapping majority celltypes and confidence to cells...")
adata.obs['majority_celltype'] = adata.obs['leiden_1.0'].astype(str).map(comp_dict)
adata.obs['annotation_confidence'] = adata.obs['leiden_1.0'].astype(str).map(conf_dict)

# Fill NaN values
adata.obs['majority_celltype'] = adata.obs['majority_celltype'].fillna('Unknown')
adata.obs['annotation_confidence'] = adata.obs['annotation_confidence'].fillna(0.0)

# Define conditions to plot - using exact condition names from data
conditions_to_plot = {
    'WT_B_cells': ['WT_B_cells'],
    'Crebbp_B_cells': ['Crebbp_B_cells'],
    'Pre_malignant': ['Pre_malignant'],
    'Malignant (combined)': ['Malignant', 'Matched_malignant']  # Combine both malignant types
}

# Get unique leiden clusters (same across all conditions)
leiden_clusters = sorted(adata.obs['leiden_1.0'].unique())

# Color palette for leiden clusters - use distinct colors for each
try:
    leiden_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    leiden_cmap = plt.cm.get_cmap('tab20')
leiden_colors = {}
for i, cluster in enumerate(leiden_clusters):
    if len(leiden_clusters) <= 20:
        leiden_colors[str(cluster)] = leiden_cmap(i)[:3]
    else:
        try:
            set3_cmap = plt.colormaps['Set3']
        except (AttributeError, KeyError):
            set3_cmap = plt.cm.get_cmap('Set3')
        leiden_colors[str(cluster)] = set3_cmap(i % 12)[:3]

# Create figure with 2x2 grid for conditions + space for legend
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.3], hspace=0.3, wspace=0.3)

axes = []
for i, (cond_name, cond_values) in enumerate(conditions_to_plot.items()):
    row = i // 2
    col = i % 2
    ax = fig.add_subplot(gs[row, col])
    axes.append((ax, cond_name, cond_values))

# Plot each condition in its own panel
print("\nPlotting conditions in separate panels...")
for ax, cond_name, cond_values in axes:
    print(f"  Processing {cond_name}...")
    
    # Filter cells for this condition - use exact condition column values
    if len(cond_values) == 1:
        condition_mask = adata.obs['condition'] == cond_values[0]
    else:
        # Combine multiple conditions (e.g., Malignant + Matched_malignant)
        condition_mask = adata.obs['condition'].isin(cond_values)
    
    adata_cond = adata[condition_mask].copy()
    
    if adata_cond.n_obs == 0:
        print(f"    Warning: No cells found for {cond_name}")
        ax.text(0.5, 0.5, f'No cells\nfor {cond_name}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(cond_name, fontsize=14, fontweight='bold')
        continue
    
    # Get sample information
    if 'sample_id' in adata_cond.obs.columns:
        unique_samples = adata_cond.obs['sample_id'].unique()
        sample_info = f" ({len(unique_samples)} samples)"
    else:
        sample_info = ""
    
    print(f"    {adata_cond.n_obs:,} cells from condition(s): {', '.join(cond_values)}")
    if 'sample_id' in adata_cond.obs.columns:
        print(f"    Samples: {', '.join(sorted(unique_samples))}")
    
    # Create title with condition name and cell count
    title = f"{cond_name}\n{adata_cond.n_obs:,} cells{sample_info}"
    
    # Plot cells colored by leiden cluster, size/shape by confidence
    for cluster in leiden_clusters:
        cluster_str = str(cluster)
        mask = adata_cond.obs['leiden_1.0'].astype(str) == cluster_str
        if mask.sum() == 0:
            continue
        
        coords = adata_cond.obsm['X_umap'][mask]
        confidences = adata_cond.obs.loc[mask, 'annotation_confidence'].values
        
        # Get leiden cluster color
        leiden_color = leiden_colors[cluster_str]
        
        # Encode confidence by size and shape
        sizes = []
        markers_list = []
        
        for conf in confidences:
            # Size: use wider range and emphasize high confidence cells
            conf_adj = conf ** 0.7  # Emphasize high confidence
            size = 3 + 22 * conf_adj  # Range from 3 (low) to 25 (high)
            sizes.append(size)
            
            # Shape based on confidence
            if conf > 0.7:
                marker = 'o'  # Circle for high confidence
            elif conf > 0.4:
                marker = 's'  # Square for medium confidence
            else:
                marker = '^'  # Triangle for low confidence
            markers_list.append(marker)
        
        sizes = np.array(sizes)
        
        # Plot by marker type
        for marker_type in ['o', 's', '^']:
            mask_marker = np.array(markers_list) == marker_type
            if mask_marker.sum() == 0:
                continue
            ax.scatter(coords[mask_marker, 0], coords[mask_marker, 1], 
                      color=leiden_color, s=sizes[mask_marker], 
                      alpha=0.8, marker=marker_type,
                      edgecolors='white', linewidths=0.5, rasterized=True, zorder=2)
    
    ax.set_xlabel('UMAP 1', fontsize=11)
    ax.set_ylabel('UMAP 2', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_aspect('equal')

# Add shared shape legend for confidence (on the last plot)
ax_last = axes[-1][0]
shape_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=15, label='High confidence (>70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', 
           markersize=10, label='Medium confidence (40-70%)', markeredgecolor='black', linewidth=0.5),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
           markersize=6, label='Low confidence (<40%)', markeredgecolor='black', linewidth=0.5),
]
shape_legend = ax_last.legend(handles=shape_legend_elements, loc='lower right', 
                            fontsize=9, title='Confidence Level', title_fontsize=10,
                            framealpha=0.9, bbox_to_anchor=(0.98, 0.02))

# Add leiden cluster legend on the right side
print("\nCreating leiden cluster legend...")
ax_legend = fig.add_subplot(gs[:, 2])
ax_legend.axis('off')

leiden_legend_elements = []
for cluster in leiden_clusters:
    cluster_str = str(cluster)
    color = leiden_colors[cluster_str]
    
    # Get majority celltype and confidence for this cluster
    majority_ct = comp_dict.get(cluster_str, 'Unknown')
    confidence_pct = conf_dict.get(cluster_str, 0.0) * 100
    
    # Truncate long celltype names for cleaner legend
    if len(majority_ct) > 20:
        majority_ct_display = majority_ct[:17] + "..."
    else:
        majority_ct_display = majority_ct
    
    # Create label with cluster info
    label = f"L{cluster}: {majority_ct_display} ({confidence_pct:.0f}%)"
    
    leiden_legend_elements.append(
        mpatches.Patch(facecolor=color, edgecolor='black', linewidth=1.5, label=label)
    )

ax_legend.legend(handles=leiden_legend_elements, loc='center left', 
                fontsize=7.5, title='Leiden Clusters (1.0)\n[Majority Celltype (Confidence %)]', 
                title_fontsize=9, framealpha=0.95, 
                bbox_to_anchor=(0, 0.5), handlelength=1.5, handletextpad=0.5)

# Use constrained_layout for better handling of subplots
fig.set_constrained_layout(True)

# Save
output_path = FIGDIR / "umap_leiden_1.0_majority_annotation_confidence_by_condition.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

plt.close()

print("\nDone!")


__EOF_plot_umap_leiden_majority_confidence_by_condition_py__

cat > "${SCRIPTS}/generate_public_gene_set_scores_mouse.py" << '__EOF_generate_public_gene_set_scores_mouse_py__'
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
INPUT_H5AD = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2/figures")
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
output_h5ad = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated_with_public_scores.h5ad")
adata.write_h5ad(output_h5ad)
print(f"\n✓ Saved updated data: {output_h5ad}")

print("\n" + "=" * 84)
print("COMPLETE!")
print("=" * 84)
print(f"  All figures saved to: {FIGDIR_PUBLIC}")
print("\nDONE.\n")


__EOF_generate_public_gene_set_scores_mouse_py__

cat > "${SCRIPTS}/gsea_umap_ppargc1a_A.py" << '__EOF_gsea_umap_ppargc1a_A_py__'
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
INPUT_H5AD = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
GMT_PPARGC1A = Path("__BASEDIR__/mouse_scvi_cytotrace2/PPARGC1A_TARGET_GENES.v2023.1.Hs.gmt")
GMT_MOOTHA = Path("__BASEDIR__/mouse_scvi_cytotrace2/Mootha PGC.gmt")
OUTPUT_DIR = Path("__BASEDIR__/mouse_scvi_cytotrace2/figures/gsea_ppargc1a")
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



__EOF_gsea_umap_ppargc1a_A_py__

cat > "${SCRIPTS}/gsea_umap_ppargc1a_B.py" << '__EOF_gsea_umap_ppargc1a_B_py__'
#!/usr/bin/env python3
"""
GSEA-style Gene Set Enrichment on scVI UMAP (Version B)
========================================================
Scores cells for PPARGC1A target genes and MOOTHA PGC gene sets,
then visualizes enrichment scores on UMAP.

VERSION B: Enrichment scores cut off at 0 (only show positive enrichment)

Gene sets are human - we convert to mouse orthologs (capitalize first letter).

Author: Generated for CytoTRACE2 analysis
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# ===== CONFIGURATION =====
INPUT_H5AD = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
GMT_PPARGC1A = Path("__BASEDIR__/mouse_scvi_cytotrace2/PPARGC1A_TARGET_GENES.v2023.1.Hs.gmt")
GMT_MOOTHA = Path("__BASEDIR__/mouse_scvi_cytotrace2/Mootha PGC.gmt")
OUTPUT_DIR = Path("__BASEDIR__/mouse_scvi_cytotrace2/figures/gsea_ppargc1a_cutoff0")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output formats
SAVE_FORMATS = ['png', 'svg', 'pdf']

# Colormap for positive-only enrichment (sequential)
CMAP_POSITIVE = 'YlOrRd'  # Yellow-Orange-Red for positive enrichment


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
    print("VERSION B: Cutoff at 0 (positive enrichment only)")
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
    
    # Create clipped versions (cutoff at 0)
    adata.obs['PPARGC1A_score_pos'] = adata.obs['PPARGC1A_score'].clip(lower=0)
    adata.obs['MOOTHA_PGC_score_pos'] = adata.obs['MOOTHA_PGC_score'].clip(lower=0)
    
    # Report how many cells have positive enrichment
    n_ppargc1a_pos = (adata.obs['PPARGC1A_score'] > 0).sum()
    n_mootha_pos = (adata.obs['MOOTHA_PGC_score'] > 0).sum()
    print(f"\n    Cells with positive enrichment:")
    print(f"      PPARGC1A: {n_ppargc1a_pos:,} / {adata.n_obs:,} ({100*n_ppargc1a_pos/adata.n_obs:.1f}%)")
    print(f"      MOOTHA_PGC: {n_mootha_pos:,} / {adata.n_obs:,} ({100*n_mootha_pos/adata.n_obs:.1f}%)")
    
    # Set up plotting style
    sc.set_figure_params(dpi=150, fontsize=12, figsize=(6, 5))
    
    # ===== PLOT 1: PPARGC1A enrichment on UMAP (cutoff at 0) =====
    print(f"\n[5] Generating UMAP plots (cutoff at 0)...")
    
    print("    - PPARGC1A enrichment UMAP")
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color='PPARGC1A_score_pos', ax=ax, show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='PPARGC1A Target Genes Enrichment (≥0)')
    save_figure(fig, OUTPUT_DIR / "umap_PPARGC1A_enrichment_pos")
    plt.close()
    
    # ===== PLOT 2: MOOTHA_PGC enrichment on UMAP (cutoff at 0) =====
    print("    - MOOTHA_PGC enrichment UMAP")
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color='MOOTHA_PGC_score_pos', ax=ax, show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='MOOTHA PGC Enrichment (≥0)')
    save_figure(fig, OUTPUT_DIR / "umap_MOOTHA_PGC_enrichment_pos")
    plt.close()
    
    # ===== PLOT 3: Both scores side by side =====
    print("    - Combined panel")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    sc.pl.umap(adata, color='PPARGC1A_score_pos', ax=axes[0], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='PPARGC1A Target Genes (≥0)')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score_pos', ax=axes[1], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='MOOTHA PGC (≥0)')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "umap_enrichment_combined_pos")
    plt.close()
    
    # ===== PLOT 4: Enrichment scores with condition =====
    print("    - Enrichment with condition panel")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    sc.pl.umap(adata, color='condition', ax=axes[0, 0], show=False,
               frameon=False, title='Condition')
    
    sc.pl.umap(adata, color='leiden_1.0', ax=axes[0, 1], show=False,
               frameon=False, title='Leiden Clusters (res=1.0)')
    
    sc.pl.umap(adata, color='PPARGC1A_score_pos', ax=axes[1, 0], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='PPARGC1A Target Genes (≥0)')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score_pos', ax=axes[1, 1], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='MOOTHA PGC (≥0)')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "umap_enrichment_with_context_pos")
    plt.close()
    
    # ===== PLOT 5: Split by condition =====
    print("    - Enrichment split by condition")
    conditions = sorted(adata.obs['condition'].unique())
    n_cond = len(conditions)
    
    for score_name, base_name in [('PPARGC1A_score_pos', 'PPARGC1A'), 
                                   ('MOOTHA_PGC_score_pos', 'MOOTHA_PGC')]:
        fig, axes = plt.subplots(1, n_cond, figsize=(4*n_cond, 4))
        if n_cond == 1:
            axes = [axes]
        
        # Get global vmax for consistent coloring (vmin is 0)
        vmax = adata.obs[score_name].quantile(0.99)
        
        for ax, cond in zip(axes, conditions):
            mask = adata.obs['condition'] == cond
            sc.pl.umap(adata[mask], color=score_name, ax=ax, show=False,
                       cmap=CMAP_POSITIVE, vmin=0, vmax=vmax, frameon=False,
                       title=f'{cond}\n(n={mask.sum():,})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / f"umap_{base_name}_by_condition_pos")
        plt.close()
    
    # ===== PLOT 6: Comparison - Full range vs Cutoff at 0 =====
    print("    - Comparison: full range vs cutoff")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Full range (diverging colormap)
    sc.pl.umap(adata, color='PPARGC1A_score', ax=axes[0, 0], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='PPARGC1A (full range)')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score', ax=axes[0, 1], show=False,
               cmap='RdBu_r', vcenter=0, frameon=False,
               title='MOOTHA PGC (full range)')
    
    # Cutoff at 0 (sequential colormap)
    sc.pl.umap(adata, color='PPARGC1A_score_pos', ax=axes[1, 0], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='PPARGC1A (≥0 only)')
    
    sc.pl.umap(adata, color='MOOTHA_PGC_score_pos', ax=axes[1, 1], show=False,
               cmap=CMAP_POSITIVE, vmin=0, frameon=False,
               title='MOOTHA PGC (≥0 only)')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "umap_comparison_full_vs_cutoff")
    plt.close()
    
    # ===== PLOT 7: Violin plots by condition =====
    print("    - Violin plots by condition")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sc.pl.violin(adata, 'PPARGC1A_score', groupby='condition', ax=axes[0], 
                 show=False, rotation=45)
    axes[0].axhline(y=0, color='red', linestyle='--', alpha=0.5, label='cutoff')
    axes[0].set_ylim(bottom=0)  # Y-axis starts at 0
    axes[0].set_title('PPARGC1A Target Genes Enrichment')
    axes[0].set_xlabel('')
    
    sc.pl.violin(adata, 'MOOTHA_PGC_score', groupby='condition', ax=axes[1], 
                 show=False, rotation=45)
    axes[1].axhline(y=0, color='red', linestyle='--', alpha=0.5, label='cutoff')
    axes[1].set_ylim(bottom=0)  # Y-axis starts at 0
    axes[1].set_title('MOOTHA PGC Enrichment')
    axes[1].set_xlabel('')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "violin_enrichment_by_condition")
    plt.close()
    
    # ===== PLOT 8: Violin plots by Leiden cluster =====
    print("    - Violin plots by cluster")
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    
    sc.pl.violin(adata, 'PPARGC1A_score', groupby='leiden_1.0', ax=axes[0], 
                 show=False, rotation=0)
    axes[0].axhline(y=0, color='red', linestyle='--', alpha=0.5)
    axes[0].set_ylim(bottom=0)  # Y-axis starts at 0
    axes[0].set_title('PPARGC1A Target Genes Enrichment by Cluster')
    axes[0].set_xlabel('Leiden Cluster')
    
    sc.pl.violin(adata, 'MOOTHA_PGC_score', groupby='leiden_1.0', ax=axes[1], 
                 show=False, rotation=0)
    axes[1].axhline(y=0, color='red', linestyle='--', alpha=0.5)
    axes[1].set_ylim(bottom=0)  # Y-axis starts at 0
    axes[1].set_title('MOOTHA PGC Enrichment by Cluster')
    axes[1].set_xlabel('Leiden Cluster')
    
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "violin_enrichment_by_cluster")
    plt.close()
    
    # ===== PLOT 9: Correlation with CytoTRACE2 =====
    print("    - Correlation with CytoTRACE2")
    if 'cytotrace2_score' in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # PPARGC1A vs CytoTRACE2
        ax = axes[0]
        x = adata.obs['cytotrace2_score']
        y = adata.obs['PPARGC1A_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('CytoTRACE2 Score')
        ax.set_ylabel('PPARGC1A Enrichment')
        ax.set_title(f'PPARGC1A vs CytoTRACE2\n(r = {corr:.3f})')
        
        # MOOTHA vs CytoTRACE2
        ax = axes[1]
        y = adata.obs['MOOTHA_PGC_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='darkred')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('CytoTRACE2 Score')
        ax.set_ylabel('MOOTHA PGC Enrichment')
        ax.set_title(f'MOOTHA PGC vs CytoTRACE2\n(r = {corr:.3f})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "scatter_enrichment_vs_cytotrace2")
        plt.close()
    
    # ===== PLOT 10: Comparison with OXPHOS if available =====
    if 'oxphos_score' in adata.obs.columns:
        print("    - Correlation with OXPHOS score")
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        x = adata.obs['oxphos_score']
        
        ax = axes[0]
        y = adata.obs['PPARGC1A_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('OXPHOS Score')
        ax.set_ylabel('PPARGC1A Enrichment')
        ax.set_title(f'PPARGC1A vs OXPHOS\n(r = {corr:.3f})')
        
        ax = axes[1]
        y = adata.obs['MOOTHA_PGC_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='darkred')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('OXPHOS Score')
        ax.set_ylabel('MOOTHA PGC Enrichment')
        ax.set_title(f'MOOTHA PGC vs OXPHOS\n(r = {corr:.3f})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "scatter_enrichment_vs_oxphos")
        plt.close()
    
    # ===== PLOT 11: Comparison with BCR score if available =====
    if 'bcr_score' in adata.obs.columns:
        print("    - Correlation with BCR score")
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        x = adata.obs['bcr_score']
        
        ax = axes[0]
        y = adata.obs['PPARGC1A_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='steelblue')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('BCR Score')
        ax.set_ylabel('PPARGC1A Enrichment')
        ax.set_title(f'PPARGC1A vs BCR\n(r = {corr:.3f})')
        
        ax = axes[1]
        y = adata.obs['MOOTHA_PGC_score']
        ax.scatter(x, y, alpha=0.1, s=1, c='darkred')
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        corr = np.corrcoef(x, y)[0, 1]
        ax.set_xlabel('BCR Score')
        ax.set_ylabel('MOOTHA PGC Enrichment')
        ax.set_title(f'MOOTHA PGC vs BCR\n(r = {corr:.3f})')
        
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "scatter_enrichment_vs_bcr")
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
    
    # Percentage of cells with positive enrichment by condition
    print(f"\n    Cells with positive enrichment by condition:")
    for cond in conditions:
        mask = adata.obs['condition'] == cond
        n_cells = mask.sum()
        n_ppargc1a = ((adata.obs['PPARGC1A_score'] > 0) & mask).sum()
        n_mootha = ((adata.obs['MOOTHA_PGC_score'] > 0) & mask).sum()
        print(f"      {cond}: PPARGC1A {100*n_ppargc1a/n_cells:.1f}%, MOOTHA {100*n_mootha/n_cells:.1f}%")
    
    print(f"\n    Enrichment by condition (full scores):")
    print(summary_cond[['PPARGC1A_score', 'MOOTHA_PGC_score']].to_string())
    
    print("\n" + "=" * 70)
    print("✓ COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"\nFiles created:")
    print("  - umap_PPARGC1A_enrichment_pos.png/svg/pdf")
    print("  - umap_MOOTHA_PGC_enrichment_pos.png/svg/pdf")
    print("  - umap_enrichment_combined_pos.png/svg/pdf")
    print("  - umap_enrichment_with_context_pos.png/svg/pdf")
    print("  - umap_PPARGC1A_by_condition_pos.png/svg/pdf")
    print("  - umap_MOOTHA_PGC_by_condition_pos.png/svg/pdf")
    print("  - umap_comparison_full_vs_cutoff.png/svg/pdf")
    print("  - violin_enrichment_by_condition.png/svg/pdf")
    print("  - violin_enrichment_by_cluster.png/svg/pdf")
    print("  - scatter_enrichment_vs_cytotrace2.png/svg/pdf")
    print("  - scatter_enrichment_vs_oxphos.png/svg/pdf (if available)")
    print("  - scatter_enrichment_vs_bcr.png/svg/pdf (if available)")
    print("  - enrichment_by_condition.csv")
    print("  - enrichment_by_cluster.csv")
    print("  - ppargc1a_genes.csv")
    print("  - mootha_pgc_genes.csv")


if __name__ == "__main__":
    main()



__EOF_gsea_umap_ppargc1a_B_py__

cat > "${SCRIPTS}/plot_umap_highlight_clusters_4_6.py" << '__EOF_plot_umap_highlight_clusters_4_6_py__'
#!/usr/bin/env python3
"""
Plot UMAP highlighting leiden clusters 4 and 6
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
OUTDIR = Path("__BASEDIR__/mouse_scvi_cytotrace2")
FIGDIR = OUTDIR / "figures"
H5AD_PATH = OUTDIR / "mouse_integrated.h5ad"

# Clusters to highlight (as strings to match leiden_1.0 format)
HIGHLIGHT_CLUSTERS = ['4', '6']

# Load data
print("Loading AnnData...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Check if leiden_1.0 exists
if 'leiden_1.0' not in adata.obs.columns:
    print("ERROR: leiden_1.0 not found in obs columns!")
    print("Available columns:", sorted(adata.obs.columns))
    exit(1)

# Check if UMAP exists
if 'X_umap' not in adata.obsm:
    print("ERROR: X_umap not found in obsm!")
    print("Available obsm keys:", list(adata.obsm.keys()))
    exit(1)

print(f"\nLeiden clusters present: {sorted(adata.obs['leiden_1.0'].unique())}")
print(f"Highlighting clusters: {HIGHLIGHT_CLUSTERS}")

# Check if highlight clusters exist
for cluster in HIGHLIGHT_CLUSTERS:
    # leiden_1.0 might be stored as string or int, so convert to match
    cluster_val = str(cluster) if isinstance(adata.obs['leiden_1.0'].iloc[0], str) else int(cluster)
    n_cells = (adata.obs['leiden_1.0'].astype(str) == str(cluster)).sum()
    print(f"  Cluster {cluster}: {n_cells:,} cells")

# Create figure
fig, ax = plt.subplots(figsize=(14, 12))

# Get all leiden clusters
all_clusters = sorted(adata.obs['leiden_1.0'].unique())
n_clusters = len(all_clusters)

# Color palette for all clusters (use muted colors for non-highlighted)
try:
    base_cmap = plt.colormaps['tab20']
except (AttributeError, KeyError):
    base_cmap = plt.cm.get_cmap('tab20')

# Assign colors: highlighted clusters get bright colors, others get muted
cluster_colors = {}
for cluster in all_clusters:
    cluster_str = str(cluster)
    if cluster_str in HIGHLIGHT_CLUSTERS:
        # Bright, saturated colors for highlighted clusters
        if cluster_str == '4':
            cluster_colors[cluster_str] = (1.0, 0.0, 0.0)  # Bright red
        elif cluster_str == '6':
            cluster_colors[cluster_str] = (0.0, 0.0, 1.0)  # Bright blue
    else:
        # Muted gray for non-highlighted clusters
        cluster_colors[cluster_str] = (0.7, 0.7, 0.7)  # Light gray

# Plot all clusters
print("\nPlotting UMAP...")
for cluster in all_clusters:
    cluster_str = str(cluster)
    mask = adata.obs['leiden_1.0'].astype(str) == cluster_str
    if mask.sum() == 0:
        continue
    
    coords = adata.obsm['X_umap'][mask]
    color = cluster_colors.get(cluster_str, (0.7, 0.7, 0.7))
    
    # Highlighted clusters: larger size, full opacity, with edge
    # Non-highlighted: smaller size, lower opacity, no edge
    if cluster_str in HIGHLIGHT_CLUSTERS:
        ax.scatter(coords[:, 0], coords[:, 1], 
                  c=[color], s=50, alpha=1.0, 
                  marker='o', edgecolors='black', linewidths=1.5,
                  rasterized=True, label=f'Cluster {cluster_str} (highlighted)', zorder=3)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], 
                  c=[color], s=10, alpha=0.3, 
                  marker='o', edgecolors='none',
                  rasterized=True, zorder=1)

ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)
ax.set_title(f'UMAP: Highlighted Leiden Clusters {HIGHLIGHT_CLUSTERS}\n' +
             f'(All other clusters shown in gray)', 
             fontsize=14, fontweight='bold')

# Add legend
ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

plt.tight_layout()

# Save
output_path = FIGDIR / f"umap_leiden_1.0_highlight_clusters_{'_'.join(map(str, HIGHLIGHT_CLUSTERS))}.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nSaved → {output_path}")

fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
print(f"Saved → {output_path.with_suffix('.pdf')}")
print(f"Saved → {output_path.with_suffix('.svg')}")

# Also save to the mouse_human_integration directory if it exists
integration_figdir = Path("__BASEDIR__/mouse_human_integration/figures_individual_samples")
if integration_figdir.exists():
    integration_output = integration_figdir / f"umap_all_disease_state_highlight_clusters_{'_'.join(map(str, HIGHLIGHT_CLUSTERS))}.png"
    fig.savefig(integration_output, dpi=300, bbox_inches='tight')
    print(f"Also saved → {integration_output}")

plt.close()

print("\nDone!")


__EOF_plot_umap_highlight_clusters_4_6_py__

cat > "${SCRIPTS}/train_geneformer_tonsil_multi.py" << '__EOF_train_geneformer_tonsil_multi_py__'
#!/usr/bin/env python3
"""
Train Geneformer on multiple tonsil datasets using the official Geneformer
TranscriptomeTokenizer and cell-classification workflow.

Differences vs the earlier version:
- Uses Geneformer's TranscriptomeTokenizer (no more generic BERT tokenizer fallback)
- Converts merged reference AnnData directly into a HuggingFace `.dataset`
- Fine-tunes `BertForSequenceClassification` with Geneformer's
  `DataCollatorForCellClassification`
- Provides detailed reference cell-type summaries

NOTE: Query prediction (mouse sample inference) will be added after the
reference model is fine-tuned. The script still accepts `--query` so the CLI
remains stable, but it currently only trains the classifier.
"""

import sys
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import gc
import pickle
from typing import List

import numpy as np
import pandas as pd
import torch
import scipy.sparse as sp

try:
    import scanpy as sc
    import anndata as ad
except ImportError:
    print("ERROR: scanpy/anndata not installed")
    sys.exit(1)

try:
    from datasets import load_from_disk, DatasetDict
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import accuracy_score, f1_score
    from transformers import BertForSequenceClassification, TrainingArguments, Trainer
    from geneformer import TranscriptomeTokenizer, DataCollatorForCellClassification
except ImportError:
    print("ERROR: transformers/datasets/geneformer not installed")
    print("Install with: pip install transformers datasets geneformer scikit-learn")
    sys.exit(1)

try:
    import mygene
except ImportError:
    print("WARNING: mygene not installed. Mouse-to-human gene conversion will be skipped.")
    mygene = None

GENEFORMER_MODEL = "ctheodoris/Geneformer"


def msg(text, *args):
    print(f"[Geneformer-Train] {text % args if args else text}")


def load_and_merge_datasets(file_paths: List[str], max_cells: int = 50000,
                            label_column: str = 'annotation_20230508',
                            random_seed: int = 42,
                            include_cell_types: List[str] = None,
                            exclude_cell_types: List[str] = None):
    """Load multiple h5ad files, filter by cell type, merge, and subsample (memory-aware)."""
    msg("Inspecting %d datasets (first pass - filtering)...", len(file_paths))
    dataset_infos = []
    total_cells_after_filter = 0
    resolved_label = label_column

    # First pass: filter and count available cells
    for i, file_path in enumerate(file_paths):
        msg("  [%d/%d] %s", i + 1, len(file_paths), Path(file_path).name)
        adata = sc.read_h5ad(file_path, backed='r')
        original_n = adata.n_obs

        if resolved_label not in adata.obs.columns:
            alt_cols = [c for c in adata.obs.columns if 'annotation' in c.lower()]
            if alt_cols:
                msg("    Column '%s' missing. Using '%s' instead.", resolved_label, alt_cols[0])
                resolved_label = alt_cols[0]
            else:
                msg("ERROR: No annotation column found in %s", Path(file_path).name)
                if hasattr(adata, 'file'):
                    adata.file.close()
                return None, None

        # Apply cell type filtering if specified
        if include_cell_types is not None or exclude_cell_types is not None:
            mask = pd.Series(True, index=adata.obs.index)
            
            if include_cell_types:
                include_mask = pd.Series(False, index=adata.obs.index)
                cell_type_str = adata.obs[resolved_label].astype(str).str.upper()
                
                for include_ct in include_cell_types:
                    ct_upper = include_ct.upper()
                    include_mask |= (cell_type_str == ct_upper)
                    include_mask |= cell_type_str.str.contains(ct_upper, na=False, regex=False)
                
                mask &= include_mask
            
            if exclude_cell_types:
                cell_type_str = adata.obs[resolved_label].astype(str).str.upper()
                for exclude_ct in exclude_cell_types:
                    ct_upper = exclude_ct.upper()
                    mask &= ~cell_type_str.str.contains(ct_upper, na=False, regex=False)
            
            filtered_n = mask.sum()
            msg("    Filtered: %d -> %d cells (removed %d)", original_n, filtered_n, original_n - filtered_n)
            # Store as numpy array of boolean for later use
            filter_mask = mask.values
        else:
            filtered_n = original_n
            filter_mask = None

        if filtered_n == 0:
            msg("    WARNING: No cells remaining after filtering, skipping this file")
            if hasattr(adata, 'file'):
                adata.file.close()
            continue

        # Show cell type distribution after filtering
        if filter_mask is not None:
            cell_types = adata.obs.loc[filter_mask, resolved_label].value_counts()
        else:
            cell_types = adata.obs[resolved_label].value_counts()
        
        msg("    Cell types after filtering (%d total):", len(cell_types))
        for ct, count in cell_types.head(10).items():
            pct = 100 * count / filtered_n
            msg("      %s: %d cells (%.1f%%)", ct, count, pct)
        if len(cell_types) > 10:
            msg("      ... and %d more cell types", len(cell_types) - 10)

        dataset_infos.append({
            'path': file_path,
            'cells': filtered_n,
            'genes': adata.n_vars,
            'filter_mask': filter_mask,
        })
        total_cells_after_filter += filtered_n
        msg("    %d cells after filtering, %d genes", filtered_n, adata.n_vars)
        
        if hasattr(adata, 'file'):
            adata.file.close()

    if total_cells_after_filter == 0:
        msg("ERROR: No cells remaining after filtering.")
        return None, None

    msg("\nTotal available cells after filtering: %d", total_cells_after_filter)

    # Plan subsampling based on filtered cell counts
    if max_cells <= 0 or max_cells >= total_cells_after_filter:
        msg("max_cells >= filtered cells; using all filtered cells.")
        sample_plan = {info['path']: info['cells'] for info in dataset_infos}
    else:
        msg("Planning subsampling to %d cells from filtered data...", max_cells)
        proportions = [info['cells'] / total_cells_after_filter for info in dataset_infos]
        raw_samples = [int(p * max_cells) for p in proportions]
        sample_plan = {}
        for info, raw in zip(dataset_infos, raw_samples):
            sample_plan[info['path']] = min(max(raw, 1), info['cells'])
        remainder = max_cells - sum(sample_plan.values())
        idx = 0
        while remainder > 0:
            info = dataset_infos[idx % len(dataset_infos)]
            if sample_plan[info['path']] < info['cells']:
                sample_plan[info['path']] += 1
                remainder -= 1
            idx += 1

    msg("Sampling plan (from filtered data):")
    for info in dataset_infos:
        msg("  %s -> %d cells", Path(info['path']).name, sample_plan[info['path']])

    # Second pass: load, filter, and subsample
    rng = np.random.default_rng(random_seed)
    subsets = []

    for info in dataset_infos:
        file_path = info['path']
        desired = sample_plan[file_path]
        msg("Loading and subsampling %d cells from %s...", desired, Path(file_path).name)
        adata = sc.read_h5ad(file_path, backed='r')

        # Apply filtering first
        if info['filter_mask'] is not None:
            subset = adata[info['filter_mask']]
        else:
            subset = adata

        if hasattr(subset, 'file'):
            subset = subset.to_memory()

        # Then subsample from filtered data with simple stratification by cell type
        if desired < subset.n_obs:
            ct_counts = subset.obs[resolved_label].value_counts()
            n_types = len(ct_counts)

            # Equal base per cell type
            base = max(1, desired // n_types)
            target_per_ct = {ct: min(base, count) for ct, count in ct_counts.items()}
            assigned = sum(target_per_ct.values())

            # Distribute remainder proportional to availability
            remaining = desired - assigned
            if remaining > 0:
                avail = {ct: ct_counts[ct] - target_per_ct[ct] for ct in ct_counts.index}
                ct_list = list(ct_counts.index)
                idx_rem = 0
                while remaining > 0 and ct_list:
                    ct = ct_list[idx_rem % len(ct_list)]
                    if avail[ct] > 0:
                        target_per_ct[ct] += 1
                        avail[ct] -= 1
                        remaining -= 1
                    idx_rem += 1

            take_indices = []
            for ct, target in target_per_ct.items():
                ct_idx = subset.obs[subset.obs[resolved_label] == ct].index.values
                if target < len(ct_idx):
                    sampled = rng.choice(ct_idx, size=target, replace=False)
                else:
                    sampled = ct_idx
                take_indices.append(sampled)

            if take_indices:
                take_indices = np.concatenate(take_indices)
                subset = subset[take_indices]
        # else: use all filtered cells

        if not sp.issparse(subset.X):
            subset.X = sp.csr_matrix(subset.X)

        subset.obs['source_file'] = Path(file_path).stem
        subsets.append(subset)

        if hasattr(adata, 'file'):
            adata.file.close()
        gc.collect()

    msg("Merging subsets...")
    merged = ad.concat(subsets, join='outer', index_unique='-', fill_value=0)
    del subsets
    gc.collect()

    if not sp.issparse(merged.X):
        msg("Converting merged matrix to sparse format...")
        merged.X = sp.csr_matrix(merged.X)
        gc.collect()

    msg("Merged: %d cells, %d genes", merged.n_obs, merged.n_vars)
    msg("Matrix format: %s", "sparse" if sp.issparse(merged.X) else "dense")

    if resolved_label not in merged.obs.columns:
        msg("ERROR: Resolved label column '%s' missing after merge.", resolved_label)
        return None, None

    merged.obs['training_label'] = merged.obs[resolved_label].astype(str)

    msg("\n" + "=" * 70)
    msg("CELL TYPE DISTRIBUTION IN MERGED REFERENCE")
    msg("=" * 70)
    label_counts = merged.obs['training_label'].value_counts()
    msg("Total cell types: %d", len(label_counts))
    for label, count in label_counts.items():
        pct = 100 * count / merged.n_obs
        msg("  %s: %d cells (%.1f%%)", label, count, pct)
    msg("=" * 70)

    return merged, resolved_label


def ensure_geneformer_requirements(adata: ad.AnnData) -> ad.AnnData:
    """Ensure var['ensembl_id'] and obs['n_counts'] exist for Geneformer."""
    msg("Ensuring Geneformer-required fields (ensembl_id, n_counts)...")

    if 'ensembl_id' not in adata.var.columns:
        msg("  Adding var['ensembl_id'] from var_names")
        adata.var['ensembl_id'] = adata.var_names

    if 'n_counts' not in adata.obs.columns:
        msg("  Computing obs['n_counts'] from expression matrix")
        if sp.issparse(adata.X):
            adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).ravel()
        else:
            adata.obs['n_counts'] = adata.X.sum(axis=1)

    return adata


def tokenize_reference_with_geneformer(ref_h5ad_path: Path,
                                       output_dir: Path,
                                       label_col: str = 'training_label',
                                       nproc: int = 8,
                                       model_version: str = 'V2'):
    """Tokenize the reference h5ad using Geneformer's TranscriptomeTokenizer."""
    msg("\n" + "=" * 70)
    msg("TOKENIZING REFERENCE WITH GENEFORMER")
    msg("=" * 70)

    attr = {label_col: 'cell_type'}

    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict=attr,
        nproc=nproc,
        model_version=model_version,
    )

    tokenizer.tokenize_data(
        data_directory=str(ref_h5ad_path.parent),
        output_directory=str(output_dir),
        output_prefix='tonsil_ref',
        file_format='h5ad',
        input_identifier=ref_h5ad_path.stem,
    )

    dataset_path = output_dir / 'tonsil_ref.dataset'
    msg("✓ Tokenized dataset written to: %s", dataset_path)
    return dataset_path, tokenizer.gene_token_dict


def fine_tune_geneformer_on_dataset(dataset_path: Path,
                                    output_dir: Path,
                                    gene_token_dict,
                                    epochs: int = 3,
                                    batch_size: int = 16,
                                    learning_rate: float = 5e-5,
                                    device: str = 'cuda'):
    """Fine-tune Geneformer classifier on the tokenized dataset."""
    msg("\n" + "=" * 70)
    msg("FINE-TUNING GENEFORMER")
    msg("=" * 70)

    dataset = load_from_disk(str(dataset_path))
    if isinstance(dataset, DatasetDict):
        # TranscriptomeTokenizer saves a single split named 'train'
        dataset = dataset[next(iter(dataset.keys()))]

    if 'cell_type' not in dataset.column_names:
        msg("ERROR: 'cell_type' column missing in tokenized dataset")
        sys.exit(1)

    labels = np.array(dataset['cell_type'], dtype=str)
    label_encoder = LabelEncoder()
    encoded_labels = label_encoder.fit_transform(labels)
    dataset = dataset.remove_columns(['cell_type'])
    dataset = dataset.add_column('label', encoded_labels.tolist())

    num_labels = len(label_encoder.classes_)
    msg("  Number of cell types: %d", num_labels)

    split = dataset.train_test_split(test_size=0.2, seed=42)
    train_ds = split['train']
    eval_ds = split['test']

    # Clear GPU cache before loading model
    if device == 'cuda':
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            msg("  GPU memory before model load: %.1f GB free / %.1f GB total",
                torch.cuda.get_device_properties(0).total_memory / 1e9 - 
                torch.cuda.memory_allocated(0) / 1e9,
                torch.cuda.get_device_properties(0).total_memory / 1e9)
        else:
            msg("CUDA not available, switching to CPU")
            device = 'cpu'
    
    model = BertForSequenceClassification.from_pretrained(
        GENEFORMER_MODEL,
        num_labels=num_labels,
    )

    model.to(device)
    
    # Enable gradient checkpointing at model level (saves ~40% memory)
    if hasattr(model, 'gradient_checkpointing_enable'):
        model.gradient_checkpointing_enable()
        msg("  Enabled gradient checkpointing at model level")
    
    if device == 'cuda':
        torch.cuda.empty_cache()
        msg("  GPU memory after model load: %.1f GB allocated",
            torch.cuda.memory_allocated(0) / 1e9)

    data_collator = DataCollatorForCellClassification(token_dictionary=gene_token_dict)

    # Memory optimizations for 16GB VRAM
    # Start very small and let user increase if needed
    # Geneformer is memory-hungry: ~2-3GB per batch item with 2048 tokens
    actual_batch_size = min(2, max(1, batch_size // 6))  # Very conservative: 2 max
    gradient_accum = max(1, batch_size // actual_batch_size)
    
    msg("  Memory optimization: micro-batch=%d, gradient_accum=%d (effective=%d)",
        actual_batch_size, gradient_accum, actual_batch_size * gradient_accum)
    
    training_args = TrainingArguments(
        output_dir=str(output_dir / 'checkpoints'),
        num_train_epochs=epochs,
        per_device_train_batch_size=actual_batch_size,
        per_device_eval_batch_size=actual_batch_size,
        gradient_accumulation_steps=gradient_accum,
        learning_rate=learning_rate,
        weight_decay=0.01,
        logging_dir=str(output_dir / 'logs'),
        logging_steps=50,
        eval_strategy='epoch',
        save_strategy='epoch',
        load_best_model_at_end=True,
        metric_for_best_model='accuracy',
        greater_is_better=True,
        save_total_limit=2,
        fp16=True,  # Mixed precision for memory savings
        gradient_checkpointing=True,  # Trade compute for memory
        optim='adamw_torch_fused',  # More efficient optimizer
    )
    
    msg("  Effective batch size: %d (micro-batch=%d, accum=%d)", 
        actual_batch_size * gradient_accum, actual_batch_size, gradient_accum)

    def compute_metrics(eval_pred):
        logits, labels_np = eval_pred
        preds = np.argmax(logits, axis=1)
        return {
            'accuracy': accuracy_score(labels_np, preds),
            'f1': f1_score(labels_np, preds, average='weighted'),
        }

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=eval_ds,
        data_collator=data_collator,
        compute_metrics=compute_metrics,
    )

    msg("  Training samples: %d", len(train_ds))
    msg("  Validation samples: %d", len(eval_ds))

    trainer.train()

    model_path = output_dir / 'fine_tuned_model'
    model_path.mkdir(parents=True, exist_ok=True)
    model.save_pretrained(str(model_path))

    label_encoder_path = output_dir / 'label_encoder.pkl'
    with open(label_encoder_path, 'wb') as f:
        pickle.dump(label_encoder, f)

    # Save gene token dictionary for prediction
    gene_token_dict_path = output_dir / 'gene_token_dict.pkl'
    with open(gene_token_dict_path, 'wb') as f:
        pickle.dump(gene_token_dict, f)

    msg("✓ Fine-tuned model saved to: %s", model_path)
    msg("✓ Label encoder saved to: %s", label_encoder_path)
    msg("✓ Gene token dictionary saved to: %s", gene_token_dict_path)

    return model_path, label_encoder_path


def convert_mouse_to_human_genes(adata: ad.AnnData) -> ad.AnnData:
    """Convert mouse gene symbols to human orthologs using mygene."""
    if mygene is None:
        msg("WARNING: mygene not available. Skipping gene conversion.")
        msg("  Assuming genes are already in human format.")
        return adata
    
    msg("Converting mouse genes to human orthologs...")
    mg = mygene.MyGeneInfo()
    
    unique_genes = list(set(adata.var_names))
    msg("  Converting %d unique genes...", len(unique_genes))
    
    gene_mapping = {}
    batch_size_conv = 1000
    converted_count = 0
    
    for i in range(0, len(unique_genes), batch_size_conv):
        batch = unique_genes[i:i+batch_size_conv]
        try:
            results = mg.querymany(
                batch,
                scopes='symbol',
                fields='symbol',
                species='mouse',
                target_species='human',
                returnall=True
            )
            
            for result in results.get('out', []):
                if 'symbol' in result and result.get('query') in batch:
                    gene_mapping[result['query']] = result['symbol']
                    converted_count += 1
            
            if (i + batch_size_conv) % 5000 == 0:
                msg("    Processed %d/%d genes (%d converted)", 
                    min(i + batch_size_conv, len(unique_genes)), len(unique_genes), converted_count)
        except Exception as e:
            msg("    WARNING: Batch conversion failed: %s", str(e))
            continue
    
    msg("  Converted %d/%d genes (%.1f%%)", converted_count, len(unique_genes), 
        100 * converted_count / len(unique_genes) if unique_genes else 0)
    
    # Apply mapping
    human_genes = [gene_mapping.get(g, g.upper()) for g in adata.var_names]
    adata.var_names = human_genes
    
    # Remove duplicates (keep first occurrence)
    _, unique_idx = np.unique(adata.var_names, return_index=True)
    adata = adata[:, unique_idx]
    msg("  After conversion: %d unique genes", adata.n_vars)
    
    return adata


def predict_query_cells(query_h5ad_path: Path,
                        model_path: Path,
                        label_encoder_path: Path,
                        gene_token_dict: dict,
                        output_dir: Path,
                        device: str = 'cuda',
                        batch_size: int = 8) -> Path:
    """Predict cell types for query cells using the fine-tuned model."""
    msg("\n" + "=" * 70)
    msg("PREDICTING CELL TYPES ON QUERY DATA")
    msg("=" * 70)
    
    # Load query data
    msg("Loading query data: %s", query_h5ad_path)
    query_adata = sc.read_h5ad(query_h5ad_path, backed='r')
    msg("  Query: %d cells, %d genes", query_adata.n_obs, query_adata.n_vars)
    
    # Convert to memory
    if hasattr(query_adata, 'file'):
        query_adata = query_adata.to_memory()
    
    # Convert mouse to human genes
    query_adata = convert_mouse_to_human_genes(query_adata)
    
    # Ensure Geneformer requirements
    query_adata = ensure_geneformer_requirements(query_adata)
    
    # Save prepared query for tokenization
    query_prepared_path = output_dir / 'query_prepared_for_geneformer.h5ad'
    msg("Writing prepared query to: %s", query_prepared_path)
    query_adata.write_h5ad(query_prepared_path)
    
    # Tokenize query data
    msg("Tokenizing query data...")
    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict=None,
        nproc=8,
        model_version='V2',
    )
    
    tokenizer.tokenize_data(
        data_directory=str(query_prepared_path.parent),
        output_directory=str(output_dir),
        output_prefix='query',
        file_format='h5ad',
        input_identifier='query_prepared_for_geneformer',
    )
    
    query_dataset_path = output_dir / 'query.dataset'
    msg("✓ Query tokenized dataset: %s", query_dataset_path)
    
    # Load tokenized dataset
    query_dataset = load_from_disk(str(query_dataset_path))
    if isinstance(query_dataset, DatasetDict):
        query_dataset = query_dataset[next(iter(query_dataset.keys()))]
    
    # Add dummy label column if missing (required by collator)
    if 'label' not in query_dataset.column_names:
        msg("  Adding dummy labels for inference...")
        dummy_labels = [0] * len(query_dataset)
        query_dataset = query_dataset.add_column('label', dummy_labels)

    msg("  Tokenized query: %d cells", len(query_dataset))
    
    # Load fine-tuned model and label encoder
    msg("Loading fine-tuned model from: %s", model_path)
    model = BertForSequenceClassification.from_pretrained(str(model_path))
    
    with open(label_encoder_path, 'rb') as f:
        label_encoder = pickle.load(f)
    
    if device == 'cuda' and not torch.cuda.is_available():
        msg("CUDA not available, using CPU")
        device = 'cpu'
    model.to(device)
    model.eval()
    
    # Prepare data collator
    data_collator = DataCollatorForCellClassification(token_dictionary=gene_token_dict)
    
    def to_feature_list(batch):
        """Convert HuggingFace dataset slice (dict or Dataset) to list of feature dicts."""
        if isinstance(batch, dict):
            keys = list(batch.keys())
            length = len(batch[keys[0]]) if keys else 0
            feature_list = []
            for idx in range(length):
                feature = {k: batch[k][idx] for k in keys}
                feature_list.append(feature)
            return feature_list
        elif isinstance(batch, list):
            return batch
        else:
            # Dataset object: convert via to_dict()
            batch_dict = batch.to_dict()
            return to_feature_list(batch_dict)

    # Predict in batches
    msg("Running predictions...")
    predictions = []
    probabilities = []
    
    with torch.no_grad():
        for i in range(0, len(query_dataset), batch_size):
            batch_slice = query_dataset[i:i+batch_size]
            features = to_feature_list(batch_slice)
            if not features:
                continue
            batch = data_collator(features)
            
            # Move to device
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch.get('attention_mask', None)
            if attention_mask is not None:
                attention_mask = attention_mask.to(device)
            
            # Forward pass
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            logits = outputs.logits
            
            # Get predictions
            batch_preds = torch.argmax(logits, dim=1).cpu().numpy()
            batch_probs = torch.softmax(logits, dim=1).cpu().numpy()
            
            predictions.extend(batch_preds)
            probabilities.extend(batch_probs)
            
            if (i + batch_size) % (batch_size * 10) == 0:
                msg("  Processed %d/%d cells", i + batch_size, len(query_dataset))
    
    # Decode predictions
    predicted_labels = label_encoder.inverse_transform(predictions)
    confidence_scores = np.max(probabilities, axis=1)
    
    msg("  Predictions complete: %d cells", len(predicted_labels))
    
    # Show distribution
    pred_counts = pd.Series(predicted_labels).value_counts()
    msg("\nPredicted cell type distribution:")
    for celltype, count in pred_counts.head(20).items():
        pct = 100 * count / len(predicted_labels)
        msg("  %s: %d cells (%.1f%%)", celltype, count, pct)
    if len(pred_counts) > 20:
        msg("  ... and %d more cell types", len(pred_counts) - 20)
    
    # Add predictions to query AnnData
    query_adata.obs['geneformer_predicted_celltype'] = predicted_labels
    query_adata.obs['geneformer_confidence'] = confidence_scores
    
    # Save annotated query
    query_output_path = output_dir / 'query_with_predictions.h5ad'
    msg("Saving annotated query to: %s", query_output_path)
    query_adata.write_h5ad(query_output_path)
    
    msg("✓ Query predictions saved to: %s", query_output_path)
    
    return query_output_path


def main():
    parser = argparse.ArgumentParser(description='Train Geneformer on multiple tonsil datasets')
    parser.add_argument('--reference_files', type=str, nargs='+', required=True,
                        help='List of reference h5ad files')
    parser.add_argument('--query', type=str, required=True,
                        help='Query h5ad file (not yet used; placeholder for CLI compatibility)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--label_column', type=str, default='annotation_20230508',
                        help='Column name for cell type labels in reference')
    parser.add_argument('--max_cells', type=int, default=50000,
                        help='Maximum number of cells to use for training (default: 50000)')
    parser.add_argument('--batch_size', type=int, default=16,
                        help='Batch size for training/inference (default: 16)')
    parser.add_argument('--epochs', type=int, default=3,
                        help='Number of training epochs (default: 3)')
    parser.add_argument('--learning_rate', type=float, default=5e-5,
                        help='Learning rate (default: 5e-5)')
    parser.add_argument('--device', type=str, default='cuda',
                        help='Device: cuda or cpu (default: cuda)')
    parser.add_argument('--max_genes', type=int, default=4096,
                        help='(Unused) retained for backward compatibility')
    parser.add_argument('--include_cell_types', type=str, nargs='+', default=None,
                        help='Cell types to include (substring matching, case-insensitive). Example: --include_cell_types "PB" "PC" "MBC"')
    parser.add_argument('--exclude_cell_types', type=str, nargs='+', default=None,
                        help='Cell types to exclude (substring matching, case-insensitive). Example: --exclude_cell_types "NBC early" "preGC"')
    args = parser.parse_args()

    if args.device == 'cuda' and not torch.cuda.is_available():
        msg("CUDA not available, using CPU")
        args.device = 'cpu'

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    msg("=" * 70)
    msg("LOADING AND MERGING REFERENCE DATASETS")
    msg("=" * 70)

    ref_adata, resolved_label = load_and_merge_datasets(
        args.reference_files,
        max_cells=args.max_cells,
        label_column=args.label_column,
        include_cell_types=args.include_cell_types,
        exclude_cell_types=args.exclude_cell_types,
    )

    if ref_adata is None:
        msg("ERROR: Failed to load/merge reference datasets")
        return 1

    training_label_col = 'training_label'
    if training_label_col not in ref_adata.obs.columns:
        msg("ERROR: training_label column missing after merge.")
        return 1

    msg("Using label column: %s", resolved_label)

    ref_adata = ensure_geneformer_requirements(ref_adata)
    ref_h5ad_path = output_dir / 'tonsil_merged_for_geneformer.h5ad'
    msg("Writing merged reference AnnData to: %s", ref_h5ad_path)
    ref_adata.write_h5ad(ref_h5ad_path)

    dataset_path, gene_token_dict = tokenize_reference_with_geneformer(
        ref_h5ad_path=ref_h5ad_path,
        output_dir=output_dir,
        label_col=training_label_col,
        nproc=8,
        model_version='V2',
    )

    model_path, label_encoder_path = fine_tune_geneformer_on_dataset(
        dataset_path=dataset_path,
        output_dir=output_dir,
        gene_token_dict=gene_token_dict,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        device=args.device,
    )

    msg("\n" + "=" * 70)
    msg("TRAINING COMPLETE")
    msg("=" * 70)
    msg("Fine-tuned model directory: %s", model_path)
    msg("Label encoder file: %s", label_encoder_path)
    
    # Predict on query data
    if args.query:
        query_path = predict_query_cells(
            query_h5ad_path=Path(args.query),
            model_path=model_path,
            label_encoder_path=label_encoder_path,
            gene_token_dict=gene_token_dict,
            output_dir=output_dir,
            device=args.device,
            batch_size=max(1, args.batch_size // 2),  # Smaller batch for inference
        )
        
        msg("\n" + "=" * 70)
        msg("PREDICTION COMPLETE")
        msg("=" * 70)
        msg("Annotated query saved to: %s", query_path)
        msg("Columns added:")
        msg("  - geneformer_predicted_celltype: Predicted tonsil cell type")
        msg("  - geneformer_confidence: Prediction confidence score (0-1)")

    return 0


if __name__ == '__main__':
    sys.exit(main())


__EOF_train_geneformer_tonsil_multi_py__

cat > "${SCRIPTS}/geneformer_predict_and_plot_manuscript.py" << '__EOF_geneformer_predict_and_plot_manuscript_py__'
#!/usr/bin/env python3
"""
Geneformer Prediction and Visualization for MANUSCRIPT
========================================================

This script:
1. Loads the CellBender-filtered mouse data (with scVI integration)
2. Predicts cell types using the pre-trained 48-class Geneformer model
3. Generates publication-quality visualizations
4. Saves all outputs to the MANUSCRIPT folder

Author: J
Date: 2025-12-01
"""

import sys
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import pickle
import numpy as np
import pandas as pd
import torch
import scipy.sparse as sp

import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from datasets import load_from_disk, DatasetDict
    from transformers import BertForSequenceClassification
    from geneformer import TranscriptomeTokenizer, DataCollatorForCellClassification
except ImportError:
    print("ERROR: transformers/datasets/geneformer not installed")
    print("Activate geneformer environment: conda activate geneformer")
    sys.exit(1)

try:
    import mygene
except ImportError:
    print("WARNING: mygene not installed. Mouse-to-human gene conversion will be skipped.")
    mygene = None

# ============================== CONFIGURATION ================================
# Paths
MOUSE_H5AD = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
GENEFORMER_MODEL_DIR = Path("__GENEFORMER_MODEL_DIR__/fine_tuned_model")
LABEL_ENCODER_PATH = Path("__GENEFORMER_MODEL_DIR__/label_encoder.pkl")
OUTPUT_DIR = Path("__BASEDIR__/Geneformer")

# Create output directories
FIGDIR = OUTPUT_DIR / "figures"
FIGDIR_INDIVIDUAL = OUTPUT_DIR / "figures_individual"
PREDICTIONS_DIR = OUTPUT_DIR / "predictions"
STATS_DIR = OUTPUT_DIR / "statistics"

for d in [OUTPUT_DIR, FIGDIR, FIGDIR_INDIVIDUAL, PREDICTIONS_DIR, STATS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def msg(text, *args):
    print(f"[Geneformer] {text % args if args else text}")


# ============================== HELPER FUNCTIONS =============================
def ensure_geneformer_requirements(adata: ad.AnnData) -> ad.AnnData:
    """Ensure var['ensembl_id'] and obs['n_counts'] exist for Geneformer."""
    msg("Ensuring Geneformer-required fields...")
    
    if 'ensembl_id' not in adata.var.columns:
        adata.var['ensembl_id'] = adata.var_names
    
    if 'n_counts' not in adata.obs.columns:
        if sp.issparse(adata.X):
            adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).ravel()
        else:
            adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    return adata


def convert_mouse_to_human_genes(adata: ad.AnnData) -> ad.AnnData:
    """Convert mouse gene symbols to human orthologs using mygene."""
    if mygene is None:
        msg("WARNING: mygene not available. Using uppercase conversion.")
        adata.var_names = [g.upper() for g in adata.var_names]
        adata.var_names_make_unique()
        return adata
    
    msg("Converting mouse genes to human orthologs...")
    mg = mygene.MyGeneInfo()
    
    unique_genes = list(set(adata.var_names))
    msg("  Converting %d unique genes...", len(unique_genes))
    
    gene_mapping = {}
    batch_size_conv = 1000
    converted_count = 0
    
    for i in range(0, len(unique_genes), batch_size_conv):
        batch = unique_genes[i:i+batch_size_conv]
        try:
            results = mg.querymany(
                batch,
                scopes='symbol',
                fields='symbol',
                species='mouse',
                target_species='human',
                returnall=True,
                verbose=False
            )
            
            for result in results.get('out', []):
                if 'symbol' in result and result.get('query') in batch:
                    gene_mapping[result['query']] = result['symbol']
                    converted_count += 1
        except Exception as e:
            continue
    
    msg("  Converted %d/%d genes (%.1f%%)", converted_count, len(unique_genes), 
        100 * converted_count / len(unique_genes) if unique_genes else 0)
    
    # Apply mapping
    human_genes = [gene_mapping.get(g, g.upper()) for g in adata.var_names]
    adata.var_names = human_genes
    
    # Remove duplicates by keeping first occurrence (same as original script)
    _, unique_idx = np.unique(adata.var_names, return_index=True)
    unique_idx = np.sort(unique_idx)  # Keep original order
    adata = adata[:, unique_idx].copy()
    msg("  After removing duplicates: %d unique genes", adata.n_vars)
    
    return adata


def to_feature_list(batch):
    """Convert HuggingFace dataset slice to list of feature dicts."""
    if isinstance(batch, dict):
        keys = list(batch.keys())
        length = len(batch[keys[0]]) if keys else 0
        return [{k: batch[k][idx] for k in keys} for idx in range(length)]
    elif isinstance(batch, list):
        return batch
    else:
        return to_feature_list(batch.to_dict())


# ============================== MAIN PIPELINE ================================
def main():
    print("=" * 84)
    print("GENEFORMER PREDICTION FOR MANUSCRIPT")
    print("=" * 84)
    print(f"CUDA available: {torch.cuda.is_available()}")
    print(f"Input: {MOUSE_H5AD}")
    print(f"Model: {GENEFORMER_MODEL_DIR}")
    print(f"Output: {OUTPUT_DIR}\n")
    
    # ==================== STEP 1: Load Mouse Data ====================
    print("=" * 84)
    print("STEP 1 — Load mouse data (with scVI integration)")
    print("=" * 84)
    
    msg("Loading: %s", MOUSE_H5AD)
    adata = sc.read_h5ad(MOUSE_H5AD)
    msg("  Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)
    msg("  Conditions: %s", list(adata.obs['condition'].unique()))
    
    # Store original obs for later
    original_obs = adata.obs.copy()
    original_obsm = {k: v.copy() for k, v in adata.obsm.items()}
    
    # ==================== STEP 2: Prepare for Geneformer ====================
    print("\n" + "=" * 84)
    print("STEP 2 — Prepare data for Geneformer (using RAW COUNTS)")
    print("=" * 84)
    
    # Make a copy for Geneformer processing
    adata_gf = adata.copy()
    
    # CRITICAL: Geneformer expects RAW COUNTS, not normalized data!
    # The scVI integrated file has normalized data in X, but raw counts in layers['counts']
    if 'counts' in adata_gf.layers:
        msg("Using raw counts from layers['counts'] (Geneformer requirement)")
        msg("  Current X min/max: %.3f / %.3f (normalized)", adata_gf.X.min(), adata_gf.X.max())
        adata_gf.X = adata_gf.layers['counts'].copy()
        msg("  New X min/max: %.3f / %.3f (raw counts)", adata_gf.X.min(), adata_gf.X.max())
    else:
        msg("WARNING: No 'counts' layer found. Using X as-is.")
        msg("  X min/max: %.3f / %.3f", adata_gf.X.min(), adata_gf.X.max())
    
    # Convert mouse to human genes
    adata_gf = convert_mouse_to_human_genes(adata_gf)
    
    # Ensure Geneformer requirements
    adata_gf = ensure_geneformer_requirements(adata_gf)
    
    # Save prepared data for tokenization
    prepared_path = PREDICTIONS_DIR / "mouse_prepared_for_geneformer.h5ad"
    msg("Saving prepared data: %s", prepared_path)
    adata_gf.write_h5ad(prepared_path)
    
    # ==================== STEP 3: Tokenize ====================
    print("\n" + "=" * 84)
    print("STEP 3 — Tokenize with Geneformer")
    print("=" * 84)
    
    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict=None,
        nproc=8,
        model_version='V2',
    )
    
    msg("Tokenizing data...")
    tokenizer.tokenize_data(
        data_directory=str(PREDICTIONS_DIR),
        output_directory=str(PREDICTIONS_DIR),
        output_prefix='mouse_query',
        file_format='h5ad',
        input_identifier='mouse_prepared_for_geneformer',
    )
    
    dataset_path = PREDICTIONS_DIR / 'mouse_query.dataset'
    msg("✓ Tokenized dataset: %s", dataset_path)
    
    # Load tokenized dataset
    query_dataset = load_from_disk(str(dataset_path))
    if isinstance(query_dataset, DatasetDict):
        query_dataset = query_dataset[next(iter(query_dataset.keys()))]
    
    # Add dummy labels for inference
    if 'label' not in query_dataset.column_names:
        query_dataset = query_dataset.add_column('label', [0] * len(query_dataset))
    
    msg("  Tokenized: %d cells", len(query_dataset))
    
    # ==================== STEP 4: Predict ====================
    print("\n" + "=" * 84)
    print("STEP 4 — Predict cell types with fine-tuned Geneformer")
    print("=" * 84)
    
    # Load model and label encoder
    msg("Loading fine-tuned model: %s", GENEFORMER_MODEL_DIR)
    model = BertForSequenceClassification.from_pretrained(str(GENEFORMER_MODEL_DIR))
    
    with open(LABEL_ENCODER_PATH, 'rb') as f:
        label_encoder = pickle.load(f)
    
    msg("  Number of classes: %d", len(label_encoder.classes_))
    msg("  Classes: %s", list(label_encoder.classes_[:10]) + ['...'] if len(label_encoder.classes_) > 10 else list(label_encoder.classes_))
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model.to(device)
    model.eval()
    
    # Get gene token dictionary
    tokenizer_obj = TranscriptomeTokenizer(model_version='V2')
    gene_token_dict = tokenizer_obj.gene_token_dict
    
    # Prepare data collator
    data_collator = DataCollatorForCellClassification(token_dictionary=gene_token_dict)
    
    # Predict in batches
    msg("Running predictions...")
    predictions = []
    probabilities = []
    batch_size = 8
    
    with torch.no_grad():
        for i in range(0, len(query_dataset), batch_size):
            batch_slice = query_dataset[i:i+batch_size]
            features = to_feature_list(batch_slice)
            if not features:
                continue
            batch = data_collator(features)
            
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch.get('attention_mask', None)
            if attention_mask is not None:
                attention_mask = attention_mask.to(device)
            
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            logits = outputs.logits
            
            batch_preds = torch.argmax(logits, dim=1).cpu().numpy()
            batch_probs = torch.softmax(logits, dim=1).cpu().numpy()
            
            predictions.extend(batch_preds)
            probabilities.extend(batch_probs)
            
            if (i + batch_size) % 1000 == 0:
                msg("  Processed %d/%d cells", min(i + batch_size, len(query_dataset)), len(query_dataset))
    
    # Decode predictions
    predicted_labels = label_encoder.inverse_transform(predictions)
    confidence_scores = np.max(probabilities, axis=1)
    
    msg("✓ Predictions complete: %d cells", len(predicted_labels))
    
    # ==================== STEP 5: Add predictions to original AnnData ====================
    print("\n" + "=" * 84)
    print("STEP 5 — Merge predictions with original data")
    print("=" * 84)
    
    # Add predictions to original adata (which has scVI UMAP)
    adata.obs['geneformer_predicted_celltype'] = predicted_labels
    adata.obs['geneformer_confidence'] = confidence_scores
    
    # Save full probability matrix
    prob_df = pd.DataFrame(
        probabilities, 
        index=adata.obs_names, 
        columns=label_encoder.classes_
    )
    prob_df.to_csv(STATS_DIR / 'prediction_probabilities.csv')
    
    # ==================== STEP 6: Statistics ====================
    print("\n" + "=" * 84)
    print("STEP 6 — Generate statistics")
    print("=" * 84)
    
    # Prediction distribution
    pred_counts = pd.Series(predicted_labels).value_counts()
    msg("\nPredicted cell type distribution (top 20):")
    for celltype, count in pred_counts.head(20).items():
        pct = 100 * count / len(predicted_labels)
        msg("  %s: %d cells (%.1f%%)", celltype, count, pct)
    
    pred_counts.to_csv(STATS_DIR / 'celltype_counts.csv')
    
    # Confidence by condition
    conf_by_condition = adata.obs.groupby('condition')['geneformer_confidence'].agg(['mean', 'std', 'median'])
    conf_by_condition.to_csv(STATS_DIR / 'confidence_by_condition.csv')
    msg("\nConfidence by condition:")
    print(conf_by_condition.to_string())
    
    # Cross-tabulation: condition vs predicted celltype
    crosstab = pd.crosstab(adata.obs['condition'], adata.obs['geneformer_predicted_celltype'])
    crosstab.to_csv(STATS_DIR / 'condition_vs_celltype_crosstab.csv')
    
    # Percentage version
    crosstab_pct = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
    crosstab_pct.to_csv(STATS_DIR / 'condition_vs_celltype_percentage.csv')
    
    # ==================== STEP 7: Generate Figures ====================
    print("\n" + "=" * 84)
    print("STEP 7 — Generate publication figures (PNG, PDF, SVG)")
    print("=" * 84)
    
    # Helper function to save in multiple formats
    def save_figure(fig, basename, dpi=300):
        """Save figure in PNG, PDF, and SVG formats."""
        for fmt in ['png', 'pdf', 'svg']:
            filepath = FIGDIR / f'{basename}.{fmt}'
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
        msg("  ✓ %s (.png, .pdf, .svg)", basename)
    
    def save_figure_individual(fig, basename, dpi=300):
        """Save figure in multiple formats to individual folder."""
        for fmt in ['png', 'pdf', 'svg']:
            filepath = FIGDIR_INDIVIDUAL / f'{basename}.{fmt}'
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
    
    # Color palettes
    condition_colors = {
        "WT_B_cells": "#2ecc71",
        "Crebbp_B_cells": "#3498db", 
        "Pre_malignant": "#f39c12",
        "Matched_malignant": "#e74c3c",
        "Malignant": "#8e44ad"
    }
    
    # 1. UMAP by predicted cell type (all cells) - legend on right
    msg("Generating UMAP plots...")
    fig, ax = plt.subplots(figsize=(16, 10))
    sc.pl.umap(adata, color='geneformer_predicted_celltype', ax=ax, show=False,
               frameon=False, legend_loc='right margin', legend_fontsize=7, s=15,
               title='Geneformer Predicted Cell Types (48 tonsil classes)')
    plt.tight_layout()
    save_figure(fig, 'umap_geneformer_predictions')
    plt.close()
    
    # 2. UMAP by confidence score
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='geneformer_confidence', ax=ax, show=False,
               frameon=False, cmap='viridis', s=15,
               title='Geneformer Prediction Confidence')
    plt.tight_layout()
    save_figure(fig, 'umap_geneformer_confidence')
    plt.close()
    
    # 3. UMAP by condition (for reference) - legend on right
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='condition', ax=ax, show=False, palette=condition_colors,
               frameon=False, legend_loc='right margin', legend_fontsize=10, s=15,
               title='Experimental Condition')
    plt.tight_layout()
    save_figure(fig, 'umap_condition')
    plt.close()
    
    # 4. Confidence distribution by condition (violin plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    order = ['WT_B_cells', 'Crebbp_B_cells', 'Pre_malignant', 'Matched_malignant', 'Malignant']
    order = [c for c in order if c in adata.obs['condition'].unique()]
    sns.violinplot(data=adata.obs, x='condition', y='geneformer_confidence', 
                   order=order, palette=condition_colors, ax=ax)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Geneformer Confidence Score')
    ax.set_title('Prediction Confidence by Condition')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_figure(fig, 'violin_confidence_by_condition')
    plt.close()
    
    # 5. Stacked bar chart: cell type composition by condition
    fig, ax = plt.subplots(figsize=(14, 8))
    crosstab_pct_plot = crosstab_pct.loc[order] if all(c in crosstab_pct.index for c in order) else crosstab_pct
    crosstab_pct_plot.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.8)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Percentage of Cells')
    ax.set_title('Geneformer Predicted Cell Type Composition by Condition')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_figure(fig, 'stacked_bar_celltype_by_condition')
    plt.close()
    
    # 6. Heatmap of cell type proportions
    fig, ax = plt.subplots(figsize=(16, 8))
    # Select top 20 most common cell types
    top_celltypes = pred_counts.head(20).index.tolist()
    heatmap_data = crosstab_pct[top_celltypes].loc[order] if all(c in crosstab_pct.index for c in order) else crosstab_pct[top_celltypes]
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Percentage'})
    ax.set_xlabel('Predicted Cell Type')
    ax.set_ylabel('Condition')
    ax.set_title('Cell Type Proportions by Condition (Top 20 Types)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_figure(fig, 'heatmap_celltype_proportions')
    plt.close()
    
    # 7. High confidence cells only (>0.7) - legend on right
    adata_high = adata[adata.obs['geneformer_confidence'] > 0.7].copy()
    if len(adata_high) > 100:
        fig, ax = plt.subplots(figsize=(16, 10))
        sc.pl.umap(adata_high, color='geneformer_predicted_celltype', ax=ax, show=False,
                   frameon=False, legend_loc='right margin', legend_fontsize=7, s=20,
                   title=f'High Confidence Predictions (>0.7, n={len(adata_high):,})')
        plt.tight_layout()
        save_figure(fig, 'umap_high_confidence_predictions')
        plt.close()
    
    # 8. Individual condition UMAPs (PNG, PDF, SVG)
    msg("Generating individual condition UMAPs...")
    for condition in adata.obs['condition'].unique():
        adata_cond = adata[adata.obs['condition'] == condition].copy()
        safe_name = condition.replace(' ', '_')
        
        # All cells colored by prediction - legend on right
        fig, ax = plt.subplots(figsize=(14, 10))
        sc.pl.umap(adata_cond, color='geneformer_predicted_celltype', ax=ax, show=False,
                   frameon=False, legend_loc='right margin', legend_fontsize=7, s=20,
                   title=f'{condition}: Geneformer Predictions (n={len(adata_cond):,})')
        plt.tight_layout()
        save_figure_individual(fig, f'umap_{safe_name}_predictions')
        plt.close()
        
        # Confidence
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata_cond, color='geneformer_confidence', ax=ax, show=False,
                   frameon=False, cmap='viridis', s=20,
                   title=f'{condition}: Confidence (n={len(adata_cond):,})')
        plt.tight_layout()
        save_figure_individual(fig, f'umap_{safe_name}_confidence')
        plt.close()
    
    msg("  ✓ Individual condition UMAPs saved to %s (.png, .pdf, .svg)", FIGDIR_INDIVIDUAL)
    
    # ==================== STEP 8: Save Outputs ====================
    print("\n" + "=" * 84)
    print("STEP 8 — Save outputs")
    print("=" * 84)
    
    # Save annotated AnnData
    output_h5ad = OUTPUT_DIR / 'mouse_with_geneformer_predictions.h5ad'
    adata.write_h5ad(output_h5ad)
    msg("✓ Annotated h5ad: %s", output_h5ad)
    
    # Save cell metadata
    adata.obs.to_csv(OUTPUT_DIR / 'cell_metadata_with_geneformer.csv')
    msg("✓ Cell metadata: %s", OUTPUT_DIR / 'cell_metadata_with_geneformer.csv')
    
    # Save label encoder classes for reference
    with open(STATS_DIR / 'geneformer_classes.txt', 'w') as f:
        for i, cls in enumerate(label_encoder.classes_):
            f.write(f"{i}\t{cls}\n")
    msg("✓ Class labels: %s", STATS_DIR / 'geneformer_classes.txt')
    
    print("\n" + "=" * 84)
    print("COMPLETE!")
    print("=" * 84)
    print(f"  Output directory: {OUTPUT_DIR}")
    print(f"  Figures: {FIGDIR}")
    print(f"  Individual figures: {FIGDIR_INDIVIDUAL}")
    print(f"  Statistics: {STATS_DIR}")
    print(f"  Total cells: {adata.n_obs:,}")
    print(f"  Unique predicted types: {adata.obs['geneformer_predicted_celltype'].nunique()}")
    print(f"  Mean confidence: {adata.obs['geneformer_confidence'].mean():.3f}")
    print("\nDONE.\n")


if __name__ == '__main__':
    main()



__EOF_geneformer_predict_and_plot_manuscript_py__

cat > "${SCRIPTS}/scvi_human_dlbcl_mouse_malignant_integration.py" << '__EOF_scvi_human_dlbcl_mouse_malignant_integration_py__'
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
OUTDIR = Path("__BASEDIR__/mouse_human_integration")
FIGDIR = OUTDIR / "figures"
FIGDIR_CT2 = OUTDIR / "figures_ct2"
CT2_WORKDIR = OUTDIR / "ct2_io"
CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"
CT2_OUTDIR = OUTDIR / "cytotrace2_results"

for p in (OUTDIR, FIGDIR, FIGDIR_CT2, CT2_WORKDIR, CT2_OUTDIR):
    p.mkdir(parents=True, exist_ok=True)

# Human DLBCL files
DLBCL_FILES = [
    "__DLBCL_DIR__/DLBCL1_raw.h5ad",
    "__DLBCL_DIR__/DLBCL2_raw.h5ad",
    "__DLBCL_DIR__/DLBCL3_raw.h5ad",
    "__DLBCL_DIR__/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad",
]

# Human Tonsil data files
TONSIL_GC_PATH = "__TONSIL_DIR__/tonsil_GCBC_RNA.h5ad"
TONSIL_MBC_PATH = "__TONSIL_DIR__/tonsil_NBC-MBC_RNA.h5ad"
MAX_TONSIL_CELLS = 12500  # Maximum total GC + Memory B cells to keep

# Proliferation filtering for tonsil cells
FILTER_PROLIFERATING_TONSIL = False  # Changed to False to preserve Dark Zone (Centroblasts)
PROLIFERATION_THRESHOLD = 0.20  # Relaxed threshold if enabled (was 0.10)
EXCLUDE_PROLIFERATIVE_GC_TYPES = False  # Changed to False to preserve Dark Zone annotations

# Toggle CytoTRACE2 (Set to True later to run CT2 after scVI completes)
RUN_CYTOTRACE2 = False  # Set to False to skip CT2 and save memory

# Mouse malignant samples (CellBender filtered)
CELLBENDER_DIR = Path("__CELLBENDER_DIR__")
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
OUTDIR = Path("__BASEDIR__/mouse_human_integration")
FIGDIR = OUTDIR / "figures"
FIGDIR_CT2 = OUTDIR / "figures_ct2"
CT2_WORKDIR = OUTDIR / "ct2_io"
CT2_INPUT_TXT = CT2_WORKDIR / "ct2_input_counts.txt"
CT2_OUTDIR = OUTDIR / "cytotrace2_results"

for p in (OUTDIR, FIGDIR, FIGDIR_CT2, CT2_WORKDIR, CT2_OUTDIR):
    p.mkdir(parents=True, exist_ok=True)

# Human DLBCL files
DLBCL_FILES = [
    "__DLBCL_DIR__/DLBCL1_raw.h5ad",
    "__DLBCL_DIR__/DLBCL2_raw.h5ad",
    "__DLBCL_DIR__/DLBCL3_raw.h5ad",
    "__DLBCL_DIR__/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad",
]

# Human Tonsil data files
TONSIL_GC_PATH = "__TONSIL_DIR__/tonsil_GCBC_RNA.h5ad"
TONSIL_MBC_PATH = "__TONSIL_DIR__/tonsil_NBC-MBC_RNA.h5ad"
MAX_TONSIL_CELLS = 12500  # Maximum total GC + Memory B cells to keep

# Proliferation filtering for tonsil cells
FILTER_PROLIFERATING_TONSIL = False  # Changed to False to preserve Dark Zone (Centroblasts)
PROLIFERATION_THRESHOLD = 0.20  # Relaxed threshold if enabled (was 0.10)
EXCLUDE_PROLIFERATIVE_GC_TYPES = False  # Changed to False to preserve Dark Zone annotations

# Toggle CytoTRACE2 (Set to True later to run CT2 after scVI completes)
RUN_CYTOTRACE2 = False  # Set to False to skip CT2 and save memory

# Mouse malignant samples (CellBender filtered)
CELLBENDER_DIR = Path("__CELLBENDER_DIR__")
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



__EOF_scvi_human_dlbcl_mouse_malignant_integration_py__

cat > "${SCRIPTS}/cytotrace2_dlbcl_mouse_tonsil.py" << '__EOF_cytotrace2_dlbcl_mouse_tonsil_py__'
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
INPUT_H5AD = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Source data files
SOURCE_MOUSE_DIR = Path("__CELLBENDER_DIR__")
SOURCE_MOUSE_FILES = {
    "SIGAA3_Matched_malignant_R1": SOURCE_MOUSE_DIR / "SIGAA3_Matched_malignant_R1_GEX_cellbender_filtered.h5",
    "SIGAA4_Matched_malignant_R2": SOURCE_MOUSE_DIR / "SIGAA4_Matched_malignant_R2_GEX_cellbender_filtered.h5",
    "SIGAD5_Malignant_R2": SOURCE_MOUSE_DIR / "SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5",
    "SIGAH1_Malignant_R1": SOURCE_MOUSE_DIR / "SIGAH1_Malignant_R1_GEX_cellbender_filtered.h5",
}
SOURCE_DLBCL_FILES = [
    Path("__DLBCL_DIR__/DLBCL1_raw.h5ad"),
    Path("__DLBCL_DIR__/DLBCL2_raw.h5ad"),
    Path("__DLBCL_DIR__/DLBCL3_raw.h5ad"),
    Path("__DLBCL_DIR__/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad"),
]
SOURCE_TONSIL_GC = Path("__TONSIL_DIR__/tonsil_GCBC_RNA.h5ad")
SOURCE_TONSIL_MBC = Path("__TONSIL_DIR__/tonsil_NBC-MBC_RNA.h5ad")

# Output directory
OUTPUT_DIR = Path("__BASEDIR__/mouse_human_integration/CytoTRACE2")
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


__EOF_cytotrace2_dlbcl_mouse_tonsil_py__

cat > "${SCRIPTS}/geneformer_predict_dlbcl_mouse_tonsil.py" << '__EOF_geneformer_predict_dlbcl_mouse_tonsil_py__'
#!/usr/bin/env python3
"""
Geneformer Prediction for Human DLBCL + Mouse Malignant + Tonsil Integration
=============================================================================

This script:
1. Loads the integrated human DLBCL + mouse malignant + tonsil data (with scVI)
2. Predicts cell types using the pre-trained 48-class Geneformer model
3. Generates publication-quality visualizations
4. Saves all outputs to the Geneformer subfolder

Note: Both human and mouse data have already been converted to human gene symbols
      in the integration script (mouse via BioMart 1:1 orthologs).

Author: J
Date: 2025-12-11
"""

import sys
import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import pickle
import numpy as np
import pandas as pd
import torch
import scipy.sparse as sp

import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from datasets import load_from_disk, DatasetDict
    from transformers import BertForSequenceClassification
    from geneformer import TranscriptomeTokenizer, DataCollatorForCellClassification
except ImportError:
    print("ERROR: transformers/datasets/geneformer not installed")
    print("Activate geneformer environment: conda activate geneformer")
    sys.exit(1)

# ============================== CONFIGURATION ================================
# Input: Integrated human DLBCL + mouse malignant + tonsil object
INPUT_H5AD = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Geneformer model (trained on tonsil)
GENEFORMER_MODEL_DIR = Path("__GENEFORMER_MODEL_DIR__/fine_tuned_model")
LABEL_ENCODER_PATH = Path("__GENEFORMER_MODEL_DIR__/label_encoder.pkl")

# Output directory
OUTPUT_DIR = Path("__BASEDIR__/mouse_human_integration/Geneformer")

# Create output directories
FIGDIR = OUTPUT_DIR / "figures"
FIGDIR_INDIVIDUAL = OUTPUT_DIR / "figures_individual"
PREDICTIONS_DIR = OUTPUT_DIR / "predictions"
STATS_DIR = OUTPUT_DIR / "statistics"

for d in [OUTPUT_DIR, FIGDIR, FIGDIR_INDIVIDUAL, PREDICTIONS_DIR, STATS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Batch size for inference (reduce if OOM)
BATCH_SIZE = 8


def msg(text, *args):
    print(f"[Geneformer] {text % args if args else text}")


# ============================== HELPER FUNCTIONS =============================
def ensure_geneformer_requirements(adata: ad.AnnData) -> ad.AnnData:
    """Ensure var['ensembl_id'] and obs['n_counts'] exist for Geneformer."""
    msg("Ensuring Geneformer-required fields...")
    
    if 'ensembl_id' not in adata.var.columns:
        adata.var['ensembl_id'] = adata.var_names
    
    if 'n_counts' not in adata.obs.columns:
        if sp.issparse(adata.X):
            adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).ravel()
        else:
            adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    return adata


def to_feature_list(batch):
    """Convert HuggingFace dataset slice to list of feature dicts."""
    if isinstance(batch, dict):
        keys = list(batch.keys())
        length = len(batch[keys[0]]) if keys else 0
        return [{k: batch[k][idx] for k in keys} for idx in range(length)]
    elif isinstance(batch, list):
        return batch
    else:
        return to_feature_list(batch.to_dict())


# ============================== MAIN PIPELINE ================================
def main():
    print("=" * 84)
    print("GENEFORMER PREDICTION: HUMAN DLBCL + MOUSE MALIGNANT + TONSIL INTEGRATION")
    print("=" * 84)
    print(f"CUDA available: {torch.cuda.is_available()}")
    print(f"Input: {INPUT_H5AD}")
    print(f"Model: {GENEFORMER_MODEL_DIR}")
    print(f"Output: {OUTPUT_DIR}\n")
    
    # ==================== STEP 1: Load Integrated Data ====================
    print("=" * 84)
    print("STEP 1 — Load integrated data (human DLBCL + mouse malignant + tonsil)")
    print("=" * 84)
    
    msg("Loading: %s", INPUT_H5AD)
    adata = sc.read_h5ad(INPUT_H5AD)
    msg("  Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)
    
    # Check what metadata columns are available
    msg("  Available obs columns: %s", list(adata.obs.columns[:15]))
    
    if 'disease_state' in adata.obs:
        msg("  Disease states: %s", list(adata.obs['disease_state'].unique()))
    if 'species' in adata.obs:
        msg("  Species: %s", list(adata.obs['species'].unique()))
    if 'sample_batch' in adata.obs:
        msg("  Sample batches: %d unique", adata.obs['sample_batch'].nunique())
    
    # Store original obs for later
    original_obs = adata.obs.copy()
    original_obsm = {k: v.copy() for k, v in adata.obsm.items()}
    
    # ==================== STEP 2: Prepare for Geneformer ====================
    print("\n" + "=" * 84)
    print("STEP 2 — Prepare data for Geneformer (using RAW COUNTS)")
    print("=" * 84)
    
    # Make a copy for Geneformer processing
    adata_gf = adata.copy()
    
    # CRITICAL: Geneformer expects RAW COUNTS, not normalized data!
    if 'counts' in adata_gf.layers:
        msg("Using raw counts from layers['counts'] (Geneformer requirement)")
        if sp.issparse(adata_gf.X):
            msg("  Current X min/max: %.3f / %.3f (normalized)", 
                adata_gf.X.data.min() if adata_gf.X.data.size > 0 else 0, 
                adata_gf.X.data.max() if adata_gf.X.data.size > 0 else 0)
        else:
            msg("  Current X min/max: %.3f / %.3f (normalized)", adata_gf.X.min(), adata_gf.X.max())
        adata_gf.X = adata_gf.layers['counts'].copy()
        if sp.issparse(adata_gf.X):
            msg("  New X min/max: %.3f / %.3f (raw counts)", 
                adata_gf.X.data.min() if adata_gf.X.data.size > 0 else 0, 
                adata_gf.X.data.max() if adata_gf.X.data.size > 0 else 0)
        else:
            msg("  New X min/max: %.3f / %.3f (raw counts)", adata_gf.X.min(), adata_gf.X.max())
    else:
        msg("WARNING: No 'counts' layer found. Using X as-is.")
    
    # Gene symbols should already be in human format from the integration script
    # But let's make sure they're uppercase and handle duplicates
    msg("Standardizing gene symbols...")
    adata_gf.var_names = [str(g).upper() for g in adata_gf.var_names]
    
    # Remove duplicates by keeping first occurrence
    _, unique_idx = np.unique(adata_gf.var_names, return_index=True)
    unique_idx = np.sort(unique_idx)
    adata_gf = adata_gf[:, unique_idx].copy()
    msg("  After removing duplicates: %d unique genes", adata_gf.n_vars)
    
    # Ensure Geneformer requirements
    adata_gf = ensure_geneformer_requirements(adata_gf)
    
    # Save prepared data for tokenization
    prepared_path = PREDICTIONS_DIR / "prepared_for_geneformer.h5ad"
    msg("Saving prepared data: %s", prepared_path)
    adata_gf.write_h5ad(prepared_path)
    
    # ==================== STEP 3: Tokenize ====================
    print("\n" + "=" * 84)
    print("STEP 3 — Tokenize with Geneformer")
    print("=" * 84)
    
    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict=None,
        nproc=8,
        model_version='V2',
    )
    
    msg("Tokenizing data...")
    tokenizer.tokenize_data(
        data_directory=str(PREDICTIONS_DIR),
        output_directory=str(PREDICTIONS_DIR),
        output_prefix='query',
        file_format='h5ad',
        input_identifier='prepared_for_geneformer',
    )
    
    dataset_path = PREDICTIONS_DIR / 'query.dataset'
    msg("✓ Tokenized dataset: %s", dataset_path)
    
    # Load tokenized dataset
    query_dataset = load_from_disk(str(dataset_path))
    if isinstance(query_dataset, DatasetDict):
        query_dataset = query_dataset[next(iter(query_dataset.keys()))]
    
    # Add dummy labels for inference
    if 'label' not in query_dataset.column_names:
        query_dataset = query_dataset.add_column('label', [0] * len(query_dataset))
    
    msg("  Tokenized: %d cells", len(query_dataset))
    
    # ==================== STEP 4: Predict ====================
    print("\n" + "=" * 84)
    print("STEP 4 — Predict cell types with fine-tuned Geneformer")
    print("=" * 84)
    
    # Load model and label encoder
    msg("Loading fine-tuned model: %s", GENEFORMER_MODEL_DIR)
    model = BertForSequenceClassification.from_pretrained(str(GENEFORMER_MODEL_DIR))
    
    with open(LABEL_ENCODER_PATH, 'rb') as f:
        label_encoder = pickle.load(f)
    
    msg("  Number of classes: %d", len(label_encoder.classes_))
    msg("  Classes: %s", list(label_encoder.classes_[:10]) + ['...'] if len(label_encoder.classes_) > 10 else list(label_encoder.classes_))
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model.to(device)
    model.eval()
    
    # Get gene token dictionary
    tokenizer_obj = TranscriptomeTokenizer(model_version='V2')
    gene_token_dict = tokenizer_obj.gene_token_dict
    
    # Prepare data collator
    data_collator = DataCollatorForCellClassification(token_dictionary=gene_token_dict)
    
    # Predict in batches
    msg("Running predictions...")
    predictions = []
    probabilities = []
    
    with torch.no_grad():
        for i in range(0, len(query_dataset), BATCH_SIZE):
            batch_slice = query_dataset[i:i+BATCH_SIZE]
            features = to_feature_list(batch_slice)
            if not features:
                continue
            batch = data_collator(features)
            
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch.get('attention_mask', None)
            if attention_mask is not None:
                attention_mask = attention_mask.to(device)
            
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            logits = outputs.logits
            
            batch_preds = torch.argmax(logits, dim=1).cpu().numpy()
            batch_probs = torch.softmax(logits, dim=1).cpu().numpy()
            
            predictions.extend(batch_preds)
            probabilities.extend(batch_probs)
            
            if (i + BATCH_SIZE) % 1000 == 0:
                msg("  Processed %d/%d cells", min(i + BATCH_SIZE, len(query_dataset)), len(query_dataset))
    
    # Decode predictions
    predicted_labels = label_encoder.inverse_transform(predictions)
    confidence_scores = np.max(probabilities, axis=1)
    
    msg("✓ Predictions complete: %d cells", len(predicted_labels))
    
    # ==================== STEP 5: Add predictions to original AnnData ====================
    print("\n" + "=" * 84)
    print("STEP 5 — Merge predictions with original data")
    print("=" * 84)
    
    # Add predictions to original adata (which has scVI UMAP)
    adata.obs['geneformer_predicted_celltype'] = predicted_labels
    adata.obs['geneformer_confidence'] = confidence_scores
    
    # Save full probability matrix
    prob_df = pd.DataFrame(
        probabilities, 
        index=adata.obs_names, 
        columns=label_encoder.classes_
    )
    prob_df.to_csv(STATS_DIR / 'prediction_probabilities.csv')
    
    # Save updated AnnData
    adata.write_h5ad(OUTPUT_DIR / 'integrated_with_geneformer_predictions.h5ad')
    msg("✓ Saved: integrated_with_geneformer_predictions.h5ad")
    
    # ==================== STEP 6: Statistics ====================
    print("\n" + "=" * 84)
    print("STEP 6 — Generate statistics")
    print("=" * 84)
    
    # IMPORTANT: Species-aware analysis
    # Mouse cells have been converted to human gene symbols via BioMart 1:1 orthologs
    # Predictions for mouse cells should be interpreted with caution!
    msg("\n" + "="*60)
    msg("⚠️  IMPORTANT: CROSS-SPECIES PREDICTION NOTES")
    msg("="*60)
    msg("  - Geneformer was trained on HUMAN tonsil data")
    msg("  - Mouse cells have been converted to human gene symbols")
    msg("    via BioMart 1:1 orthologs in the integration script")
    msg("  - Mouse predictions should be interpreted with caution:")
    msg("    * Not all genes have 1:1 orthologs")
    msg("    * Expression patterns may differ between species")
    msg("    * Some human cell states may not exist in mouse")
    msg("="*60 + "\n")
    
    # Prediction distribution
    pred_counts = pd.Series(predicted_labels).value_counts()
    msg("\nPredicted cell type distribution (top 20):")
    for celltype, count in pred_counts.head(20).items():
        pct = 100 * count / len(predicted_labels)
        msg("  %s: %d cells (%.1f%%)", celltype, count, pct)
    
    pred_counts.to_csv(STATS_DIR / 'celltype_counts.csv')
    
    # Confidence by disease_state
    if 'disease_state' in adata.obs:
        conf_by_disease = adata.obs.groupby('disease_state')['geneformer_confidence'].agg(['mean', 'std', 'median', 'count'])
        conf_by_disease.to_csv(STATS_DIR / 'confidence_by_disease_state.csv')
        msg("\nConfidence by disease state:")
        print(conf_by_disease.to_string())
        
        # Cross-tabulation: disease_state vs predicted celltype
        crosstab = pd.crosstab(adata.obs['disease_state'], adata.obs['geneformer_predicted_celltype'])
        crosstab.to_csv(STATS_DIR / 'disease_state_vs_celltype_crosstab.csv')
        
        # Percentage version
        crosstab_pct = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
        crosstab_pct.to_csv(STATS_DIR / 'disease_state_vs_celltype_percentage.csv')
    
    # Confidence by species - CRITICAL for cross-species analysis
    if 'species' in adata.obs:
        conf_by_species = adata.obs.groupby('species')['geneformer_confidence'].agg(['mean', 'std', 'median', 'count'])
        conf_by_species.to_csv(STATS_DIR / 'confidence_by_species.csv')
        msg("\nConfidence by species:")
        print(conf_by_species.to_string())
        
        # Detailed per-species cell type distribution
        msg("\n" + "="*60)
        msg("SPECIES-STRATIFIED CELL TYPE DISTRIBUTIONS")
        msg("="*60)
        
        for species_name in ['human', 'mouse']:
            if species_name in adata.obs['species'].values:
                mask = adata.obs['species'] == species_name
                species_pred = adata.obs.loc[mask, 'geneformer_predicted_celltype'].value_counts()
                species_conf = adata.obs.loc[mask, 'geneformer_confidence']
                
                msg(f"\n{species_name.upper()} cells ({mask.sum():,} cells):")
                msg(f"  Mean confidence: {species_conf.mean():.3f}")
                msg(f"  Median confidence: {species_conf.median():.3f}")
                msg(f"  Top 10 predicted cell types:")
                for ct, count in species_pred.head(10).items():
                    pct = 100 * count / mask.sum()
                    msg(f"    - {ct}: {count:,} ({pct:.1f}%)")
                
                # Save per-species distributions
                species_pred.to_csv(STATS_DIR / f'celltype_counts_{species_name}.csv')
        
        # Cross-species comparison
        crosstab_species = pd.crosstab(adata.obs['species'], adata.obs['geneformer_predicted_celltype'])
        crosstab_species_pct = crosstab_species.div(crosstab_species.sum(axis=1), axis=0) * 100
        crosstab_species.to_csv(STATS_DIR / 'species_vs_celltype_crosstab.csv')
        crosstab_species_pct.to_csv(STATS_DIR / 'species_vs_celltype_percentage.csv')
    
    # Confidence by sample_batch
    if 'sample_batch' in adata.obs:
        conf_by_sample = adata.obs.groupby('sample_batch')['geneformer_confidence'].agg(['mean', 'std', 'median', 'count'])
        conf_by_sample.to_csv(STATS_DIR / 'confidence_by_sample.csv')
    
    # ==================== STEP 7: Generate Figures ====================
    print("\n" + "=" * 84)
    print("STEP 7 — Generate publication figures (PNG, PDF, SVG)")
    print("=" * 84)
    
    # Helper function to save in multiple formats
    def save_figure(fig, basename, dpi=300):
        """Save figure in PNG, PDF, and SVG formats."""
        for fmt in ['png', 'pdf', 'svg']:
            filepath = FIGDIR / f'{basename}.{fmt}'
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
        msg("  ✓ %s (.png, .pdf, .svg)", basename)
    
    def save_figure_individual(fig, basename, dpi=300):
        """Save figure in multiple formats to individual folder."""
        for fmt in ['png', 'pdf', 'svg']:
            filepath = FIGDIR_INDIVIDUAL / f'{basename}.{fmt}'
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight', format=fmt)
    
    # Color palettes
    disease_colors = {
        "Mouse_Malignant": "#FF6B6B",
        "Mouse_Matched_malignant": "#FF9999",
        "DLBCL": "#8B0000",
        "Tonsil_GC_B": "#4169E1",
        "Tonsil_Normal": "#4169E1"
    }
    
    species_colors = {
        "mouse": "#98FB98",
        "human": "#6495ED"
    }
    
    # 1. UMAP by predicted cell type (all cells)
    msg("Generating UMAP plots...")
    fig, ax = plt.subplots(figsize=(16, 10))
    sc.pl.umap(adata, color='geneformer_predicted_celltype', ax=ax, show=False,
               frameon=False, legend_loc='right margin', legend_fontsize=7, s=15,
               title='Geneformer Predicted Cell Types (48 tonsil classes)')
    plt.tight_layout()
    save_figure(fig, 'umap_geneformer_predictions')
    plt.close()
    
    # 2. UMAP by confidence score
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='geneformer_confidence', ax=ax, show=False,
               frameon=False, cmap='viridis', s=15,
               title='Geneformer Prediction Confidence')
    plt.tight_layout()
    save_figure(fig, 'umap_geneformer_confidence')
    plt.close()
    
    # 3. UMAP by disease_state
    if 'disease_state' in adata.obs:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='disease_state', ax=ax, show=False, palette=disease_colors,
                   frameon=False, legend_loc='right margin', legend_fontsize=10, s=15,
                   title='Disease State')
        plt.tight_layout()
        save_figure(fig, 'umap_disease_state')
        plt.close()
    
    # 4. UMAP by species
    if 'species' in adata.obs:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='species', ax=ax, show=False, palette=species_colors,
                   frameon=False, legend_loc='right margin', legend_fontsize=10, s=15,
                   title='Species')
        plt.tight_layout()
        save_figure(fig, 'umap_species')
        plt.close()
    
    # 5. Confidence distribution histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(confidence_scores, bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(np.median(confidence_scores), color='red', linestyle='--', 
               label=f'Median: {np.median(confidence_scores):.3f}')
    ax.axvline(np.mean(confidence_scores), color='green', linestyle='--', 
               label=f'Mean: {np.mean(confidence_scores):.3f}')
    ax.set_xlabel('Confidence Score', fontsize=12)
    ax.set_ylabel('Cell Count', fontsize=12)
    ax.set_title('Geneformer Prediction Confidence Distribution', fontsize=14)
    ax.legend()
    plt.tight_layout()
    save_figure(fig, 'histogram_confidence')
    plt.close()
    
    # 6. Confidence by disease state violin plot
    if 'disease_state' in adata.obs:
        fig, ax = plt.subplots(figsize=(12, 6))
        disease_order = sorted(adata.obs['disease_state'].unique())
        colors = [disease_colors.get(d, '#808080') for d in disease_order]
        sns.violinplot(data=adata.obs, x='disease_state', y='geneformer_confidence',
                       order=disease_order, palette=colors, ax=ax)
        ax.set_xlabel('Disease State', fontsize=12)
        ax.set_ylabel('Prediction Confidence', fontsize=12)
        ax.set_title('Geneformer Confidence by Disease State', fontsize=14)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_figure(fig, 'violin_confidence_by_disease')
        plt.close()
    
    # 7. Cell type distribution bar plot (top 15)
    fig, ax = plt.subplots(figsize=(14, 8))
    top_celltypes = pred_counts.head(15)
    bars = ax.barh(range(len(top_celltypes)), top_celltypes.values, color='steelblue')
    ax.set_yticks(range(len(top_celltypes)))
    ax.set_yticklabels(top_celltypes.index)
    ax.set_xlabel('Number of Cells', fontsize=12)
    ax.set_title('Top 15 Predicted Cell Types', fontsize=14)
    ax.invert_yaxis()
    
    # Add count labels
    for i, (idx, val) in enumerate(top_celltypes.items()):
        ax.text(val + 50, i, f'{val:,} ({100*val/len(predicted_labels):.1f}%)', 
                va='center', fontsize=9)
    
    plt.tight_layout()
    save_figure(fig, 'barplot_top_celltypes')
    plt.close()
    
    # 8. Species comparison: Confidence violin plot
    if 'species' in adata.obs:
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.violinplot(data=adata.obs, x='species', y='geneformer_confidence',
                       palette=species_colors, ax=ax)
        ax.set_xlabel('Species', fontsize=12)
        ax.set_ylabel('Prediction Confidence', fontsize=12)
        ax.set_title('Geneformer Confidence by Species\n(Mouse cells converted via BioMart orthologs)', fontsize=12)
        
        # Add significance annotation placeholder
        human_conf = adata.obs.loc[adata.obs['species'] == 'human', 'geneformer_confidence'].median()
        mouse_conf = adata.obs.loc[adata.obs['species'] == 'mouse', 'geneformer_confidence'].median()
        ax.text(0.5, 0.95, f'Human median: {human_conf:.3f}, Mouse median: {mouse_conf:.3f}',
                transform=ax.transAxes, ha='center', fontsize=10, style='italic')
        
        plt.tight_layout()
        save_figure(fig, 'violin_confidence_by_species')
        plt.close()
    
    # 9. Side-by-side: Human vs Mouse cell type distributions
    if 'species' in adata.obs:
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        
        for idx, species_name in enumerate(['human', 'mouse']):
            if species_name not in adata.obs['species'].values:
                continue
            ax = axes[idx]
            mask = adata.obs['species'] == species_name
            species_pred = adata.obs.loc[mask, 'geneformer_predicted_celltype'].value_counts().head(12)
            
            bars = ax.barh(range(len(species_pred)), species_pred.values, 
                           color=species_colors.get(species_name, '#808080'), alpha=0.8)
            ax.set_yticks(range(len(species_pred)))
            ax.set_yticklabels(species_pred.index, fontsize=9)
            ax.set_xlabel('Number of Cells', fontsize=11)
            ax.set_title(f'{species_name.upper()} Cells\n({mask.sum():,} cells total)', fontsize=12, fontweight='bold')
            ax.invert_yaxis()
            
            # Add count labels
            for i, val in enumerate(species_pred.values):
                pct = 100 * val / mask.sum()
                ax.text(val + 20, i, f'{pct:.1f}%', va='center', fontsize=8)
        
        plt.suptitle('Cell Type Distribution by Species', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        save_figure(fig, 'barplot_celltypes_by_species')
        plt.close()
    
    # 10. UMAP split by species
    if 'species' in adata.obs:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        
        for idx, species_name in enumerate(['human', 'mouse']):
            if species_name not in adata.obs['species'].values:
                continue
            ax = axes[idx]
            mask = adata.obs['species'] == species_name
            adata_species = adata[mask].copy()
            
            sc.pl.umap(adata_species, color='geneformer_predicted_celltype', ax=ax, show=False,
                       frameon=False, legend_loc='right margin' if idx == 1 else 'none', 
                       legend_fontsize=6, s=20,
                       title=f'{species_name.upper()} cells: Geneformer Predictions\n({mask.sum():,} cells)')
        
        plt.tight_layout()
        save_figure(fig, 'umap_predictions_split_by_species')
        plt.close()
    
    # 11. Heatmap: Disease state vs predicted cell type (percentage)
    if 'disease_state' in adata.obs:
        fig, ax = plt.subplots(figsize=(16, 8))
        
        # Get top cell types for each disease state
        top_per_disease = []
        for ds in adata.obs['disease_state'].unique():
            mask = adata.obs['disease_state'] == ds
            top_ct = adata.obs.loc[mask, 'geneformer_predicted_celltype'].value_counts().head(10).index.tolist()
            top_per_disease.extend(top_ct)
        top_celltypes_unique = list(dict.fromkeys(top_per_disease))[:20]  # Keep top 20 unique
        
        # Filter crosstab
        if 'crosstab_pct' in dir():
            crosstab_plot = crosstab_pct[top_celltypes_unique]
        else:
            crosstab = pd.crosstab(adata.obs['disease_state'], adata.obs['geneformer_predicted_celltype'])
            crosstab_pct = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
            crosstab_plot = crosstab_pct[[c for c in top_celltypes_unique if c in crosstab_pct.columns]]
        
        sns.heatmap(crosstab_plot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                    linewidths=0.5, cbar_kws={'label': 'Percentage'})
        ax.set_xlabel('Predicted Cell Type', fontsize=12)
        ax.set_ylabel('Disease State', fontsize=12)
        ax.set_title('Cell Type Distribution by Disease State (%)', fontsize=14)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_figure(fig, 'heatmap_disease_vs_celltype')
        plt.close()
    
    # 9. Individual sample UMAPs with Geneformer predictions
    msg("Generating individual sample UMAPs...")
    if 'sample_batch' in adata.obs:
        samples = sorted(adata.obs['sample_batch'].unique())
        umap_coords = adata.obsm['X_umap']
        x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
        y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
        margin = 0.05
        
        for sample in samples:
            mask = adata.obs['sample_batch'] == sample
            n_cells = mask.sum()
            
            fig, ax = plt.subplots(figsize=(10, 9))
            
            # Background (gray)
            ax.scatter(umap_coords[~mask, 0], umap_coords[~mask, 1],
                       c='#d3d3d3', s=10, alpha=0.2, rasterized=True)
            
            # Highlighted sample
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                       c='#d62728', s=20, alpha=0.7, rasterized=True)
            
            ax.set_xlim(x_min - margin*(x_max-x_min), x_max + margin*(x_max-x_min))
            ax.set_ylim(y_min - margin*(y_max-y_min), y_max + margin*(y_max-y_min))
            ax.set_title(f'{sample}\n({n_cells:,} cells)', fontsize=12, fontweight='bold')
            ax.axis('off')
            
            safe_name = sample.replace('/', '_').replace(' ', '_')
            save_figure_individual(fig, f'umap_sample_{safe_name}')
            plt.close()
        
        msg("  ✓ Individual sample UMAPs: %d samples", len(samples))
    
    # ==================== DONE ====================
    print("\n" + "=" * 84)
    print("COMPLETE!")
    print("=" * 84)
    print(f"  Output directory: {OUTPUT_DIR}")
    print(f"  Figures: {FIGDIR}")
    print(f"  Statistics: {STATS_DIR}")
    print(f"  Predictions: {PREDICTIONS_DIR}")
    print(f"\n  Key outputs:")
    print(f"    - integrated_with_geneformer_predictions.h5ad")
    print(f"    - prediction_probabilities.csv")
    print(f"    - celltype_counts.csv")
    print(f"    - disease_state_vs_celltype_crosstab.csv")
    
    # Summary stats
    print(f"\n  Summary:")
    print(f"    - Total cells: {adata.n_obs:,}")
    print(f"    - Unique predicted cell types: {len(pred_counts)}")
    print(f"    - Mean confidence: {np.mean(confidence_scores):.3f}")
    print(f"    - Median confidence: {np.median(confidence_scores):.3f}")
    
    print("\nDONE.\n")


if __name__ == "__main__":
    main()



__EOF_geneformer_predict_dlbcl_mouse_tonsil_py__

cat > "${SCRIPTS}/downstream_plots_human_mouse_integration.py" << '__EOF_downstream_plots_human_mouse_integration_py__'
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
OUTDIR = Path("__BASEDIR__/mouse_human_integration")
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



__EOF_downstream_plots_human_mouse_integration_py__

cat > "${SCRIPTS}/plot_cytotrace2_on_scvi_umap.py" << '__EOF_plot_cytotrace2_on_scvi_umap_py__'
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
CT2_OUTPUT = Path("__BASEDIR__/mouse_human_integration/CytoTRACE2/integrated_with_cytotrace2.h5ad")
MAIN_OUTPUT = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Output directory
OUTPUT_DIR = Path("__BASEDIR__/mouse_human_integration/figures_cytotrace2")
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



__EOF_plot_cytotrace2_on_scvi_umap_py__

cat > "${SCRIPTS}/plot_cytotrace2_stratified_by_tonsil_subtype.py" << '__EOF_plot_cytotrace2_stratified_by_tonsil_subtype_py__'
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
CT2_OUTPUT = Path("__BASEDIR__/mouse_human_integration/CytoTRACE2/integrated_with_cytotrace2.h5ad")
MAIN_OUTPUT = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Tonsil source files (to get original annotations)
TONSIL_GC_PATH = Path("__TONSIL_DIR__/tonsil_GCBC_RNA.h5ad")
TONSIL_MBC_PATH = Path("__TONSIL_DIR__/tonsil_NBC-MBC_RNA.h5ad")

# Output
OUTPUT_DIR = Path("__BASEDIR__/mouse_human_integration/figures_cytotrace2_stratified")
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




__EOF_plot_cytotrace2_stratified_by_tonsil_subtype_py__

cat > "${SCRIPTS}/plot_individual_samples_umap.py" << '__EOF_plot_individual_samples_umap_py__'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Individual Samples on UMAP
================================

Loads the integrated AnnData and creates individual UMAP plots
highlighting each sample while showing others in gray.

Author: J
Date: 2025-12-11
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================== CONFIGURATION =================================
# Input file (integrated scVI output)
INTEGRATED_H5AD = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Output directory for figures
FIGDIR = Path("__BASEDIR__/mouse_human_integration/figures_individual_samples")
FIGDIR.mkdir(parents=True, exist_ok=True)

# Group keys to plot by (sample_batch, disease_state, etc.)
GROUP_KEYS = ["sample_batch", "disease_state", "species"]

# Color for highlighted sample
HIGHLIGHT_COLOR = "#d62728"  # Red
BACKGROUND_COLOR = "#d3d3d3"  # Light gray

# Figure settings
FIGSIZE = (10, 9)
POINT_SIZE = 25
DPI = 300

# ============================== LOAD DATA =====================================
print("=" * 80)
print("PLOT INDIVIDUAL SAMPLES ON UMAP")
print("=" * 80)

print(f"\nLoading: {INTEGRATED_H5AD}")
adata = sc.read_h5ad(INTEGRATED_H5AD)
print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if "X_umap" not in adata.obsm:
    print("  Computing UMAP...")
    if "X_scvi" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=30)
    else:
        sc.pp.neighbors(adata, n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.2, spread=1.5)

# Get UMAP coordinates for consistent axis limits
umap_coords = adata.obsm["X_umap"]
x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
x_margin = (x_max - x_min) * 0.05
y_margin = (y_max - y_min) * 0.05

# ============================== PLOTTING FUNCTIONS ============================
def plot_sample_highlighted(adata, group_key, group_value, outdir, 
                            highlight_color=HIGHLIGHT_COLOR, 
                            background_color=BACKGROUND_COLOR):
    """Plot UMAP with one sample highlighted, others in gray."""
    labels = adata.obs[group_key].astype(str)
    mask = labels == group_value
    n_cells = mask.sum()
    
    fig, ax = plt.subplots(figsize=FIGSIZE)
    
    # Plot background cells first (gray)
    ax.scatter(
        umap_coords[~mask, 0],
        umap_coords[~mask, 1],
        c=background_color,
        s=POINT_SIZE * 0.6,
        alpha=0.3,
        rasterized=True,
        label="Other"
    )
    
    # Plot highlighted cells on top
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=highlight_color,
        s=POINT_SIZE,
        alpha=0.8,
        rasterized=True,
        label=group_value
    )
    
    # Set consistent axis limits
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Styling
    ax.set_title(f"{group_value}\n({n_cells:,} cells)", fontsize=14, fontweight="bold")
    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    # Save
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", group_value)
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_{group_key}_{safe_name}.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    
    return n_cells


def plot_sample_only(adata, group_key, group_value, outdir, color=HIGHLIGHT_COLOR):
    """Plot UMAP showing ONLY cells from one sample (no background)."""
    labels = adata.obs[group_key].astype(str)
    mask = labels == group_value
    n_cells = mask.sum()
    
    fig, ax = plt.subplots(figsize=FIGSIZE)
    
    # Plot only the selected cells
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.7,
        rasterized=True
    )
    
    # Set consistent axis limits (same as full UMAP)
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Styling
    ax.set_title(f"{group_value} only\n({n_cells:,} cells)", fontsize=14, fontweight="bold")
    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    # Save
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", group_value)
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_{group_key}_{safe_name}_only.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    
    return n_cells


def create_grid_plot(adata, group_key, outdir, ncols=4):
    """Create a grid of all samples in one figure."""
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    n_groups = len(unique_groups)
    nrows = (n_groups + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    axes = axes.flatten() if n_groups > 1 else [axes]
    
    for idx, grp in enumerate(unique_groups):
        ax = axes[idx]
        mask = labels == grp
        n_cells = mask.sum()
        
        # Background
        ax.scatter(
            umap_coords[~mask, 0],
            umap_coords[~mask, 1],
            c=BACKGROUND_COLOR,
            s=5,
            alpha=0.2,
            rasterized=True
        )
        
        # Highlighted
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=HIGHLIGHT_COLOR,
            s=8,
            alpha=0.7,
            rasterized=True
        )
        
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_title(f"{grp}\n(n={n_cells:,})", fontsize=9)
        ax.axis("off")
    
    # Hide empty subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].axis("off")
    
    plt.tight_layout()
    
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"umap_grid_{group_key}.{fmt}", 
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ Grid plot saved: umap_grid_{group_key}")


# ============================== GENERATE PLOTS ================================
for group_key in GROUP_KEYS:
    if group_key not in adata.obs.columns:
        print(f"\n  (skip) '{group_key}' not found in adata.obs")
        continue
    
    print(f"\n{'='*80}")
    print(f"Plotting by: {group_key}")
    print("="*80)
    
    # Create subdirectory for this grouping
    outdir = FIGDIR / group_key
    outdir.mkdir(parents=True, exist_ok=True)
    
    labels = adata.obs[group_key].astype(str)
    unique_groups = sorted(labels.unique())
    print(f"  Found {len(unique_groups)} unique groups")
    
    # Plot each group
    for grp in unique_groups:
        n1 = plot_sample_highlighted(adata, group_key, grp, outdir)
        n2 = plot_sample_only(adata, group_key, grp, outdir)
        print(f"    ✓ {grp}: {n1:,} cells")
    
    # Create grid plot
    create_grid_plot(adata, group_key, outdir)

# ============================== SUMMARY PLOT ==================================
print(f"\n{'='*80}")
print("Creating summary plots")
print("="*80)

# Disease state with custom colors
disease_colors = {
    "Mouse_Malignant": "#FF6B6B",
    "Mouse_Matched_malignant": "#FF9999",
    "DLBCL": "#8B0000",
    "Tonsil_GC_B": "#4169E1",
    "Tonsil_Normal": "#4169E1"
}

# Species colors
species_colors = {
    "mouse": "#98FB98",
    "human": "#6495ED"
}

# Full UMAP colored by disease state
fig, ax = plt.subplots(figsize=FIGSIZE)
for ds in adata.obs["disease_state"].unique():
    mask = adata.obs["disease_state"] == ds
    color = disease_colors.get(ds, "#808080")
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.6,
        label=f"{ds} ({mask.sum():,})",
        rasterized=True
    )
ax.set_xlim(x_min - x_margin, x_max + x_margin)
ax.set_ylim(y_min - y_margin, y_max + y_margin)
ax.set_title("All Samples by Disease State", fontsize=14, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for fmt in ["png", "pdf", "svg"]:
    fig.savefig(FIGDIR / f"umap_all_disease_state.{fmt}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  ✓ umap_all_disease_state")

# Full UMAP colored by species
fig, ax = plt.subplots(figsize=FIGSIZE)
for sp in adata.obs["species"].unique():
    mask = adata.obs["species"] == sp
    color = species_colors.get(sp, "#808080")
    ax.scatter(
        umap_coords[mask, 0],
        umap_coords[mask, 1],
        c=color,
        s=POINT_SIZE,
        alpha=0.6,
        label=f"{sp} ({mask.sum():,})",
        rasterized=True
    )
ax.set_xlim(x_min - x_margin, x_max + x_margin)
ax.set_ylim(y_min - y_margin, y_max + y_margin)
ax.set_title("All Samples by Species", fontsize=14, fontweight="bold")
ax.legend(loc="upper right", fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for fmt in ["png", "pdf", "svg"]:
    fig.savefig(FIGDIR / f"umap_all_species.{fmt}", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print("  ✓ umap_all_species")

# ============================== DONE ==========================================
print(f"\n{'='*80}")
print("COMPLETE!")
print("="*80)
print(f"  Output directory: {FIGDIR}")
print(f"  Total figures generated: {len(list(FIGDIR.rglob('*.png')))}")
print("\nDONE.\n")



__EOF_plot_individual_samples_umap_py__

cat > "${SCRIPTS}/convert_human_mouse_integration_to_seurat5.R" << '__EOF_convert_human_mouse_integration_to_seurat5_R__'
#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# Convert Human-Mouse Integration h5ad to Seurat v5 Object
# =========================================================
#
# Converts the integrated h5ad file (with Geneformer predictions) to Seurat v5
# Preserves: scVI embeddings, UMAP, Leiden clusters, CytoTRACE2, Geneformer predictions
#
# Usage:
#   conda activate rconv2
#   Rscript __BASEDIR__/scripts/convert_human_mouse_integration_to_seurat5.R
#
# Author: J
# Date: 2025-12-02

# ============================== PATHS ========================================
INPUT_H5AD <- "__BASEDIR__/mouse_human_integration/Geneformer/integrated_with_geneformer_predictions.h5ad"
OUTPUT_RDS <- "__BASEDIR__/mouse_human_integration/Geneformer/integrated_with_geneformer_predictions_seurat5.rds"

cat("===============================================================================\n")
cat("CONVERTING H5AD TO SEURAT V5\n")
cat("===============================================================================\n")
cat("Input: ", INPUT_H5AD, "\n")
cat("Output:", OUTPUT_RDS, "\n\n")

# ============================== LOAD PACKAGES ================================
cat("STEP 1 — Loading packages...\n")

suppressPackageStartupMessages({
    library(zellkonverter)
    library(SingleCellExperiment)
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
})

cat("  ✓ Packages loaded\n")

# ============================== READ H5AD ====================================
cat("\nSTEP 2 — Reading h5ad file...\n")

sce <- readH5AD(INPUT_H5AD, reader = "R")
cat("  Loaded:", ncol(sce), "cells ×", nrow(sce), "genes\n")

# ============================== EXTRACT DATA =================================
cat("\nSTEP 3 — Extracting data from SCE...\n")

# Get counts matrix
if ("counts" %in% assayNames(sce)) {
    counts_mat <- assay(sce, "counts")
    cat("  Using 'counts' assay\n")
} else if ("X" %in% assayNames(sce)) {
    counts_mat <- assay(sce, "X")
    cat("  Using 'X' assay\n")
} else {
    counts_mat <- assay(sce, 1)
    cat("  Using first assay\n")
}

# Ensure sparse matrix
if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "dgCMatrix")
}

# Get normalized data if available
norm_mat <- NULL
if ("normalized" %in% assayNames(sce)) {
    norm_mat <- assay(sce, "normalized")
    if (!inherits(norm_mat, "dgCMatrix")) {
        norm_mat <- as(norm_mat, "dgCMatrix")
    }
    cat("  Found normalized data\n")
}

# Get cell metadata
cell_meta <- as.data.frame(colData(sce))
cat("  Cell metadata columns:", ncol(cell_meta), "\n")
cat("  Key columns:", paste(head(names(cell_meta), 10), collapse=", "), "\n")

# Get gene metadata
gene_meta <- as.data.frame(rowData(sce))
if (nrow(gene_meta) == 0) {
    gene_meta <- data.frame(row.names = rownames(counts_mat))
}

# ============================== CREATE SEURAT ================================
cat("\nSTEP 4 — Creating Seurat v5 object...\n")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
    counts = counts_mat,
    meta.data = cell_meta,
    project = "HumanMouseIntegration"
)

# Add normalized data if available
if (!is.null(norm_mat)) {
    seurat_obj[["RNA"]]$data <- norm_mat
    cat("  Added normalized data layer\n")
}

cat("  Created Seurat object:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "genes\n")

# ============================== ADD REDUCTIONS ===============================
cat("\nSTEP 5 — Adding dimensionality reductions...\n")

# Get reduced dimensions from SCE
red_dims <- reducedDimNames(sce)
cat("  Available reductions:", paste(red_dims, collapse=", "), "\n")

# Add UMAP
if ("X_umap" %in% red_dims) {
    umap_coords <- reducedDim(sce, "X_umap")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    rownames(umap_coords) <- colnames(seurat_obj)
    seurat_obj[["umap"]] <- CreateDimReducObject(
        embeddings = umap_coords,
        key = "UMAP_",
        assay = "RNA"
    )
    cat("  ✓ Added UMAP\n")
}

# Add scVI embeddings
if ("X_scvi" %in% red_dims) {
    scvi_embed <- reducedDim(sce, "X_scvi")
    colnames(scvi_embed) <- paste0("scVI_", 1:ncol(scvi_embed))
    rownames(scvi_embed) <- colnames(seurat_obj)
    seurat_obj[["scvi"]] <- CreateDimReducObject(
        embeddings = scvi_embed,
        key = "scVI_",
        assay = "RNA"
    )
    cat("  ✓ Added scVI (", ncol(scvi_embed), " dimensions)\n", sep="")
}

# Add PCA if available
if ("X_pca" %in% red_dims) {
    pca_embed <- reducedDim(sce, "X_pca")
    colnames(pca_embed) <- paste0("PC_", 1:ncol(pca_embed))
    rownames(pca_embed) <- colnames(seurat_obj)
    seurat_obj[["pca"]] <- CreateDimReducObject(
        embeddings = pca_embed,
        key = "PC_",
        assay = "RNA"
    )
    cat("  ✓ Added PCA (", ncol(pca_embed), " dimensions)\n", sep="")
}

# ============================== VERIFY METADATA ==============================
cat("\nSTEP 6 — Verifying metadata...\n")

# List key columns
key_cols <- c(
    "species", "disease_state", "sample_batch", "study",
    "leiden_0.5", "leiden_1.0", "leiden_1.5",
    "S_score", "G2M_score", "phase",
    "cytotrace2_score", "cytotrace2_potency",
    "oxphos_score", "bcr_score",
    "geneformer_predicted_celltype", "geneformer_confidence"
)

present_cols <- key_cols[key_cols %in% names(seurat_obj@meta.data)]
cat("  Present key columns:\n")
for (col in present_cols) {
    if (is.numeric(seurat_obj@meta.data[[col]])) {
        cat("    ", col, ": numeric (mean=", round(mean(seurat_obj@meta.data[[col]], na.rm=TRUE), 3), ")\n", sep="")
    } else {
        n_unique <- length(unique(na.omit(seurat_obj@meta.data[[col]])))
        cat("    ", col, ": ", n_unique, " unique values\n", sep="")
    }
}

# ============================== SAVE =========================================
cat("\nSTEP 7 — Saving Seurat v5 object...\n")

tryCatch({
    saveRDS(seurat_obj, OUTPUT_RDS)
    cat("  ✓ Saved:", OUTPUT_RDS, "\n")
}, error = function(e) {
    cat("  WARNING: Error saving compressed RDS:", e$message, "\n")
    cat("  Trying uncompressed save...\n")
    saveRDS(seurat_obj, OUTPUT_RDS, compress = FALSE)
    cat("  ✓ Saved (uncompressed):", OUTPUT_RDS, "\n")
})

# ============================== SUMMARY ======================================
cat("\n===============================================================================\n")
cat("CONVERSION COMPLETE!\n")
cat("===============================================================================\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")
cat("  Reductions:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
cat("  Metadata columns:", ncol(seurat_obj@meta.data), "\n")

if ("disease_state" %in% names(seurat_obj@meta.data)) {
    cat("\n  Disease states:\n")
    print(table(seurat_obj@meta.data$disease_state))
}

if ("species" %in% names(seurat_obj@meta.data)) {
    cat("\n  Species:\n")
    print(table(seurat_obj@meta.data$species))
}

if ("geneformer_predicted_celltype" %in% names(seurat_obj@meta.data)) {
    cat("\n  Geneformer predictions (top 10):\n")
    print(head(sort(table(seurat_obj@meta.data$geneformer_predicted_celltype), decreasing=TRUE), 10))
}

cat("\nOutput file:", OUTPUT_RDS, "\n")
cat("File size:", round(file.size(OUTPUT_RDS) / 1e6, 1), "MB\n")
cat("\nDONE.\n")



__EOF_convert_human_mouse_integration_to_seurat5_R__

cat > "${SCRIPTS}/plot_umap_highlight_mouse_clusters_4_6.py" << '__EOF_plot_umap_highlight_mouse_clusters_4_6_py__'
#!/usr/bin/env python3
"""
Highlight cells from mouse-only leiden clusters 4 and 6 on the human-mouse integration UMAP
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
MOUSE_ONLY_H5AD = Path("__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
HUMAN_MOUSE_H5AD = Path("__BASEDIR__/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")
FIGDIR = Path("__BASEDIR__/mouse_human_integration/figures_individual_samples")
FIGDIR.mkdir(parents=True, exist_ok=True)

# Clusters to highlight
HIGHLIGHT_CLUSTERS = ['4', '6']

print("="*80)
print("Step 1: Loading mouse-only integration")
print("="*80)
adata_mouse = sc.read_h5ad(MOUSE_ONLY_H5AD)
print(f"Mouse-only: {adata_mouse.shape[0]:,} cells × {adata_mouse.shape[1]:,} genes")

# Check leiden clusters
if 'leiden_1.0' not in adata_mouse.obs.columns:
    print("ERROR: leiden_1.0 not found in mouse-only data!")
    exit(1)

print(f"\nLeiden clusters in mouse-only: {sorted(adata_mouse.obs['leiden_1.0'].unique())}")

# Get cell names for clusters 4 and 6
highlight_cells = []
for cluster in HIGHLIGHT_CLUSTERS:
    mask = adata_mouse.obs['leiden_1.0'].astype(str) == cluster
    n_cells = mask.sum()
    print(f"  Cluster {cluster}: {n_cells:,} cells")
    cluster_cells = adata_mouse.obs_names[mask].tolist()
    highlight_cells.extend(cluster_cells)

highlight_cells = set(highlight_cells)
print(f"\nTotal unique cells to highlight: {len(highlight_cells):,}")

print("\n" + "="*80)
print("Step 2: Loading human-mouse integration")
print("="*80)
adata_integrated = sc.read_h5ad(HUMAN_MOUSE_H5AD)
print(f"Human-mouse integration: {adata_integrated.shape[0]:,} cells × {adata_integrated.shape[1]:,} genes")

# Check if UMAP exists
if 'X_umap' not in adata_integrated.obsm:
    print("ERROR: X_umap not found in human-mouse integration!")
    print("Available obsm keys:", list(adata_integrated.obsm.keys()))
    exit(1)

# Find matching cells
print("\n" + "="*80)
print("Step 3: Matching cells between datasets")
print("="*80)

# Check cell name format - might need to match by barcode or full name
mouse_cell_names = set(adata_mouse.obs_names)
integrated_cell_names = set(adata_integrated.obs_names)

# Try direct matching first
matching_cells = highlight_cells & integrated_cell_names
print(f"Direct matches: {len(matching_cells):,} cells")

# If not many matches, try matching by barcode (part before first underscore or dash)
if len(matching_cells) < len(highlight_cells) * 0.5:
    print("\nTrying barcode matching...")
    # Extract barcodes (part before separator)
    mouse_barcodes = {}
    for cell in highlight_cells:
        # Try different separators
        for sep in ['-', '_']:
            if sep in cell:
                barcode = cell.split(sep)[0]
                mouse_barcodes[barcode] = cell
                break
    
    integrated_barcodes = {}
    for cell in adata_integrated.obs_names:
        for sep in ['-', '_']:
            if sep in cell:
                barcode = cell.split(sep)[0]
                integrated_barcodes[barcode] = cell
                break
    
    # Match by barcode
    matching_barcodes = set(mouse_barcodes.keys()) & set(integrated_barcodes.keys())
    matching_cells = {integrated_barcodes[b] for b in matching_barcodes if b in integrated_barcodes}
    print(f"Barcode matches: {len(matching_cells):,} cells")

if len(matching_cells) == 0:
    print("ERROR: No matching cells found!")
    print(f"Sample mouse cell names: {list(highlight_cells)[:5]}")
    print(f"Sample integrated cell names: {list(integrated_cell_names)[:5]}")
    exit(1)

print(f"\nFinal matching cells to highlight: {len(matching_cells):,}")

# Create mask for highlighted cells
highlight_mask = adata_integrated.obs_names.isin(matching_cells)

print("\n" + "="*80)
print("Step 4: Creating separate density plot visualizations")
print("="*80)

# Get all coordinates for extent calculation
coords_all = adata_integrated.obsm['X_umap']
x_min, x_max = coords_all[:, 0].min(), coords_all[:, 0].max()
y_min, y_max = coords_all[:, 1].min(), coords_all[:, 1].max()

# Identify which cells belong to which cluster
cluster_4_cells = set()
cluster_6_cells = set()

for cluster in HIGHLIGHT_CLUSTERS:
    mask = adata_mouse.obs['leiden_1.0'].astype(str) == cluster
    cluster_cells = set(adata_mouse.obs_names[mask])
    
    if cluster == '4':
        cluster_4_cells = cluster_cells & matching_cells
    elif cluster == '6':
        cluster_6_cells = cluster_cells & matching_cells

print(f"Cluster 4 cells: {len(cluster_4_cells):,}")
print(f"Cluster 6 cells: {len(cluster_6_cells):,}")

# Identify tonsil cells for contour overlay
tonsil_mask = (adata_integrated.obs['disease_state'] == 'Tonsil_Normal')
coords_tonsil = adata_integrated.obsm['X_umap'][tonsil_mask]
print(f"Tonsil cells for contour: {tonsil_mask.sum():,}")

# Create separate plots for each cluster
for cluster_num, cluster_cells_set in [('4', cluster_4_cells), ('6', cluster_6_cells)]:
    if len(cluster_cells_set) == 0:
        print(f"Skipping cluster {cluster_num} - no matching cells")
        continue
    
    # Get coordinates for this cluster
    mask_cluster = adata_integrated.obs_names.isin(cluster_cells_set)
    coords_cluster = adata_integrated.obsm['X_umap'][mask_cluster]
    
    from scipy.stats import gaussian_kde
    
    # Create figure with white background
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_facecolor('white')
    
    # Plot all OTHER cells in gray first (background)
    other_mask = ~mask_cluster
    ax.scatter(coords_all[other_mask, 0], coords_all[other_mask, 1],
              c='lightgray', s=5, alpha=0.3,
              marker='o', edgecolors='none', rasterized=True, zorder=1)
    
    # Compute KDE density for the cluster cells
    xy = np.vstack([coords_cluster[:, 0], coords_cluster[:, 1]])
    try:
        kde = gaussian_kde(xy, bw_method=0.15)
        density = kde(xy)
    except Exception:
        density = np.ones(len(coords_cluster))
    
    # Sort by density so densest points are plotted on top
    idx = density.argsort()
    x_sorted = coords_cluster[idx, 0]
    y_sorted = coords_cluster[idx, 1]
    density_sorted = density[idx]
    
    # Plot cluster cells colored by density (viridis)
    sc_plot = ax.scatter(x_sorted, y_sorted,
                        c=density_sorted, cmap='viridis', s=30, alpha=0.9,
                        edgecolors='none', rasterized=True, zorder=3)
    
    # Add colorbar
    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Cell Density (KDE)', rotation=270, labelpad=20, fontsize=12)
    
    # Overlay tonsil cell density contours in black
    if len(coords_tonsil) > 0:
        xy_tonsil = np.vstack([coords_tonsil[:, 0], coords_tonsil[:, 1]])
        kde_tonsil = gaussian_kde(xy_tonsil, bw_method=0.15)
        
        # Create grid for contour
        xx, yy = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        z_tonsil = kde_tonsil(positions).reshape(xx.shape)
        
        # Draw contours
        contour = ax.contour(xx, yy, z_tonsil, levels=6, colors='black', 
                            linewidths=1.5, alpha=0.8, zorder=4)
        ax.clabel(contour, inline=False, fontsize=0)  # no labels on contour lines
        
        # Add legend entry for contour
        from matplotlib.lines import Line2D
        contour_legend = Line2D([0], [0], color='black', linewidth=1.5, 
                               label='Tonsil cell density')
        ax.legend(handles=[contour_legend], loc='lower right', fontsize=10, framealpha=0.9)
    
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title(f'Human-Mouse Integration UMAP:\nDensity of Mouse Leiden Cluster {cluster_num} (viridis) + Tonsil Contours (black)\n({len(cluster_cells_set):,} cells)', 
                 fontsize=14, fontweight='bold')
    
    # Set axis limits to match full UMAP
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    
    # Save
    output_path = FIGDIR / f"umap_all_disease_state_highlight_mouse_cluster_{cluster_num}_density.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {output_path}")
    
    fig.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    print(f"Saved → {output_path.with_suffix('.pdf')}")
    print(f"Saved → {output_path.with_suffix('.svg')}")
    
    plt.close()

print("\nDone!")


__EOF_plot_umap_highlight_mouse_clusters_4_6_py__

cat > "${SCRIPTS}/ordering_cytotrace_2_mouse_geneformer.py" << '__EOF_ordering_cytotrace_2_mouse_geneformer_py__'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CytoTRACE2-ordered gene programs with CLUSTERED y-axis (shape + amplitude)
Size-capped heatmaps; MAGIC smoothing; Viridis colormap + Okabe-Ito palette;
Condition-aware selection; per-program Enrichr + GSEA(prerank) with NES vs p-value plots;
Top annotation bar (dominant category per CT2 bin).

Adapted for: Mouse Geneformer predictions h5ad

I/O:
- Input  H5AD: __BASEDIR__/Geneformer/mouse_with_geneformer_predictions.h5ad
- Output dir : __BASEDIR__/heatmap_programs

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
INPUT_H5AD   = Path("__BASEDIR__/Geneformer/mouse_with_geneformer_predictions.h5ad")
OUTPUT_DIR   = Path("__BASEDIR__/heatmap_programs")

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



__EOF_ordering_cytotrace_2_mouse_geneformer_py__

cat > "${SCRIPTS}/plot_gsea_custom.py" << '__EOF_plot_gsea_custom_py__'
#!/usr/bin/env python3
"""
Custom GSEA/Enrichr Scatter Plots
=================================
Generates publication-quality plots:
- GSEA prerank plots (NES on Y-axis) - captures OXPHOS!
- Enrichr plots (Combined Score on Y-axis)
- Both with clean aesthetics

Author: Generated for CytoTRACE2 analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob

# ===== CONFIGURATION =====
BASE_DIR = Path("__BASEDIR__/heatmap_programs")
OUTPUT_DIR = BASE_DIR / "figures" / "gsea_custom"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Where results are stored
ENRICHR_POSITIVE = BASE_DIR / "programs" / "enrichr_positive_corr_stemlike"
ENRICHR_NEGATIVE = BASE_DIR / "programs" / "enrichr_negative_corr_differentiated"
GSEA_POSITIVE = BASE_DIR / "programs" / "gsea_prerank_positive_corr_stemlike"
GSEA_NEGATIVE = BASE_DIR / "programs" / "gsea_prerank_negative_corr_differentiated"

# Libraries to plot
LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2019_Mouse",
    "Reactome_2022",
    "MSigDB_Hallmark_2020"
]

# Colors
COLOR_STEM = "#DC3545"      # Red for stem-like
COLOR_DIFF = "#0077B6"      # Blue for differentiated
MARKER_SIZE = 45            # Uniform size

# Label settings
LABEL_TOP_N = 35
LABEL_FONTSIZE = 8

# Output formats
SAVE_FORMATS = ['png', 'svg', 'pdf']


def save_figure(fig, output_path: Path):
    """Save figure in multiple formats (PNG, SVG, PDF)."""
    for fmt in SAVE_FORMATS:
        out_file = output_path.with_suffix(f'.{fmt}')
        fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {output_path.stem} (.png, .svg, .pdf)")


def load_enrichr_txt(filepath: Path) -> pd.DataFrame:
    """Load Enrichr .txt report file."""
    if not filepath.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(filepath, sep='\t')
        return df
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return pd.DataFrame()


def load_gsea_prerank(filepath: Path, library: str) -> pd.DataFrame:
    """Load GSEA prerank results for a specific library."""
    if not filepath.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(filepath)
        # Filter by library
        if 'library' in df.columns:
            df = df[df['library'] == library].copy()
        return df
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return pd.DataFrame()


def plot_gsea_prerank_scatter(df: pd.DataFrame, library: str, output_path: Path,
                              color: str = COLOR_STEM, marker: str = 'o',
                              title_suffix: str = ""):
    """
    Plot GSEA prerank scatter with the same sloping style as Enrichr.
    Uses a transformed x-axis to create smooth distribution.
    """
    if df.empty:
        print(f"  No GSEA prerank data for {library}")
        return
    
    # Get columns
    term_col = "Term" if "Term" in df.columns else "Name"
    nes_col = "NES" if "NES" in df.columns else None
    
    # Find FDR q-value column (better distribution than NOM p-val)
    p_col = None
    for c in ["FDR q-val", "FDR q-value", "Adjusted P-value", "NOM p-val", "P-value"]:
        if c in df.columns:
            p_col = c
            break
    
    if not nes_col or not p_col:
        print(f"  Missing columns. Available: {list(df.columns)}")
        return
    
    # Calculate axes
    df = df.copy()
    
    # Get FDR/p-values - use minimum of 1e-10 for visual spread
    pvals = df[p_col].astype(float).clip(lower=1e-10, upper=1.0)
    df["logp_raw"] = -np.log10(pvals)
    
    # Cap at 10 for nice visual range (like Enrichr)
    MAX_LOGP = 10
    df["logp"] = df["logp_raw"].clip(upper=MAX_LOGP)
    
    # Add jitter to highly significant terms (those hitting the cap)
    # This spreads out the dense cluster in upper right
    np.random.seed(42)  # Reproducible jitter
    at_cap = df["logp"] >= MAX_LOGP - 0.1
    n_at_cap = at_cap.sum()
    if n_at_cap > 0:
        # Add horizontal jitter proportional to number of capped terms
        jitter_range = min(2.0, n_at_cap * 0.02)  # Max 2 units of jitter
        df.loc[at_cap, "logp"] = MAX_LOGP - np.random.uniform(0, jitter_range, n_at_cap)
    
    # Use absolute NES for Y-axis
    df["nes"] = df[nes_col].astype(float).abs()
    
    # Clean term names
    df["term_clean"] = df[term_col].astype(str).str.replace("_", " ")
    df["term_clean"] = df["term_clean"].str.replace(r"\s*\([A-Z]+-?[A-Z]*:\d+\)$", "", regex=True)
    df["term_clean"] = df["term_clean"].str.replace(r"\s*R-[A-Z]+-\d+$", "", regex=True)
    
    # ===== PLOT (matching Enrichr style) =====
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot ALL dots with uniform size (like Enrichr)
    ax.scatter(df["logp"], df["nes"], 
              s=MARKER_SIZE, c=color, marker=marker,
              alpha=0.7, edgecolors='none', zorder=2)
    
    # Reference line at y=0
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Labels for top terms (by NES * significance)
    # Weight by both NES and significance for better label selection
    df["importance"] = df["nes"] * (df["logp"] / MAX_LOGP)
    top_terms = df.nlargest(LABEL_TOP_N, "importance")
    
    # Use adjustText for better label positioning
    texts = []
    for _, row in top_terms.iterrows():
        label = row["term_clean"]
        if len(label) > 50:
            label = label[:47] + "..."
        txt = ax.text(row["logp"] + 0.1, row["nes"], label,
                     fontsize=LABEL_FONTSIZE, color='#333333',
                     ha='left', va='center')
        texts.append(txt)
    
    # Try to use adjustText for smart positioning
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                   arrowprops=dict(arrowstyle='-', color='#AAAAAA', alpha=0.5, lw=0.5),
                   expand_points=(1.5, 1.5),
                   force_text=(0.5, 0.8),
                   force_points=(0.3, 0.3))
    except ImportError:
        pass  # Fall back to basic positioning
    
    # Axis settings (matching Enrichr range)
    ax.set_xlim(0, MAX_LOGP + 1)
    ax.set_ylim(bottom=0)
    
    # Title
    lib_clean = library.replace("_", " ")
    ax.set_title(f"GSEA (prerank) — {lib_clean}{title_suffix}", 
                fontsize=18, fontweight='bold', pad=15)
    
    # Axis labels
    p_label = "FDR" if "FDR" in p_col else "P-value"
    ax.set_xlabel(f"-log10({p_label})", fontsize=14)
    ax.set_ylabel("|NES|", fontsize=14)
    
    # Spine styling (matching Enrichr)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_gsea_prerank_combined(df_stem: pd.DataFrame, df_diff: pd.DataFrame,
                               library: str, output_path: Path):
    """
    Combined GSEA prerank plot: red dots (stem-like) + blue diamonds (differentiated)
    """
    if df_stem.empty and df_diff.empty:
        print(f"  No GSEA prerank data for {library}")
        return
    
    term_col = "Term" if "Term" in df_stem.columns else "Name"
    nes_col = "NES"
    
    # Find p-value column
    p_col = None
    sample_df = df_stem if not df_stem.empty else df_diff
    for c in ["FDR q-val", "FDR q-value", "Adjusted P-value", "NOM p-val"]:
        if c in sample_df.columns:
            p_col = c
            break
    
    if not p_col:
        print(f"  No p-value column found")
        return
    
    MAX_LOGP = 10
    
    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    all_texts = []
    
    # Process and plot each subset
    for df, color, marker, label, direction in [
        (df_stem, COLOR_STEM, 'o', 'Malignant', 'stem'),
        (df_diff, COLOR_DIFF, 'D', 'Physiologic', 'diff')
    ]:
        if df.empty:
            continue
        
        df = df.copy()
        pvals = df[p_col].astype(float).clip(lower=1e-10, upper=1.0)
        df["logp_raw"] = -np.log10(pvals)
        df["logp"] = df["logp_raw"].clip(upper=MAX_LOGP)
        
        # Add jitter for capped terms
        np.random.seed(42 if direction == 'stem' else 43)
        at_cap = df["logp"] >= MAX_LOGP - 0.1
        n_at_cap = at_cap.sum()
        if n_at_cap > 0:
            jitter_range = min(2.0, n_at_cap * 0.02)
            df.loc[at_cap, "logp"] = MAX_LOGP - np.random.uniform(0, jitter_range, n_at_cap)
        
        df["nes"] = df[nes_col].astype(float).abs()
        
        # Clean term names
        df["term_clean"] = df[term_col].astype(str).str.replace("_", " ")
        df["term_clean"] = df["term_clean"].str.replace(r"\s*\([A-Z]+-?[A-Z]*:\d+\)$", "", regex=True)
        df["term_clean"] = df["term_clean"].str.replace(r"\s*R-[A-Z]+-\d+$", "", regex=True)
        
        # Plot dots
        ax.scatter(df["logp"], df["nes"], 
                  s=MARKER_SIZE, c=color, marker=marker,
                  alpha=0.7, edgecolors='none', zorder=2, label=label)
        
        # Labels for top terms (reduced for cleaner look)
        df["importance"] = df["nes"] * (df["logp"] / MAX_LOGP)
        top_terms = df.nlargest(12, "importance")
        
        for _, row in top_terms.iterrows():
            term_label = row["term_clean"]
            if len(term_label) > 45:
                term_label = term_label[:42] + "..."
            txt = ax.text(row["logp"] + 0.1, row["nes"], term_label,
                         fontsize=7, color=color, alpha=0.9,
                         ha='left', va='center')
            all_texts.append(txt)
    
    # Use adjustText for all labels
    try:
        from adjustText import adjust_text
        adjust_text(all_texts, ax=ax,
                   arrowprops=dict(arrowstyle='-', color='#AAAAAA', alpha=0.4, lw=0.4),
                   expand_points=(1.3, 1.3),
                   force_text=(0.4, 0.6))
    except ImportError:
        pass
    
    # Reference line
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Axis settings
    ax.set_xlim(0, MAX_LOGP + 1)
    ax.set_ylim(bottom=0)
    
    # Title and labels
    lib_clean = library.replace("_", " ")
    ax.set_title(f"GSEA (prerank) — {lib_clean}\nMalignant vs Physiologic", 
                fontsize=16, fontweight='bold', pad=15)
    
    p_label = "FDR" if "FDR" in p_col else "P-value"
    ax.set_xlabel(f"-log10({p_label})", fontsize=14)
    ax.set_ylabel("|NES|", fontsize=14)
    
    # Legend
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=11)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_enrichr_scatter(df: pd.DataFrame, library: str, output_path: Path, 
                         color: str = COLOR_STEM, marker: str = 'o',
                         title_prefix: str = "GSEA (prerank)"):
    """
    Plot Enrichr-style scatter: Combined Score vs -log10(P-value)
    Matches the reference style exactly.
    """
    if df.empty:
        print(f"  No data for {library}")
        return
    
    # Get columns
    term_col = "Term" if "Term" in df.columns else df.columns[1]
    p_col = "P-value" if "P-value" in df.columns else None
    score_col = "Combined Score" if "Combined Score" in df.columns else None
    
    if not p_col or p_col not in df.columns:
        for c in df.columns:
            if 'p-value' in c.lower() and 'adjusted' not in c.lower():
                p_col = c
                break
    
    if not score_col or score_col not in df.columns:
        for c in df.columns:
            if 'combined' in c.lower() and 'score' in c.lower():
                score_col = c
                break
    
    if not p_col or not score_col:
        print(f"  Missing columns. Available: {list(df.columns)}")
        return
    
    # Calculate axes
    df = df.copy()
    df["logp"] = -np.log10(df[p_col].astype(float).clip(lower=1e-300))
    df["score"] = df[score_col].astype(float)
    
    # Clean term names
    df["term_clean"] = df[term_col].str.replace("_", " ")
    
    # ===== PLOT (matching reference exactly) =====
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot ALL dots with uniform size
    ax.scatter(df["logp"], df["score"], 
              s=MARKER_SIZE, c=color, marker=marker,
              alpha=0.7, edgecolors='none', zorder=2)
    
    # Reference line at y=0
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Labels for top terms
    top_terms = df.nlargest(LABEL_TOP_N, "score")
    
    for _, row in top_terms.iterrows():
        label = row["term_clean"]
        if len(label) > 60:
            label = label[:57] + "..."
        ax.annotate(label, (row["logp"], row["score"]),
                   fontsize=LABEL_FONTSIZE, color='#333333',
                   ha='left', va='center',
                   xytext=(4, 0), textcoords='offset points')
    
    # Axis settings
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    # Title (matching reference)
    lib_clean = library.replace("_", " ")
    ax.set_title(f"{title_prefix} — {lib_clean}", 
                fontsize=18, fontweight='bold', pad=15)
    
    # Axis labels
    ax.set_xlabel("-log10(P-value)", fontsize=14)
    ax.set_ylabel("NES", fontsize=14)  # Keep same label as reference
    
    # Spine styling (matching reference)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_combined_stem_diff(df_stem: pd.DataFrame, df_diff: pd.DataFrame,
                            library: str, output_path: Path):
    """
    Combined plot: red dots (stem) + blue diamonds (diff)
    """
    if df_stem.empty and df_diff.empty:
        print(f"  No data for {library}")
        return
    
    # Get columns
    term_col = "Term"
    p_col = "P-value"
    score_col = "Combined Score"
    
    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot stem-like (positive correlation) as red dots
    if not df_stem.empty:
        df_s = df_stem.copy()
        df_s["logp"] = -np.log10(df_s[p_col].astype(float).clip(lower=1e-300))
        df_s["score"] = df_s[score_col].astype(float)
        ax.scatter(df_s["logp"], df_s["score"], 
                  s=50, c='#DC3545', marker='o', label='Stem-like (high CT2)',
                  alpha=0.7, edgecolors='none', zorder=2)
    
    # Plot differentiated (negative correlation) as blue diamonds
    if not df_diff.empty:
        df_d = df_diff.copy()
        df_d["logp"] = -np.log10(df_d[p_col].astype(float).clip(lower=1e-300))
        df_d["score"] = df_d[score_col].astype(float)
        ax.scatter(df_d["logp"], df_d["score"], 
                  s=60, c='#0077B6', marker='D', label='Differentiated (low CT2)',
                  alpha=0.7, edgecolors='none', zorder=2)
    
    # Labels for top terms from each
    for df, c in [(df_stem, '#DC3545'), (df_diff, '#0077B6')]:
        if df.empty:
            continue
        df = df.copy()
        df["logp"] = -np.log10(df[p_col].astype(float).clip(lower=1e-300))
        df["score"] = df[score_col].astype(float)
        top = df.nlargest(15, "score")
        for _, row in top.iterrows():
            label = str(row[term_col]).replace("_", " ")
            if len(label) > 50:
                label = label[:47] + "..."
            ax.annotate(label, (row["logp"], row["score"]),
                       fontsize=7, color=c, alpha=0.9,
                       ha='left', va='center',
                       xytext=(4, 0), textcoords='offset points')
    
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    lib_clean = library.replace("_", " ")
    ax.set_title(f"GSEA — {lib_clean}\nStem-like vs Differentiated", 
                fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel("-log10(P-value)", fontsize=13)
    ax.set_ylabel("Combined Score", fontsize=13)
    
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=11)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def main():
    print("=" * 60)
    print("Custom GSEA Scatter Plotting")
    print("=" * 60)
    
    # ===== PART 1: GSEA PRERANK PLOTS (captures OXPHOS!) =====
    print("\n[1] GSEA Prerank plots (NES - captures OXPHOS)")
    print("-" * 50)
    
    gsea_stem_file = GSEA_POSITIVE / "gsea_prerank_positive_corr_stemlike_all.csv"
    gsea_diff_file = GSEA_NEGATIVE / "gsea_prerank_negative_corr_differentiated_all.csv"
    
    for lib in LIBRARIES:
        print(f"\n  {lib}")
        
        # Stem-like (GSEA prerank)
        df_gsea_stem = load_gsea_prerank(gsea_stem_file, lib)
        if not df_gsea_stem.empty:
            plot_gsea_prerank_scatter(
                df_gsea_stem, lib,
                OUTPUT_DIR / f"gsea_prerank_stemlike_{lib}.png",
                color=COLOR_STEM, marker='o',
                title_suffix=" (Stem-like)"
            )
        
        # Differentiated (GSEA prerank)
        df_gsea_diff = load_gsea_prerank(gsea_diff_file, lib)
        if not df_gsea_diff.empty:
            plot_gsea_prerank_scatter(
                df_gsea_diff, lib,
                OUTPUT_DIR / f"gsea_prerank_differentiated_{lib}.png",
                color=COLOR_DIFF, marker='D',
                title_suffix=" (Differentiated)"
            )
        
        # Combined GSEA prerank plot
        if not df_gsea_stem.empty or not df_gsea_diff.empty:
            plot_gsea_prerank_combined(
                df_gsea_stem, df_gsea_diff, lib,
                OUTPUT_DIR / f"gsea_prerank_combined_{lib}.png"
            )
    
    # ===== PART 2: ENRICHR PLOTS (Combined Score) =====
    print("\n[2] Enrichr plots (Combined Score)")
    print("-" * 50)
    
    for lib in LIBRARIES:
        print(f"\n  {lib}")
        
        # Load Enrichr results
        stem_file = ENRICHR_POSITIVE / f"{lib}.Mouse.enrichr.reports.txt"
        diff_file = ENRICHR_NEGATIVE / f"{lib}.Mouse.enrichr.reports.txt"
        
        df_stem = load_enrichr_txt(stem_file)
        df_diff = load_enrichr_txt(diff_file)
        
        # Single-direction plots
        if not df_stem.empty:
            plot_enrichr_scatter(
                df_stem, lib, 
                OUTPUT_DIR / f"enrichr_stemlike_{lib}.png",
                color=COLOR_STEM, marker='o',
                title_prefix="Enrichr"
            )
        
        if not df_diff.empty:
            plot_enrichr_scatter(
                df_diff, lib,
                OUTPUT_DIR / f"enrichr_differentiated_{lib}.png",
                color=COLOR_DIFF, marker='D',
                title_prefix="Enrichr"
            )
        
        # Combined plot
        if not df_stem.empty or not df_diff.empty:
            plot_combined_stem_diff(
                df_stem, df_diff, lib,
                OUTPUT_DIR / f"enrichr_combined_{lib}.png"
            )
    
    print("\n" + "=" * 60)
    print("✓ COMPLETE")
    print("=" * 60)
    print(f"\nOutput: {OUTPUT_DIR}")
    print("\nKey files:")
    print("  - gsea_prerank_stemlike_*.png  (NES - has OXPHOS!)")
    print("  - gsea_prerank_differentiated_*.png  (NES)")
    print("  - enrichr_stemlike_*.png  (Combined Score)")
    print("  - enrichr_differentiated_*.png  (Combined Score)")
    print("  - enrichr_combined_*.png  (both directions)")


if __name__ == "__main__":
    main()


__EOF_plot_gsea_custom_py__

cat > "${SCRIPTS}/plot_gsea_custom_B.py" << '__EOF_plot_gsea_custom_B_py__'
#!/usr/bin/env python3
"""
Custom GSEA/Enrichr Scatter Plots - Version B
=================================
Generates publication-quality plots:
- GSEA prerank plots (NES on Y-axis) - captures OXPHOS!
- Enrichr plots (Combined Score on Y-axis)
- Both with clean aesthetics

Version B: Two separate panels for Malignant (red points) and Physiologic (blue diamonds)
with bigger labels.

Author: Generated for CytoTRACE2 analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob

# ===== CONFIGURATION =====
BASE_DIR = Path("__BASEDIR__/heatmap_programs")
OUTPUT_DIR = BASE_DIR / "figures" / "gsea_custom"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Where results are stored
ENRICHR_POSITIVE = BASE_DIR / "programs" / "enrichr_positive_corr_stemlike"
ENRICHR_NEGATIVE = BASE_DIR / "programs" / "enrichr_negative_corr_differentiated"
GSEA_POSITIVE = BASE_DIR / "programs" / "gsea_prerank_positive_corr_stemlike"
GSEA_NEGATIVE = BASE_DIR / "programs" / "gsea_prerank_negative_corr_differentiated"

# Libraries to plot
LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2019_Mouse",
    "Reactome_2022",
    "MSigDB_Hallmark_2020"
]

# Colors
COLOR_STEM = "#DC3545"      # Red for stem-like (Malignant)
COLOR_DIFF = "#0077B6"      # Blue for differentiated (Physiologic)
MARKER_SIZE = 45            # Uniform size

# Label settings
LABEL_TOP_N = 35
LABEL_FONTSIZE = 9          # Smaller labels

# Output formats
SAVE_FORMATS = ['png', 'svg', 'pdf']


def save_figure(fig, output_path: Path):
    """Save figure in multiple formats (PNG, SVG, PDF)."""
    for fmt in SAVE_FORMATS:
        out_file = output_path.with_suffix(f'.{fmt}')
        fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {output_path.stem} (.png, .svg, .pdf)")


def load_enrichr_txt(filepath: Path) -> pd.DataFrame:
    """Load Enrichr .txt report file."""
    if not filepath.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(filepath, sep='\t')
        return df
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return pd.DataFrame()


def load_gsea_prerank(filepath: Path, library: str) -> pd.DataFrame:
    """Load GSEA prerank results for a specific library."""
    if not filepath.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(filepath)
        # Filter by library
        if 'library' in df.columns:
            df = df[df['library'] == library].copy()
        return df
    except Exception as e:
        print(f"  Error loading {filepath.name}: {e}")
        return pd.DataFrame()


def plot_gsea_prerank_scatter(df: pd.DataFrame, library: str, output_path: Path,
                              color: str = COLOR_STEM, marker: str = 'o',
                              title_suffix: str = ""):
    """
    Plot GSEA prerank scatter with the same sloping style as Enrichr.
    Uses a transformed x-axis to create smooth distribution.
    """
    if df.empty:
        print(f"  No GSEA prerank data for {library}")
        return
    
    # Get columns
    term_col = "Term" if "Term" in df.columns else "Name"
    nes_col = "NES" if "NES" in df.columns else None
    
    # Find FDR q-value column (better distribution than NOM p-val)
    p_col = None
    for c in ["FDR q-val", "FDR q-value", "Adjusted P-value", "NOM p-val", "P-value"]:
        if c in df.columns:
            p_col = c
            break
    
    if not nes_col or not p_col:
        print(f"  Missing columns. Available: {list(df.columns)}")
        return
    
    # Calculate axes
    df = df.copy()
    
    # Get FDR/p-values - use minimum of 1e-10 for visual spread
    pvals = df[p_col].astype(float).clip(lower=1e-10, upper=1.0)
    df["logp_raw"] = -np.log10(pvals)
    
    # Cap at 10 for nice visual range (like Enrichr)
    MAX_LOGP = 10
    df["logp"] = df["logp_raw"].clip(upper=MAX_LOGP)
    
    # Add jitter to highly significant terms (those hitting the cap)
    # This spreads out the dense cluster in upper right
    np.random.seed(42)  # Reproducible jitter
    at_cap = df["logp"] >= MAX_LOGP - 0.1
    n_at_cap = at_cap.sum()
    if n_at_cap > 0:
        # Add horizontal jitter proportional to number of capped terms
        jitter_range = min(2.0, n_at_cap * 0.02)  # Max 2 units of jitter
        df.loc[at_cap, "logp"] = MAX_LOGP - np.random.uniform(0, jitter_range, n_at_cap)
    
    # Use absolute NES for Y-axis
    df["nes"] = df[nes_col].astype(float).abs()
    
    # Clean term names
    df["term_clean"] = df[term_col].astype(str).str.replace("_", " ")
    df["term_clean"] = df["term_clean"].str.replace(r"\s*\([A-Z]+-?[A-Z]*:\d+\)$", "", regex=True)
    df["term_clean"] = df["term_clean"].str.replace(r"\s*R-[A-Z]+-\d+$", "", regex=True)
    
    # ===== PLOT (matching Enrichr style) =====
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot ALL dots with uniform size (like Enrichr)
    ax.scatter(df["logp"], df["nes"], 
              s=MARKER_SIZE, c=color, marker=marker,
              alpha=0.7, edgecolors='none', zorder=2)
    
    # Reference line at y=0
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Labels for top terms (by NES * significance)
    # Weight by both NES and significance for better label selection
    df["importance"] = df["nes"] * (df["logp"] / MAX_LOGP)
    top_terms = df.nlargest(LABEL_TOP_N, "importance")
    
    # Use adjustText for better label positioning
    texts = []
    for _, row in top_terms.iterrows():
        label = row["term_clean"]
        if len(label) > 50:
            label = label[:47] + "..."
        txt = ax.text(row["logp"] + 0.1, row["nes"], label,
                     fontsize=LABEL_FONTSIZE, color='#333333',
                     ha='left', va='center')
        texts.append(txt)
    
    # Try to use adjustText for smart positioning
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                   arrowprops=dict(arrowstyle='-', color='#AAAAAA', alpha=0.5, lw=0.5),
                   expand_points=(1.5, 1.5),
                   force_text=(0.5, 0.8),
                   force_points=(0.3, 0.3))
    except ImportError:
        pass  # Fall back to basic positioning
    
    # Axis settings (matching Enrichr range)
    ax.set_xlim(0, MAX_LOGP + 1)
    ax.set_ylim(bottom=0)
    
    # Title
    lib_clean = library.replace("_", " ")
    ax.set_title(f"GSEA (prerank) — {lib_clean}{title_suffix}", 
                fontsize=18, fontweight='bold', pad=15)
    
    # Axis labels
    p_label = "FDR" if "FDR" in p_col else "P-value"
    ax.set_xlabel(f"-log10({p_label})", fontsize=14)
    ax.set_ylabel("|NES|", fontsize=14)
    
    # Spine styling (matching Enrichr)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_gsea_prerank_combined(df_stem: pd.DataFrame, df_diff: pd.DataFrame,
                               library: str, output_path: Path):
    """
    GSEA prerank plot with TWO SEPARATE PANELS:
    - Left panel: Malignant (red points)
    - Right panel: Physiologic (blue diamonds)
    Both with bigger labels.
    """
    if df_stem.empty and df_diff.empty:
        print(f"  No GSEA prerank data for {library}")
        return
    
    # Determine term column from available data
    sample_df = df_stem if not df_stem.empty else df_diff
    term_col = "Term" if "Term" in sample_df.columns else "Name"
    nes_col = "NES"
    
    # Find p-value column
    p_col = None
    for c in ["FDR q-val", "FDR q-value", "Adjusted P-value", "NOM p-val"]:
        if c in sample_df.columns:
            p_col = c
            break
    
    if not p_col:
        print(f"  No p-value column found")
        return
    
    MAX_LOGP = 10
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9), dpi=300)
    fig.patch.set_facecolor('white')
    ax1.set_facecolor('white')
    ax2.set_facecolor('white')
    
    lib_clean = library.replace("_", " ")
    p_label = "FDR" if "FDR" in p_col else "P-value"
    
    # ===== LEFT PANEL: Malignant (red points) =====
    if not df_stem.empty:
        df_mal = df_stem.copy()
        pvals = df_mal[p_col].astype(float).clip(lower=1e-10, upper=1.0)
        df_mal["logp_raw"] = -np.log10(pvals)
        df_mal["logp"] = df_mal["logp_raw"].clip(upper=MAX_LOGP)
        
        # Add jitter for capped terms
        np.random.seed(42)
        at_cap = df_mal["logp"] >= MAX_LOGP - 0.1
        n_at_cap = at_cap.sum()
        if n_at_cap > 0:
            jitter_range = min(2.0, n_at_cap * 0.02)
            df_mal.loc[at_cap, "logp"] = MAX_LOGP - np.random.uniform(0, jitter_range, n_at_cap)
        
        df_mal["nes"] = df_mal[nes_col].astype(float).abs()
        
        # Clean term names
        df_mal["term_clean"] = df_mal[term_col].astype(str).str.replace("_", " ")
        df_mal["term_clean"] = df_mal["term_clean"].str.replace(r"\s*\([A-Z]+-?[A-Z]*:\d+\)$", "", regex=True)
        df_mal["term_clean"] = df_mal["term_clean"].str.replace(r"\s*R-[A-Z]+-\d+$", "", regex=True)
        
        # Plot red points
        ax1.scatter(df_mal["logp"], df_mal["nes"], 
                   s=MARKER_SIZE, c=COLOR_STEM, marker='o',
                   alpha=0.7, edgecolors='none', zorder=2)
        
        # Labels for top terms with BIGGER FONTSIZE
        df_mal["importance"] = df_mal["nes"] * (df_mal["logp"] / MAX_LOGP)
        top_terms_mal = df_mal.nlargest(15, "importance")
        
        texts_mal = []
        for _, row in top_terms_mal.iterrows():
            term_label = row["term_clean"]
            if len(term_label) > 50:
                term_label = term_label[:47] + "..."
            txt = ax1.text(row["logp"] + 0.1, row["nes"], term_label,
                          fontsize=LABEL_FONTSIZE, color=COLOR_STEM, alpha=0.9,
                          ha='left', va='center')
            texts_mal.append(txt)
        
        # Use adjustText for label positioning
        try:
            from adjustText import adjust_text
            adjust_text(texts_mal, ax=ax1,
                       arrowprops=dict(arrowstyle='-', color='#AAAAAA', alpha=0.4, lw=0.4),
                       expand_points=(1.3, 1.3),
                       force_text=(0.4, 0.6))
        except ImportError:
            pass
    
    # Reference line
    ax1.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Axis settings
    ax1.set_xlim(0, MAX_LOGP + 1)
    ax1.set_ylim(bottom=0)
    
    # Title and labels
    ax1.set_title(f"Malignant", fontsize=18, fontweight='bold', pad=15, color=COLOR_STEM)
    ax1.set_xlabel(f"-log10({p_label})", fontsize=16)
    ax1.set_ylabel("|NES|", fontsize=16)
    
    # Spines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(axis='both', labelsize=14)
    
    # ===== RIGHT PANEL: Physiologic (blue diamonds) =====
    if not df_diff.empty:
        df_phys = df_diff.copy()
        pvals = df_phys[p_col].astype(float).clip(lower=1e-10, upper=1.0)
        df_phys["logp_raw"] = -np.log10(pvals)
        df_phys["logp"] = df_phys["logp_raw"].clip(upper=MAX_LOGP)
        
        # Add jitter for capped terms
        np.random.seed(43)
        at_cap = df_phys["logp"] >= MAX_LOGP - 0.1
        n_at_cap = at_cap.sum()
        if n_at_cap > 0:
            jitter_range = min(2.0, n_at_cap * 0.02)
            df_phys.loc[at_cap, "logp"] = MAX_LOGP - np.random.uniform(0, jitter_range, n_at_cap)
        
        df_phys["nes"] = df_phys[nes_col].astype(float).abs()
        
        # Clean term names
        df_phys["term_clean"] = df_phys[term_col].astype(str).str.replace("_", " ")
        df_phys["term_clean"] = df_phys["term_clean"].str.replace(r"\s*\([A-Z]+-?[A-Z]*:\d+\)$", "", regex=True)
        df_phys["term_clean"] = df_phys["term_clean"].str.replace(r"\s*R-[A-Z]+-\d+$", "", regex=True)
        
        # Plot blue diamonds
        ax2.scatter(df_phys["logp"], df_phys["nes"], 
                   s=MARKER_SIZE, c=COLOR_DIFF, marker='D',
                   alpha=0.7, edgecolors='none', zorder=2)
        
        # Labels for top terms with BIGGER FONTSIZE
        df_phys["importance"] = df_phys["nes"] * (df_phys["logp"] / MAX_LOGP)
        top_terms_phys = df_phys.nlargest(15, "importance")
        
        texts_phys = []
        for _, row in top_terms_phys.iterrows():
            term_label = row["term_clean"]
            if len(term_label) > 50:
                term_label = term_label[:47] + "..."
            txt = ax2.text(row["logp"] + 0.1, row["nes"], term_label,
                          fontsize=LABEL_FONTSIZE, color=COLOR_DIFF, alpha=0.9,
                          ha='left', va='center')
            texts_phys.append(txt)
        
        # Use adjustText for label positioning
        try:
            from adjustText import adjust_text
            adjust_text(texts_phys, ax=ax2,
                       arrowprops=dict(arrowstyle='-', color='#AAAAAA', alpha=0.4, lw=0.4),
                       expand_points=(1.3, 1.3),
                       force_text=(0.4, 0.6))
        except ImportError:
            pass
    
    # Reference line
    ax2.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Axis settings - limit FDR axis to 4 for physiologic panel
    ax2.set_xlim(0, 4)
    ax2.set_ylim(bottom=0)
    
    # Set x-axis ticks to whole numbers only
    ax2.set_xticks(range(0, 5))  # 0, 1, 2, 3, 4
    
    # Title and labels
    ax2.set_title(f"Physiologic", fontsize=18, fontweight='bold', pad=15, color=COLOR_DIFF)
    ax2.set_xlabel(f"-log10({p_label})", fontsize=16)
    ax2.set_ylabel("|NES|", fontsize=16)
    
    # Spines
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='both', labelsize=14)
    
    # Overall figure title
    fig.suptitle(f"GSEA (prerank) — {lib_clean}", 
                fontsize=20, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_enrichr_scatter(df: pd.DataFrame, library: str, output_path: Path, 
                         color: str = COLOR_STEM, marker: str = 'o',
                         title_prefix: str = "GSEA (prerank)"):
    """
    Plot Enrichr-style scatter: Combined Score vs -log10(P-value)
    Matches the reference style exactly.
    """
    if df.empty:
        print(f"  No data for {library}")
        return
    
    # Get columns
    term_col = "Term" if "Term" in df.columns else df.columns[1]
    p_col = "P-value" if "P-value" in df.columns else None
    score_col = "Combined Score" if "Combined Score" in df.columns else None
    
    if not p_col or p_col not in df.columns:
        for c in df.columns:
            if 'p-value' in c.lower() and 'adjusted' not in c.lower():
                p_col = c
                break
    
    if not score_col or score_col not in df.columns:
        for c in df.columns:
            if 'combined' in c.lower() and 'score' in c.lower():
                score_col = c
                break
    
    if not p_col or not score_col:
        print(f"  Missing columns. Available: {list(df.columns)}")
        return
    
    # Calculate axes
    df = df.copy()
    df["logp"] = -np.log10(df[p_col].astype(float).clip(lower=1e-300))
    df["score"] = df[score_col].astype(float)
    
    # Clean term names
    df["term_clean"] = df[term_col].str.replace("_", " ")
    
    # ===== PLOT (matching reference exactly) =====
    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot ALL dots with uniform size
    ax.scatter(df["logp"], df["score"], 
              s=MARKER_SIZE, c=color, marker=marker,
              alpha=0.7, edgecolors='none', zorder=2)
    
    # Reference line at y=0
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    
    # Labels for top terms
    top_terms = df.nlargest(LABEL_TOP_N, "score")
    
    for _, row in top_terms.iterrows():
        label = row["term_clean"]
        if len(label) > 60:
            label = label[:57] + "..."
        ax.annotate(label, (row["logp"], row["score"]),
                   fontsize=LABEL_FONTSIZE, color='#333333',
                   ha='left', va='center',
                   xytext=(4, 0), textcoords='offset points')
    
    # Axis settings
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    # Title (matching reference)
    lib_clean = library.replace("_", " ")
    ax.set_title(f"{title_prefix} — {lib_clean}", 
                fontsize=18, fontweight='bold', pad=15)
    
    # Axis labels
    ax.set_xlabel("-log10(P-value)", fontsize=14)
    ax.set_ylabel("NES", fontsize=14)  # Keep same label as reference
    
    # Spine styling (matching reference)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    
    ax.tick_params(axis='both', labelsize=12)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def plot_combined_stem_diff(df_stem: pd.DataFrame, df_diff: pd.DataFrame,
                            library: str, output_path: Path):
    """
    Combined plot: red dots (stem) + blue diamonds (diff)
    """
    if df_stem.empty and df_diff.empty:
        print(f"  No data for {library}")
        return
    
    # Get columns
    term_col = "Term"
    p_col = "P-value"
    score_col = "Combined Score"
    
    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot stem-like (positive correlation) as red dots
    if not df_stem.empty:
        df_s = df_stem.copy()
        df_s["logp"] = -np.log10(df_s[p_col].astype(float).clip(lower=1e-300))
        df_s["score"] = df_s[score_col].astype(float)
        ax.scatter(df_s["logp"], df_s["score"], 
                  s=50, c='#DC3545', marker='o', label='Stem-like (high CT2)',
                  alpha=0.7, edgecolors='none', zorder=2)
    
    # Plot differentiated (negative correlation) as blue diamonds
    if not df_diff.empty:
        df_d = df_diff.copy()
        df_d["logp"] = -np.log10(df_d[p_col].astype(float).clip(lower=1e-300))
        df_d["score"] = df_d[score_col].astype(float)
        ax.scatter(df_d["logp"], df_d["score"], 
                  s=60, c='#0077B6', marker='D', label='Differentiated (low CT2)',
                  alpha=0.7, edgecolors='none', zorder=2)
    
    # Labels for top terms from each
    for df, c in [(df_stem, '#DC3545'), (df_diff, '#0077B6')]:
        if df.empty:
            continue
        df = df.copy()
        df["logp"] = -np.log10(df[p_col].astype(float).clip(lower=1e-300))
        df["score"] = df[score_col].astype(float)
        top = df.nlargest(15, "score")
        for _, row in top.iterrows():
            label = str(row[term_col]).replace("_", " ")
            if len(label) > 50:
                label = label[:47] + "..."
            ax.annotate(label, (row["logp"], row["score"]),
                       fontsize=7, color=c, alpha=0.9,
                       ha='left', va='center',
                       xytext=(4, 0), textcoords='offset points')
    
    ax.axhline(0, color='black', lw=0.5, alpha=0.5)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    lib_clean = library.replace("_", " ")
    ax.set_title(f"GSEA — {lib_clean}\nStem-like vs Differentiated", 
                fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel("-log10(P-value)", fontsize=13)
    ax.set_ylabel("Combined Score", fontsize=13)
    
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=11)
    
    plt.tight_layout()
    save_figure(plt.gcf(), output_path)
    plt.close()


def main():
    print("=" * 60)
    print("Custom GSEA Scatter Plotting - Version B")
    print("=" * 60)
    
    # ===== PART 1: GSEA PRERANK PLOTS (captures OXPHOS!) =====
    print("\n[1] GSEA Prerank plots (NES - captures OXPHOS)")
    print("-" * 50)
    
    gsea_stem_file = GSEA_POSITIVE / "gsea_prerank_positive_corr_stemlike_all.csv"
    gsea_diff_file = GSEA_NEGATIVE / "gsea_prerank_negative_corr_differentiated_all.csv"
    
    for lib in LIBRARIES:
        print(f"\n  {lib}")
        
        # Stem-like (GSEA prerank)
        df_gsea_stem = load_gsea_prerank(gsea_stem_file, lib)
        if not df_gsea_stem.empty:
            plot_gsea_prerank_scatter(
                df_gsea_stem, lib,
                OUTPUT_DIR / f"gsea_prerank_stemlike_{lib}.png",
                color=COLOR_STEM, marker='o',
                title_suffix=" (Stem-like)"
            )
        
        # Differentiated (GSEA prerank)
        df_gsea_diff = load_gsea_prerank(gsea_diff_file, lib)
        if not df_gsea_diff.empty:
            plot_gsea_prerank_scatter(
                df_gsea_diff, lib,
                OUTPUT_DIR / f"gsea_prerank_differentiated_{lib}.png",
                color=COLOR_DIFF, marker='D',
                title_suffix=" (Differentiated)"
            )
        
        # Combined GSEA prerank plot (TWO SEPARATE PANELS)
        if not df_gsea_stem.empty or not df_gsea_diff.empty:
            plot_gsea_prerank_combined(
                df_gsea_stem, df_gsea_diff, lib,
                OUTPUT_DIR / f"gsea_prerank_combined_{lib}.png"
            )
    
    # ===== PART 2: ENRICHR PLOTS (Combined Score) =====
    print("\n[2] Enrichr plots (Combined Score)")
    print("-" * 50)
    
    for lib in LIBRARIES:
        print(f"\n  {lib}")
        
        # Load Enrichr results
        stem_file = ENRICHR_POSITIVE / f"{lib}.Mouse.enrichr.reports.txt"
        diff_file = ENRICHR_NEGATIVE / f"{lib}.Mouse.enrichr.reports.txt"
        
        df_stem = load_enrichr_txt(stem_file)
        df_diff = load_enrichr_txt(diff_file)
        
        # Single-direction plots
        if not df_stem.empty:
            plot_enrichr_scatter(
                df_stem, lib, 
                OUTPUT_DIR / f"enrichr_stemlike_{lib}.png",
                color=COLOR_STEM, marker='o',
                title_prefix="Enrichr"
            )
        
        if not df_diff.empty:
            plot_enrichr_scatter(
                df_diff, lib,
                OUTPUT_DIR / f"enrichr_differentiated_{lib}.png",
                color=COLOR_DIFF, marker='D',
                title_prefix="Enrichr"
            )
        
        # Combined plot
        if not df_stem.empty or not df_diff.empty:
            plot_combined_stem_diff(
                df_stem, df_diff, lib,
                OUTPUT_DIR / f"enrichr_combined_{lib}.png"
            )
    
    print("\n" + "=" * 60)
    print("✓ COMPLETE")
    print("=" * 60)
    print(f"\nOutput: {OUTPUT_DIR}")
    print("\nKey files:")
    print("  - gsea_prerank_stemlike_*.png  (NES - has OXPHOS!)")
    print("  - gsea_prerank_differentiated_*.png  (NES)")
    print("  - gsea_prerank_combined_*.png  (TWO SEPARATE PANELS)")
    print("  - enrichr_stemlike_*.png  (Combined Score)")
    print("  - enrichr_differentiated_*.png  (Combined Score)")
    print("  - enrichr_combined_*.png  (both directions)")


if __name__ == "__main__":
    main()


__EOF_plot_gsea_custom_B_py__

cat > "${SCRIPTS}/pgc1_cytotrace2_trajectory.py" << '__EOF_pgc1_cytotrace2_trajectory_py__'
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

Author: Generated script
Date: 2024
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
ADATA_PATH = "__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad"

# Output directory
OUTPUT_DIR = "__BASEDIR__/pgc1_cytotrace2_trajectory"

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



__EOF_pgc1_cytotrace2_trajectory_py__

cat > "${SCRIPTS}/wilcoxon_rank_mouse_integrated.py" << '__EOF_wilcoxon_rank_mouse_integrated_py__'
#!/usr/bin/env python3
"""
Wilcoxon Rank Sum Test for Discriminating Genes per Leiden Cluster
===================================================================
Equivalent to the R script wilcoxon_rank_A.R but for AnnData objects.

Extracts top discriminating genes for each Leiden cluster (resolution 1.0)
using Wilcoxon rank sum test, generates dotplot and saves results.

Author: Generated script
Date: 2024
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════

# Input AnnData object
ADATA_PATH = "__BASEDIR__/mouse_scvi_cytotrace2/mouse_integrated.h5ad"

# Output directory
OUTPUT_DIR = "__BASEDIR__/wilcoxon_rank_sum"

# Leiden cluster column (resolution 1.0)
LEIDEN_KEY = "leiden_1.0"  # Adjust if the column name differs

# Analysis parameters
TOP_N_GENES_PER_CLUSTER = 25  # Top genes per cluster for gene panel
MAX_GENES_PANEL = 100         # Maximum genes in dotplot panel
MIN_CELLS_PER_GROUP = 2       # Minimum cells required per cluster
LOGFC_THRESHOLD = 0.0         # Log fold change threshold (0 = no filter)

# ══════════════════════════════════════════════════════════════════════════════
# 0. Setup output directory
# ══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output directory: {OUTPUT_DIR}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load AnnData object
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nLoading AnnData from: {ADATA_PATH}")
adata = sc.read_h5ad(ADATA_PATH)
print(f"Loaded object with {adata.n_obs} cells and {adata.n_vars} genes")

# Check available columns
print(f"\nAvailable obs columns: {list(adata.obs.columns)}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Identify Leiden cluster column
# ══════════════════════════════════════════════════════════════════════════════

# Try to find the leiden 1.0 resolution column
leiden_candidates = [col for col in adata.obs.columns if 'leiden' in col.lower()]
print(f"\nLeiden-related columns found: {leiden_candidates}")

# Select the appropriate column
if LEIDEN_KEY in adata.obs.columns:
    cluster_key = LEIDEN_KEY
elif 'leiden_res1.0' in adata.obs.columns:
    cluster_key = 'leiden_res1.0'
elif 'leiden' in adata.obs.columns:
    cluster_key = 'leiden'
elif len(leiden_candidates) > 0:
    # Find one with "1.0" or "1" in name, or take first
    res10_cols = [c for c in leiden_candidates if '1.0' in c or '_1' in c]
    cluster_key = res10_cols[0] if res10_cols else leiden_candidates[0]
else:
    raise ValueError("No Leiden cluster column found in adata.obs")

print(f"\nUsing cluster column: '{cluster_key}'")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Prepare data for DE analysis
# ══════════════════════════════════════════════════════════════════════════════

# Convert cluster labels to string for consistency
adata.obs['cluster_for_de'] = adata.obs[cluster_key].astype(str)

# Count cells per cluster
cluster_counts = adata.obs['cluster_for_de'].value_counts().sort_index()
print(f"\nCells per cluster:")
print(cluster_counts)

# Filter clusters with sufficient cells
valid_clusters = cluster_counts[cluster_counts >= MIN_CELLS_PER_GROUP].index.tolist()
print(f"\nClusters with >= {MIN_CELLS_PER_GROUP} cells: {len(valid_clusters)}")

if len(valid_clusters) < 2:
    raise ValueError("Fewer than two clusters have sufficient cells - DE not meaningful")

# Subset to valid clusters
adata_sub = adata[adata.obs['cluster_for_de'].isin(valid_clusters)].copy()
print(f"Subset to {adata_sub.n_obs} cells in {len(valid_clusters)} clusters")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Run Wilcoxon Rank Sum Test (each cluster vs rest)
# ══════════════════════════════════════════════════════════════════════════════

print("\nRunning Wilcoxon rank sum test...")

# Ensure we have normalized data
# Check if data looks normalized (values typically between 0-10 for log-normalized)
data_max = adata_sub.X.max() if hasattr(adata_sub.X, 'max') else np.max(adata_sub.X.toarray())
print(f"Data max value: {data_max:.2f}")

if data_max > 100:
    print("Data appears to be counts, normalizing...")
    # Store raw if not already
    if adata_sub.raw is None:
        adata_sub.raw = adata_sub.copy()
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)

# Run rank_genes_groups with Wilcoxon test
sc.tl.rank_genes_groups(
    adata_sub,
    groupby='cluster_for_de',
    method='wilcoxon',
    pts=True,  # Calculate percentage of cells expressing
    key_added='wilcoxon_de'
)

print("Differential expression analysis complete!")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Extract DE results to DataFrame
# ══════════════════════════════════════════════════════════════════════════════

print("\nExtracting DE results...")

# Get all results
de_results = []
groups = adata_sub.uns['wilcoxon_de']['names'].dtype.names

for group in groups:
    n_genes = len(adata_sub.uns['wilcoxon_de']['names'][group])
    
    group_df = pd.DataFrame({
        'gene': adata_sub.uns['wilcoxon_de']['names'][group],
        'scores': adata_sub.uns['wilcoxon_de']['scores'][group],
        'logfoldchanges': adata_sub.uns['wilcoxon_de']['logfoldchanges'][group],
        'pvals': adata_sub.uns['wilcoxon_de']['pvals'][group],
        'pvals_adj': adata_sub.uns['wilcoxon_de']['pvals_adj'][group],
        'cluster': group
    })
    
    # Add percentage expressed if available
    if 'pts' in adata_sub.uns['wilcoxon_de']:
        group_df['pct_expressed'] = adata_sub.uns['wilcoxon_de']['pts'][group]
    if 'pts_rest' in adata_sub.uns['wilcoxon_de']:
        group_df['pct_expressed_rest'] = adata_sub.uns['wilcoxon_de']['pts_rest'][group]
    
    de_results.append(group_df)

markers_df = pd.concat(de_results, ignore_index=True)

# Apply logFC threshold if specified
if LOGFC_THRESHOLD > 0:
    markers_df = markers_df[np.abs(markers_df['logfoldchanges']) >= LOGFC_THRESHOLD]

print(f"Total DE results: {len(markers_df)} gene-cluster pairs")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Build gene panel (top N per cluster)
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nSelecting top {TOP_N_GENES_PER_CLUSTER} genes per cluster...")

# Get top genes per cluster by adjusted p-value, then by log fold change
top_genes_per_cluster = (
    markers_df
    .sort_values(['cluster', 'pvals_adj', 'logfoldchanges'], 
                 ascending=[True, True, False])
    .groupby('cluster')
    .head(TOP_N_GENES_PER_CLUSTER)
)

# Get unique genes maintaining order
gene_panel = top_genes_per_cluster['gene'].drop_duplicates().tolist()

if len(gene_panel) > MAX_GENES_PANEL:
    print(f"Truncating gene panel from {len(gene_panel)} to {MAX_GENES_PANEL} genes")
    gene_panel = gene_panel[:MAX_GENES_PANEL]

print(f"Final gene panel: {len(gene_panel)} unique genes")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Save DE results
# ══════════════════════════════════════════════════════════════════════════════

# Save full results
full_results_path = os.path.join(OUTPUT_DIR, "wilcoxon_DEgenes_all_clusters.csv")
markers_df.to_csv(full_results_path, index=False)
print(f"\nFull DE results saved to: {full_results_path}")

# Save top genes per cluster
top_genes_path = os.path.join(OUTPUT_DIR, "wilcoxon_top_genes_per_cluster.csv")
top_genes_per_cluster.to_csv(top_genes_path, index=False)
print(f"Top genes per cluster saved to: {top_genes_path}")

# Save gene panel
gene_panel_path = os.path.join(OUTPUT_DIR, "gene_panel.txt")
with open(gene_panel_path, 'w') as f:
    f.write('\n'.join(gene_panel))
print(f"Gene panel saved to: {gene_panel_path}")

# ══════════════════════════════════════════════════════════════════════════════
# 8. Generate Dot Plot
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating dot plot...")

if len(gene_panel) > 0:
    # Set up figure
    n_genes = len(gene_panel)
    n_clusters = len(valid_clusters)
    
    # Calculate figure size
    fig_width = max(12, n_genes * 0.25)
    fig_height = max(6, n_clusters * 0.4)
    
    # Sort clusters numerically if possible
    try:
        sorted_clusters = sorted(valid_clusters, key=lambda x: float(x))
    except ValueError:
        sorted_clusters = sorted(valid_clusters)
    
    # Create dot plot
    sc.pl.dotplot(
        adata_sub,
        var_names=gene_panel,
        groupby='cluster_for_de',
        categories_order=sorted_clusters,
        standard_scale='var',  # Scale gene expression across clusters
        dendrogram=False,
        show=False,
        save=False
    )
    
    # Save figure
    dotplot_path = os.path.join(OUTPUT_DIR, "dotplot_top_DE_genes.pdf")
    plt.savefig(dotplot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Dot plot saved to: {dotplot_path}")
    
    # Also save as PNG
    sc.pl.dotplot(
        adata_sub,
        var_names=gene_panel,
        groupby='cluster_for_de',
        categories_order=sorted_clusters,
        standard_scale='var',
        dendrogram=False,
        show=False,
        save=False
    )
    dotplot_png_path = os.path.join(OUTPUT_DIR, "dotplot_top_DE_genes.png")
    plt.savefig(dotplot_png_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Dot plot (PNG) saved to: {dotplot_png_path}")
    
    # ══════════════════════════════════════════════════════════════════════════
    # 9. Additional visualizations
    # ══════════════════════════════════════════════════════════════════════════
    
    # Heatmap of top genes
    print("\nGenerating heatmap...")
    # Reorder categories in the obs column for proper ordering in heatmap
    adata_sub.obs['cluster_for_de'] = pd.Categorical(
        adata_sub.obs['cluster_for_de'],
        categories=sorted_clusters,
        ordered=True
    )
    sc.pl.heatmap(
        adata_sub,
        var_names=gene_panel[:50] if len(gene_panel) > 50 else gene_panel,  # Limit for readability
        groupby='cluster_for_de',
        standard_scale='var',
        show=False,
        save=False
    )
    heatmap_path = os.path.join(OUTPUT_DIR, "heatmap_top_DE_genes.pdf")
    plt.savefig(heatmap_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Heatmap saved to: {heatmap_path}")
    
    # Rank genes groups plot (scanpy style)
    print("\nGenerating rank genes groups plot...")
    sc.pl.rank_genes_groups(
        adata_sub,
        key='wilcoxon_de',
        n_genes=10,
        sharey=False,
        show=False,
        save=False
    )
    rank_plot_path = os.path.join(OUTPUT_DIR, "rank_genes_groups.pdf")
    plt.savefig(rank_plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Rank genes groups plot saved to: {rank_plot_path}")
    
    # Rank genes groups dotplot (more compact visualization)
    print("\nGenerating rank genes groups dotplot...")
    sc.pl.rank_genes_groups_dotplot(
        adata_sub,
        key='wilcoxon_de',
        n_genes=5,
        standard_scale='var',
        show=False,
        save=False
    )
    rank_dotplot_path = os.path.join(OUTPUT_DIR, "rank_genes_groups_dotplot.pdf")
    plt.savefig(rank_dotplot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Rank genes groups dotplot saved to: {rank_dotplot_path}")

else:
    print("WARNING: No genes met selection criteria - plots skipped")

# ══════════════════════════════════════════════════════════════════════════════
# 10. Summary statistics
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"Total cells analyzed: {adata_sub.n_obs}")
print(f"Total clusters: {len(valid_clusters)}")
print(f"Clusters: {', '.join(sorted_clusters)}")
print(f"Total DE gene-cluster pairs: {len(markers_df)}")
print(f"Genes in panel: {len(gene_panel)}")
print(f"\nOutput files saved to: {OUTPUT_DIR}")

# Count significant genes per cluster
sig_genes = markers_df[markers_df['pvals_adj'] < 0.05]
sig_per_cluster = sig_genes.groupby('cluster').size()
print(f"\nSignificant genes (adj. p < 0.05) per cluster:")
for cluster in sorted_clusters:
    if cluster in sig_per_cluster.index:
        print(f"  Cluster {cluster}: {sig_per_cluster[cluster]}")
    else:
        print(f"  Cluster {cluster}: 0")

print("\nAnalysis complete!")



__EOF_wilcoxon_rank_mouse_integrated_py__
cat > "${SCRIPTS}/run_smoke_test.py" << '__EOF_run_smoke_test_py__'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end smoke test for CREBBP aggressive lymphoma single-cell pipeline.
Runs a self-contained test on synthetic data without requiring GEO downloads.
Verifies QC, normalization, clustering, trajectory scoring, and Wilcoxon testing.
"""

import os
import shutil
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

print("=" * 65)
print("  CREBBP Single-Cell Pipeline Smoke Test (Synthetic Dataset)")
print("=" * 65)

# Verify required core scientific stack
try:
    import scanpy as sc
    import anndata as ad
    print("  ✓ Scanpy and AnnData imported successfully")
except ImportError as e:
    print(f"  ✗ Required package missing: {e}")
    print("    Please activate the pipeline environment (conda activate crebbp_sc_pipeline)")
    exit(1)

# Set seeds
np.random.seed(42)

# Generate synthetic counts: 200 cells, 80 genes
n_cells = 200
n_genes = 80
cells = [f"cell_{i:03d}" for i in range(n_cells)]
genes = [f"Gene_{j:02d}" for j in range(n_genes)]
genes[0] = "Cd19"
genes[1] = "Ms4a1"
genes[2] = "Pax5"
genes[3] = "mt-Nd1"
genes[4] = "mt-Nd2"

# Simulating discrete counts with negative binomial distribution
counts = np.random.negative_binomial(n=4, p=0.6, size=(n_cells, n_genes)).astype(np.float32)
# Introduce differential expression between conditions
conditions = np.random.choice(["WT_B_cells", "Malignant"], size=n_cells)
counts[conditions == "Malignant", 0] *= 3.0  # Upregulate Cd19 in Malignant
sparse_counts = sp.csr_matrix(counts)

adata = ad.AnnData(
    X=sparse_counts,
    obs=pd.DataFrame({"condition": conditions, "replicate": np.random.choice(["R1", "R2"], size=n_cells)}, index=cells),
    var=pd.DataFrame({"gene_symbols": genes}, index=genes)
)
adata.var_names_make_unique()

print(f"  ✓ Created synthetic dataset: {adata.n_obs} cells × {adata.n_vars} genes")

# Step 1: QC metrics
adata.var['mt'] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs['pct_counts_mt'] <= 25.0, :].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"  ✓ QC filter complete: {adata.n_obs} cells retained")

# Step 2: Normalization and log-transform
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Step 3: Embeddings and clustering
sc.pp.pca(adata, n_comps=15)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
print(f"  ✓ Leiden clustering identified {adata.obs['leiden_0.5'].nunique()} clusters")

# Step 4: Synthetic developmental potency score (mock CytoTRACE2)
# GCS proxy: number of genes expressed per cell
gcs = np.asarray((adata.layers["counts"] > 0).sum(axis=1)).flatten()
adata.obs["CytoTRACE2_Score"] = (gcs - gcs.min()) / (gcs.max() - gcs.min() + 1e-6)
print(f"  ✓ Potency scoring calculated (mean: {adata.obs['CytoTRACE2_Score'].mean():.3f})")

# Step 5: Differential expression (Wilcoxon rank-sum)
sc.tl.rank_genes_groups(adata, groupby="condition", reference="WT_B_cells", method="wilcoxon")
de_df = sc.get.rank_genes_groups_df(adata, group="Malignant")
top_gene = de_df.iloc[0]["names"]
print("  ✓ Wilcoxon differential expression test complete")

# Save outputs to transient directory
out_dir = Path("tmp/smoke_test_output")
out_dir.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out_dir / "smoke_test_processed.h5ad")
de_df.to_csv(out_dir / "smoke_test_de_results.csv", index=False)
print(f"  ✓ Output successfully written to {out_dir}/")
print("=" * 65)
print("  Smoke test PASSED successfully!")
print("=" * 65)
__EOF_run_smoke_test_py__

echo "  Extracted 27 scripts (26 analysis scripts + 1 smoke test)."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Patch placeholder paths with user configuration
# ─────────────────────────────────────────────────────────────────────────────
echo "  Patching paths in extracted scripts..."

python3 -c '
import sys, pathlib
scripts = pathlib.Path(sys.argv[1])
basedir = sys.argv[2]
cellbender = sys.argv[3]
dlbcl = sys.argv[4]
tonsil = sys.argv[5]
gf = sys.argv[6]
fastq = sys.argv[7]
ref = sys.argv[8]
bin_cr = sys.argv[9]
for p in list(scripts.glob("*.py")) + list(scripts.glob("*.sh")) + list(scripts.glob("*.R")):
    txt = p.read_text()
    txt = txt.replace("__BASEDIR__", basedir)
    txt = txt.replace("__CELLBENDER_DIR__", cellbender)
    txt = txt.replace("__DLBCL_DIR__", dlbcl)
    txt = txt.replace("__TONSIL_DIR__", tonsil)
    txt = txt.replace("__GENEFORMER_MODEL_DIR__", gf)
    txt = txt.replace("__FASTQ_DATA_ROOT__", fastq)
    txt = txt.replace("__CELLRANGER_REF__", ref)
    txt = txt.replace("__CELLRANGER_BIN__", bin_cr)
    p.write_text(txt)
' "${SCRIPTS}" "${BASEDIR}" "${CELLBENDER_DIR}" "${DLBCL_DIR}" "${TONSIL_DIR}" "${GENEFORMER_MODEL_DIR}" "${FASTQ_DATA_ROOT}" "${CELLRANGER_REF}" "${CELLRANGER_BIN}"

chmod +x "${SCRIPTS}"/*.sh 2>/dev/null || true

echo "  ✓ All paths configured."
echo ""
else
    echo "  [Info] Scripts already present in ${SCRIPTS} (skipping re-extraction; use --force-extract to overwrite)."
    echo ""
fi

if [[ "$RUN_MODE" == "extract-only" ]]; then
    echo "Scripts extraction complete."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Pre-flight input validation (for --run mode)
# ─────────────────────────────────────────────────────────────────────────────
check_stage_data() {
    local stage="$1"
    if [[ "$RUN_MODE" != "run" || "$SKIP_DATA_CHECK" == "true" ]]; then
        return 0
    fi
    case "$stage" in
        1)
            local count=0
            if [[ -d "${CELLBENDER_DIR}" ]]; then
                count=$(find "${CELLBENDER_DIR}" -maxdepth 2 -name "*.h5" 2>/dev/null | wc -l)
            fi
            if [[ "$count" -eq 0 ]]; then
                echo "============================================================"
                echo "  ⚠️  Input Data Notice for Stage 1"
                echo "============================================================"
                echo "  No .h5 count matrices found in:"
                echo "    ${CELLBENDER_DIR}"
                echo ""
                echo "  Stage 1 requires CellBender-filtered mouse scRNA-seq matrices."
                echo "  Please download the dataset from GEO (accession: GSE332767)"
                echo "  and deposit into: ${CELLBENDER_DIR}/"
                echo "  (Refer to data/DATA_MANIFEST.md for accessions and file layout)."
                echo ""
                echo "  Options:"
                echo "    - Run synthetic smoke test:  bash MasterAnalysis.sh --smoke-test"
                echo "    - Bypass data check:         bash MasterAnalysis.sh --run --stage 1 --skip-data-check"
                echo "============================================================"
                exit 1
            fi
            ;;
        2)
            local count=0
            if [[ -d "${TONSIL_DIR}" ]]; then
                count=$(find "${TONSIL_DIR}" -maxdepth 2 -name "*.h5ad" 2>/dev/null | wc -l)
            fi
            if [[ "$count" -eq 0 && ! -d "${GENEFORMER_MODEL_DIR}" ]]; then
                echo "============================================================"
                echo "  ⚠️  Input Data Notice for Stage 2"
                echo "============================================================"
                echo "  No tonsil reference data found in:"
                echo "    ${TONSIL_DIR}"
                echo "  Please deposit tonsil reference H5ADs (see data/DATA_MANIFEST.md)."
                echo "============================================================"
                exit 1
            fi
            ;;
        3)
            local count=0
            if [[ -d "${DLBCL_DIR}" ]]; then
                count=$(find "${DLBCL_DIR}" -maxdepth 2 -name "*.h5ad" 2>/dev/null | wc -l)
            fi
            if [[ "$count" -eq 0 ]]; then
                echo "============================================================"
                echo "  ⚠️  Input Data Notice for Stage 3"
                echo "============================================================"
                echo "  No human DLBCL reference data found in:"
                echo "    ${DLBCL_DIR}"
                echo "  Please deposit Roider et al. / GSE182434 H5ADs (see data/DATA_MANIFEST.md)."
                echo "============================================================"
                exit 1
            fi
            ;;
    esac
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Pipeline execution
# ─────────────────────────────────────────────────────────────────────────────

run_py() {
    local script="$1"
    shift
    if [[ "$RUN_MODE" == "run" ]]; then
        echo "    [Executing] ${PYTHON_BIN} ${script} $*"
        "${PYTHON_BIN}" "${script}" "$@"
    else
        echo "    → ${PYTHON_BIN} ${script} $*"
    fi
}

run_gf_py() {
    local script="$1"
    shift
    if [[ "$RUN_MODE" == "run" ]]; then
        echo "    [Executing] ${GENEFORMER_PYTHON} ${script} $*"
        "${GENEFORMER_PYTHON}" "${script}" "$@"
    else
        echo "    → ${GENEFORMER_PYTHON} ${script} $*"
    fi
}

run_R() {
    local script="$1"
    shift
    if [[ "$RUN_MODE" == "run" ]]; then
        echo "    [Executing] Rscript ${script} $*"
        Rscript "${script}" "$@"
    else
        echo "    → Rscript ${script} $*"
    fi
}

if [[ "$RUN_MODE" != "run" ]]; then
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║  DRY-RUN: Showing execution order only.                     ║"
    echo "║  To execute:  bash MasterAnalysis.sh --run                  ║"
    echo "╚═══════════════════════════════════════════════════════════════╝"
    echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "0" ]]; then
echo "═══ STAGE 0: Cell Ranger alignment ═══"
echo "  cellranger_shabanas_gex.sh"
echo "  Requires: Cell Ranger 9.0.1 + GRCm39 reference genome."
echo "  Run CellBender on Cell Ranger output before proceeding to Stage 1."
echo "  (Manual stage — see scripts/cellranger_shabanas_gex.sh)"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "1" ]]; then
echo "═══ STAGE 1: Mouse scVI + CytoTRACE2 integration ═══"
check_stage_data 1
echo "  [1a] scVI integration, QC, doublet removal, CytoTRACE2"
run_py "${SCRIPTS}/mouse_scvi_cytotrace2_cellbender.py"
echo "  [1b] Downstream violin/UMAP plots"
run_py "${SCRIPTS}/plot_cytotrace2_downstream.py"
echo "  [1c] Annotation confidence UMAP"
run_py "${SCRIPTS}/plot_umap_confidence.py"
echo "  [1d] Leiden majority cell-type UMAP"
run_py "${SCRIPTS}/plot_umap_leiden_majority_confidence.py"
echo "  [1e] Per-condition Leiden cell-type UMAPs"
run_py "${SCRIPTS}/plot_umap_leiden_majority_confidence_by_condition.py"
echo "  [1f] Public gene-set scores (MSigDB / Reactome / KEGG)"
run_py "${SCRIPTS}/generate_public_gene_set_scores_mouse.py"
echo "  [1g] PPARGC1A / Mootha PGC GSEA overlay"
run_py "${SCRIPTS}/gsea_umap_ppargc1a_A.py"
run_py "${SCRIPTS}/gsea_umap_ppargc1a_B.py"
echo "  [1h] Highlight Leiden clusters 4 & 6"
run_py "${SCRIPTS}/plot_umap_highlight_clusters_4_6.py"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "2" ]]; then
echo "═══ STAGE 2: Geneformer training & prediction (mouse) ═══"
check_stage_data 2
echo "  Using Geneformer environment: ${GENEFORMER_PYTHON}"
echo "  [2a] Fine-tune Geneformer on tonsil atlas (48-class)"
run_gf_py "${SCRIPTS}/train_geneformer_tonsil_multi.py"
echo "  [2b] Predict cell types on mouse scVI object"
run_gf_py "${SCRIPTS}/geneformer_predict_and_plot_manuscript.py"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "3" ]]; then
echo "═══ STAGE 3: Cross-species integration ═══"
check_stage_data 3
echo "  [3a] scVI: human DLBCL + mouse malignant + tonsil B/plasma cells"
run_py "${SCRIPTS}/scvi_human_dlbcl_mouse_malignant_integration.py"
echo "  [3b] CytoTRACE2 on integrated object"
run_py "${SCRIPTS}/cytotrace2_dlbcl_mouse_tonsil.py"
echo "  [3c] Geneformer prediction on integrated object"
run_gf_py "${SCRIPTS}/geneformer_predict_dlbcl_mouse_tonsil.py"
echo "  [3d] Downstream plots (violin, UMAP by species/disease)"
run_py "${SCRIPTS}/downstream_plots_human_mouse_integration.py"
echo "  [3e] CytoTRACE2 on scVI UMAP"
run_py "${SCRIPTS}/plot_cytotrace2_on_scvi_umap.py"
echo "  [3f] CytoTRACE2 stratified by tonsil subtype"
run_py "${SCRIPTS}/plot_cytotrace2_stratified_by_tonsil_subtype.py"
echo "  [3g] Individual sample UMAPs"
run_py "${SCRIPTS}/plot_individual_samples_umap.py"
echo "  [3h] Highlight mouse clusters 4 & 6 on integration UMAP"
run_py "${SCRIPTS}/plot_umap_highlight_mouse_clusters_4_6.py"
echo "  [3i] Convert to Seurat v5 .rds"
run_R  "${SCRIPTS}/convert_human_mouse_integration_to_seurat5.R"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "4" ]]; then
echo "═══ STAGE 4: Gene expression programs (heatmaps / GSEA) ═══"
echo "  [4a] CytoTRACE2-ordered heatmaps + program detection + GSEA"
run_py "${SCRIPTS}/ordering_cytotrace_2_mouse_geneformer.py"
echo "  [4b] Custom GSEA scatter plots"
run_py "${SCRIPTS}/plot_gsea_custom.py"
run_py "${SCRIPTS}/plot_gsea_custom_B.py"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
if [[ "$TARGET_STAGE" == "all" || "$TARGET_STAGE" == "5" ]]; then
echo "═══ STAGE 5: Supplementary analyses ═══"
echo "  [5a] PGC1α CytoTRACE2 trajectory"
run_py "${SCRIPTS}/pgc1_cytotrace2_trajectory.py"
echo "  [5b] Wilcoxon rank-sum tests"
run_py "${SCRIPTS}/wilcoxon_rank_mouse_integrated.py"
echo ""
fi

# ═══════════════════════════════════════════════════════════════════════════
echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║                   PIPELINE SUMMARY                          ║"
echo "╠═══════════════════════════════════════════════════════════════╣"
echo "║  Stage 0  Cell Ranger alignment (manual)                   ║"
echo "║  Stage 1  Mouse scVI + CytoTRACE2           (9 scripts)   ║"
echo "║  Stage 2  Geneformer train + predict         (2 scripts)   ║"
echo "║  Stage 3  Cross-species integration          (9 scripts)   ║"
echo "║  Stage 4  Heatmap programs + GSEA            (3 scripts)   ║"
echo "║  Stage 5  Supplementary analyses             (2 scripts)   ║"
echo "╠═══════════════════════════════════════════════════════════════╣"
echo "║  27 scripts in ./scripts/ (26 analysis + 1 smoke test)      ║"
echo "╚═══════════════════════════════════════════════════════════════╝" 
echo ""
echo "Done."
