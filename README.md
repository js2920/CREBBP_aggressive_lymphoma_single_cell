Complete, reproducible analysis code for the single-cell RNA-seq analyses in the accompanying manuscript on *Crebbp* loss, B-cell lymphoma progression, and cross-species comparison with human DLBCL.

This is a self-contained reproducibility recipe. It embeds all 26 analysis scripts, extracts them into `./scripts/`, and can run the full pipeline end-to-end.

## What this repository provides

The code documents the statistical framework underlying the single-cell figures in the manuscript, including:

- Quality control and biological-replicate integration of mouse samples  
- Cell-state annotation using a fine-tuned Geneformer model (human tonsil atlas)  
- Cross-species integration of mouse malignant B cells with human DLBCL and tonsil germinal-centre B cells  
- CytoTRACE2 differentiation / stemness scoring along malignant progression  
- Gene-expression programs and GSEA along the CytoTRACE2 axis  
- Cluster-level marker testing (Wilcoxon rank-sum)

---

## Pipeline overview

| Stage | Description | Key tools |
|-------|-------------|-----------|
| 0 | Alignment to GRCm39 | Cell Ranger 9.0.1 |
| 1 | Mouse scVI integration, CytoTRACE2, gene-set scoring | scVI, CytoTRACE2, scanpy |
| 2 | Geneformer fine-tuning on human tonsil atlas → mouse prediction | Geneformer |
| 3 | Cross-species scVI (human DLBCL + mouse + tonsil GC B cells) | scVI, CytoTRACE2, Geneformer |
| 4 | Gene-expression programs along CytoTRACE2 + GSEA | scipy, gseapy |
| 5 | Supplementary analyses (PGC1α trajectory, Wilcoxon rank-sum) | scanpy, scipy |

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/js2920/CREBBP_aggressive_lymphoma_single_cell.git
cd CREBBP_aggressive_lymphoma_single_cell

# 2. Edit the USER CONFIGURATION section at the top of MasterAnalysis.sh
#    so that data paths match your local layout (see Data setup below)
nano MasterAnalysis.sh

# 3. Extract the 26 analysis scripts (dry-run; does not execute the pipeline)
bash MasterAnalysis.sh

# 4. Run the full pipeline, or execute individual stages from ./scripts/
bash MasterAnalysis.sh --run


---

## Mouse samples

| Condition | Genotype | Replicates | Description |
|-----------|----------|------------|-------------|
| WT_B_cells | Wild-type | 2 | Control B cells |
| Crebbp_B_cells | Crebbp-KO | 2 | Crebbp-deficient B cells |
| Pre_malignant | Crebbp-WT | 2 | Pre-malignant B cells |
| Malignant | Crebbp-KO | 2 | Malignant B cells |
| Matched_malignant | Crebbp-KO | 2 | Matched malignant B cells |

---

## Data setup

### Generated mouse data

Raw and processed mouse scRNA-seq data are deposited in GEO.


GEO accession: GSE332767


Place CellBender-filtered count matrices under `data/`:


data/
├── cellbender_filtered/        # CellBender-filtered .h5 files (10 samples)
│   ├── WT_B_cells_R1_filtered.h5
│   ├── WT_B_cells_R2_filtered.h5
│   ├── Crebbp_B_cells_R1_filtered.h5
│   ├── Crebbp_B_cells_R2_filtered.h5
│   ├── Pre_malignant_R1_filtered.h5
│   ├── Pre_malignant_R2_filtered.h5
│   ├── Malignant_R1_filtered.h5
│   ├── Malignant_R2_filtered.h5
│   ├── Matched_malignant_R1_filtered.h5
│   └── Matched_malignant_R2_filtered.h5
├── raw_fastq/                  # Optional: raw FASTQs if re-running Cell Ranger
└── refdata-gex-GRCm39-2024-A/  # Optional: 10x GRCm39 reference (Stage 0 only)


### Public datasets used for cross-species integration (Stage 3)

| Dataset | Source | Expected location |
|---------|--------|-------------------|
| Roider et al. DLBCL | Publication-associated scRNA-seq | `data/DLBCL/DLBCL1_raw.h5ad`, `DLBCL2_raw.h5ad`, `DLBCL3_raw.h5ad` |
| Alizadeh et al. DLBCL (CD20+) | GEO [GSE182434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182434) | `data/DLBCL/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad` |
| Human tonsil atlas | [HCATonsilData](https://bioconductor.org/packages/HCATonsilData) (R/Bioconductor) | `data/tonsil_export/tonsil_GCBC_RNA.h5ad`, `tonsil_NBC-MBC_RNA.h5ad` |

Export tonsil atlas subsets from R:

```r
library(HCATonsilData)
library(zellkonverter)

tonsil_gc <- HCATonsilData(assayType = "RNA", cellType = "GCBC")
writeH5AD(tonsil_gc, "data/tonsil_export/tonsil_GCBC_RNA.h5ad")

tonsil_mbc <- HCATonsilData(assayType = "RNA", cellType = "NBC")
# or NBC-MBC export as used in the manuscript methods
writeH5AD(tonsil_mbc, "data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad")
```

---

## Software requirements

### Python (≥ 3.10)

```bash
pip install "scanpy>=1.9" "scvi-tools>=1.0" anndata pandas numpy scipy \
    matplotlib seaborn scrublet mygene gseapy magic-impute
pip install cytotrace2_py   # https://github.com/digitalcytometry/cytotrace2
```

### Geneformer (recommended: separate environment)

```bash
conda create -n geneformer python=3.10
conda activate geneformer
pip install transformers datasets accelerate scikit-learn
# Install Geneformer from: https://huggingface.co/ctheodoris/Geneformer
```

### R (≥ 4.0)

```r
install.packages(c("Seurat", "SeuratObject"))
BiocManager::install(c("zellkonverter", "SingleCellExperiment", "HCATonsilData"))
```

### Cell Ranger (Stage 0 only)

Cell Ranger 9.0.1 from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger).

---

## Repository layout


├── MasterAnalysis.sh          # Self-contained pipeline (run this)
├── README.md                  # This file
└── scripts/                   # Created by: bash MasterAnalysis.sh
    ├── mouse_scvi_cytotrace2_cellbender.py
    ├── train_geneformer_tonsil_multi.py
    ├── geneformer_predict_and_plot_manuscript.py
    ├── scvi_human_dlbcl_mouse_malignant_integration.py
    ├── ordering_cytotrace_2_mouse_geneformer.py
    ├── wilcoxon_rank_mouse_integrated.py
    └── ...                    # 26 scripts in total


After a successful run, major outputs are written to:

| Directory | Stage |
|-----------|-------|
| `mouse_scvi_cytotrace2/` | Stage 1 |
| `Geneformer/` | Stage 2 |
| `mouse_human_integration/` | Stage 3 |
| `heatmap_programs/` | Stage 4 |
| `pgc1_cytotrace2_trajectory/` | Stage 5a |
| `wilcoxon_rank_sum/` | Stage 5b |

---

## Path configuration

Edit the **USER CONFIGURATION** block at the top of `MasterAnalysis.sh`:

| Variable | Description | Default |
|----------|-------------|---------|
| `BASEDIR` | Root directory (auto-detected) | Directory containing `MasterAnalysis.sh` |
| `CELLBENDER_DIR` | CellBender-filtered h5 files | `${BASEDIR}/data/cellbender_filtered` |
| `DLBCL_DIR` | Human DLBCL h5ad files | `${BASEDIR}/data/DLBCL` |
| `TONSIL_DIR` | Human tonsil atlas h5ad exports | `${BASEDIR}/data/tonsil_export` |
| `GENEFORMER_MODEL_DIR` | Fine-tuned Geneformer model | `${BASEDIR}/models/geneformer_tonsil_multi` |
| `FASTQ_DATA_ROOT` | Raw FASTQ archives (Stage 0) | `${BASEDIR}/data/raw_fastq` |
| `CELLRANGER_REF` | GRCm39 reference (Stage 0) | `${BASEDIR}/data/refdata-gex-GRCm39-2024-A` |
| `CELLRANGER_BIN` | Cell Ranger binary | `cellranger` |

When scripts are extracted, paths are patched automatically to match this configuration.

---

## Running individual stages

```bash
# Extract scripts only
bash MasterAnalysis.sh

# Stage 1 — mouse integration + CytoTRACE2
python3 scripts/mouse_scvi_cytotrace2_cellbender.py

# Stage 2 — Geneformer (use the geneformer conda env)
conda activate geneformer
python3 scripts/train_geneformer_tonsil_multi.py
python3 scripts/geneformer_predict_and_plot_manuscript.py

# Stage 3 — cross-species integration
python3 scripts/scvi_human_dlbcl_mouse_malignant_integration.py

# Stage 4 — CytoTRACE2 programs / heatmaps
python3 scripts/ordering_cytotrace_2_mouse_geneformer.py

# Stage 5 — cluster markers
python3 scripts/wilcoxon_rank_mouse_integrated.py
```

---


## Citation not availabe yet

Manuscript citation will be added upon publication / acceptance.

## Contact

For questions about the analysis code, open an issue in this repository.
