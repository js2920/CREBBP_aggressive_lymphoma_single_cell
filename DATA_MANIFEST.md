# Data Manifest & Availability

This document provides  metadata, accessions, and download instructions for all primary mouse sequencing data and external reference datasets used in the study:

> **"Epigenetic deregulation and metabolic reprogramming driven by Crebbp loss in aggressive B-cell lymphoma"**

---

## 1. Primary Mouse Single-Cell RNA-seq Datasets (GEO: GSE332767)

Ten single-cell RNA sequencing libraries from murine germinal center and transformed lymphoma states were generated using the 10x Genomics Chromium Single Cell 3' Gene Expression platform (v3 chemistry) and sequenced on an Illumina NovaSeq 6000.

| Sample ID | GEO Accession | Genotype | Biological Condition | Replicate | CellBender Filtered File |
|-----------|---------------|----------|----------------------|-----------|--------------------------|
| SIGAC6 | GSMxxxxxxx | Wild-type (*Crebbp* +/+) | WT control B cells | R1 | `SIGAC6_WT_B_cells_R1_GEX_cellbender_filtered.h5` |
| SIGAB6 | GSMxxxxxxx | Wild-type (*Crebbp* +/+) | WT control B cells | R2 | `SIGAB6_WT_B_cells_R2_GEX_cellbender_filtered.h5` |
| SIGAA6 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Crebbp-deficient B cells | R1 | `SIGAA6_Crebbp_B_cells_R1_GEX_cellbender_filtered.h5` |
| SIGAD6 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Crebbp-deficient B cells | R1_2 | `SIGAD6_Crebbp_B_cells_R1_2_GEX_cellbender_filtered.h5` |
| SIGAE6 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Crebbp-deficient B cells | R2 | `SIGAE6_Crebbp_B_cells_R2_GEX_cellbender_filtered.h5` |
| SIGAC2 | GSMxxxxxxx | *Crebbp* wild-type | Pre-malignant B cells | R2 | `SIGAC2_Pre_malignant_R2_GEX_cellbender_filtered.h5` |
| SIGAD5 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Malignant lymphoma | R2 | `SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5` |
| SIGAE5 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Malignant lymphoma | R1 | `SIGAE5_Malignant_R1_GEX_cellbender_filtered.h5` |
| SIGAA3 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Matched malignant lymphoma | R1 | `SIGAA3_Matched_malignant_R1_GEX_cellbender_filtered.h5` |
| SIGAA4 | GSMxxxxxxx | *Crebbp* knockout (*Crebbp* -/-) | Matched malignant lymphoma | R2 | `SIGAA4_Matched_malignant_R2_GEX_cellbender_filtered.h5` |

### Expected Directory Layout
Place the CellBender-filtered `.h5` files in `data/cellbender_filtered/`:
```text
data/cellbender_filtered/
├── SIGAA3_Matched_malignant_R1_GEX_cellbender_filtered.h5
├── SIGAA4_Matched_malignant_R2_GEX_cellbender_filtered.h5
├── SIGAA6_Crebbp_B_cells_R1_GEX_cellbender_filtered.h5
├── SIGAB6_WT_B_cells_R2_GEX_cellbender_filtered.h5
├── SIGAC2_Pre_malignant_R2_GEX_cellbender_filtered.h5
├── SIGAC6_WT_B_cells_R1_GEX_cellbender_filtered.h5
├── SIGAD5_Malignant_R2_GEX_cellbender_filtered.h5
├── SIGAD6_Crebbp_B_cells_R1_2_GEX_cellbender_filtered.h5
├── SIGAE5_Malignant_R1_GEX_cellbender_filtered.h5
└── SIGAE6_Crebbp_B_cells_R2_GEX_cellbender_filtered.h5
```

---

## 2. External Public Datasets (Stage 3 Cross-Species Integration)

### A. Human DLBCL Cohort 1 (Roider et al., 2021)
- **Source**: Primary patient biopsies of diffuse large B-cell lymphoma.
- **Expected files**:
  - `data/DLBCL/DLBCL1_raw.h5ad`
  - `data/DLBCL/DLBCL2_raw.h5ad`
  - `data/DLBCL/DLBCL3_raw.h5ad`

### B. Human DLBCL Cohort 2 (Alizadeh et al., GEO GSE182434)
- **Source**: CD20+ enriched single cells from DLBCL patient biopsies.
- **Expected file**: `data/DLBCL/Alizadeh/GSE182434_DLBCL_CD20pos.h5ad`

### C. Human Tonsil Atlas (Human Cell Atlas / Massoni-Badosa et al.)
- **Source**: Bioconductor package `HCATonsilData`.
- **Expected files**:
  - `data/tonsil_export/tonsil_GCBC_RNA.h5ad` (Germinal center B cells: light zone, dark zone, centrocytes, centroblasts)
  - `data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad` (Naive B cells & Memory B cells)
- **Extraction Command**:
  ```r
  library(HCATonsilData)
  library(zellkonverter)
  dir.create("data/tonsil_export", recursive = TRUE, showWarnings = FALSE)
  writeH5AD(HCATonsilData(assayType = "RNA", cellType = "GCBC"), "data/tonsil_export/tonsil_GCBC_RNA.h5ad")
  writeH5AD(HCATonsilData(assayType = "RNA", cellType = "NBC"), "data/tonsil_export/tonsil_NBC-MBC_RNA.h5ad")
  ```

---

## 3. Pretrained Foundation Models (Stage 2 Geneformer)

- **Base Architecture**: 6-layer or 12-layer Geneformer foundation model from Hugging Face (`ctheodoris/Geneformer`).
- **Fine-tuned Checkpoint**: Fine-tuned on the human tonsil atlas covering 48 immune cell phenotypes.
- **Expected Directory**: `models/geneformer_tonsil_multi/`
  - `fine_tuned_model/config.json`
  - `fine_tuned_model/pytorch_model.bin` (or `model.safetensors`)
  - `label_encoder.pkl`
