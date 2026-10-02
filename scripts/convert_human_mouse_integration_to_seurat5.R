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
#   Rscript /home/gusti/CREBBP_aggressive_lymphoma_single_cell/scripts/convert_human_mouse_integration_to_seurat5.R
#
# Author: J
# Date: 2025-12-02

# ============================== PATHS ========================================
INPUT_H5AD <- "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/Geneformer/integrated_with_geneformer_predictions.h5ad"
OUTPUT_RDS <- "/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/Geneformer/integrated_with_geneformer_predictions_seurat5.rds"

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



