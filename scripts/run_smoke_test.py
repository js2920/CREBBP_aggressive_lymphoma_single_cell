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
