#!/usr/bin/env python3
"""
Geneformer Prediction and Visualization for MANUSCRIPT
========================================================

This script:
1. Loads the CellBender-filtered mouse data (with scVI integration)
2. Predicts cell types using the pre-trained 48-class Geneformer model
3. Generates publication-quality visualizations


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
MOUSE_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_scvi_cytotrace2/mouse_integrated.h5ad")
GENEFORMER_MODEL_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/models/geneformer_tonsil_multi/fine_tuned_model")
LABEL_ENCODER_PATH = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/models/geneformer_tonsil_multi/label_encoder.pkl")
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/Geneformer")

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



