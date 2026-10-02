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
INPUT_H5AD = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/integrated_human_dlbcl_mouse_malignant_tonsil.h5ad")

# Geneformer model (trained on tonsil)
GENEFORMER_MODEL_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/models/geneformer_tonsil_multi/fine_tuned_model")
LABEL_ENCODER_PATH = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/models/geneformer_tonsil_multi/label_encoder.pkl")

# Output directory
OUTPUT_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/mouse_human_integration/Geneformer")

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



