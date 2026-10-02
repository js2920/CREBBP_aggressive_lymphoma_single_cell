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


