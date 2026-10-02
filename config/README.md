# Configuration Guide

The `config/pipeline_config.yaml` file centralizes all tunable parameters, file paths, QC thresholds, model hyperparameters, and statistical settings across the 6 pipeline stages.

## Structure

1. **`pipeline`**: Global random seeds, number of CPU threads for OpenMP/MKL/NUMEXPR.
2. **`paths`**: Configurable root and data directories. By default, paths are relative to the repository base directory.
3. **`stage0_cellranger`**: Compute resources and reference paths for Cell Ranger alignment.
4. **`stage1_mouse_scvi`**: Quality control filters (mitochondrial %, minimum genes), Scrublet doublet simulation settings, scVI generative model hyperparameters (latent dimension $d=30$, Negative Binomial likelihood, 2-layer encoder/decoder), Leiden clustering resolutions, and CytoTRACE2 configuration.
5. **`stage2_geneformer`**: Transformer hyperparameter settings (token context length 2048, learning rate $5 \times 10^{-5}$, macro-F1 early stopping).
6. **`stage3_cross_species`**: Ensembl BioMart 1:1 orthology constraints, human tonsil GC B cell subset selection, and joint multi-species scVI conditioning.
7. **`stage4_trajectory_programs`**: MAGIC graph diffusion operator ($t=3$), CytoTRACE2 binning ($N=50$), Gaussian smoothing ($\sigma=1.5$), and dual-metric trajectory distance parameters ($\alpha_{\mathrm{shape}} = 0.5$).
8. **`stage5_differential_expression`**: Non-parametric Wilcoxon rank-sum test thresholds and Benjamini-Hochberg FDR cutoffs.

To override settings without modifying code, edit `pipeline_config.yaml` or pass custom paths via environment variables (e.g., `export CREBBP_BASE_DIR=/custom/path`).
