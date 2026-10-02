# Mathematical Methods and Analytical Formulations

This document provides a concise, formal specification of the mathematical models, statistical formulations, and analytical assumptions implemented across the *Crebbp*-deficient aggressive B-cell lymphoma single-cell analysis pipeline.

---

## 1. Study Design & Cohort Architecture

The analytical framework investigates cellular heterogeneity, developmental differentiation potency, and transcriptional dynamics across murine aggressive B-cell lymphoma models and human patient cohorts.

### 1.1 Experimental Conditions and Sample Hierarchy

The murine single-cell RNA sequencing cohort comprises ten discrete libraries spanning five biological conditions:
- **`WT_B_cells`**: Wild-type physiological B-cell controls (`SIGAC6` [R1], `SIGAB6` [R2]).
- **`Crebbp_B_cells`**: Non-transformed *Crebbp*-deficient B cells (`SIGAA6` [R1], `SIGAD6` [R1_2]).
- **`Pre_malignant`**: Early-stage pre-malignant *Crebbp*-deficient B cells (`SIGAF2` [R1], `SIGAC2` [R2]).
- **`Malignant`**: Frank aggressive B-cell lymphoma specimens (`SIGAH1` [R1], `SIGAD5` [R2]).
- **`Matched_malignant`**: Serial matched malignant lymphoma recurrences (`SIGAA3` [R1], `SIGAA4` [R2]).

### 1.2 Multi-Tier Analytical Architecture

The computational workflow consists of four interconnected modeling domains:
1. **Count Decontamination & Latent Space Inference**: Ambient RNA clearance via CellBender [1] and deep generative variational inference via scVI [2, 3] to resolve non-linear transcriptional manifolds.
2. **Cross-Species Reference Alignment**: Reciprocal orthology mapping bridging murine lymphoma transcriptomes with human DLBCL patient profiles and healthy reactive tonsil B and plasma cell populations.
3. **Differentiation Potency Modeling**: Computational scoring of developmental arrest and stem-like plasticity using regularized graph diffusion via CytoTRACE2 [7, 8].
4. **Dynamic Gene Program Discovery**: Trajectory-aligned feature selection, Markov affinity graph denoising via MAGIC [10], and dual-metric hierarchical clustering to isolate condition-specific transcriptional programs.

---

## 2. Count Processing, Decontamination, and Quality Control

### 2.1 Ambient RNA Decontamination Model

Observed droplet expression profiles reflect a mixture of intracellular RNA and cell-free ambient contamination. CellBender [1] formulates observed count $y_{cg}$ for cell barcode $c$ and gene $g$ via a hierarchical Bayesian model:

$$
y_{cg} = z_{cg} + a_{cg}
$$

where $z_{cg}$ represents the true intracellular transcriptional signal and $a_{cg}$ represents ambient background contamination:

$$
a_{cg} \sim \mathrm{Poisson}(\rho_c^{\mathrm{amb}} \cdot \chi_g)
$$

with droplet ambient scale parameter $\rho_c^{\mathrm{amb}}$ and ambient gene profile distribution $\chi_g$. The clean decontaminated count matrix is derived from the posterior expectation:

$$
\hat{z}_{cg} = \mathbb{E}_{q_\phi}[z_{cg} \mid y_{cg}]
$$

These background-corrected counts populate `layers['counts']` as raw integer-like inputs for downstream generative modeling.

### 2.2 Quality Control Filtering

Cellular barcodes are filtered to remove non-viable cells, broken droplets, and multi-cell multiplets:
- **Gene complexity**: Minimum of $N_{\mathrm{genes}} \ge 200$ distinct detected genes per cell barcode.
- **Mitochondrial count fraction**: Maximum mitochondrial threshold $\le 10.0\%$.
- **Multiplet removal**: Scrublet [4] detects artificial transcriptomic doublets using simulated homotypic/heterotypic doublets, flagging droplets exceeding an expected doublet rate of 0.06.
- **Gene prevalence**: Genes must be detected in at least $N_{\mathrm{cells}} \ge 3$ cells across the cohort.

The mitochondrial count percentage for cell barcode $c$ is evaluated as:

$$
\mathrm{pct}_{\mathrm{MT}, c} = \frac{\sum_{g \in \mathcal{G}_{\mathrm{MT}}} x_{cg}}{\sum_g x_{cg}} \times 100 \le 10.0\%
$$

### 2.3 Library Size Normalization and Transformation

For downstream visualization and correlation calculations, decontaminated counts are scaled by total library depth $L_c = \sum_g x_{cg}$ and log-transformed:

$$
y_{cg} = \log\left(1 + 10^4 \cdot \frac{x_{cg}}{\max(L_c, \varepsilon)}\right)
$$

### 2.4 Feature Selection and Confounder Masking

To identify biologically informative manifold axes, 5,000 highly variable genes (HVGs) are selected using the `seurat_v3` variance-stabilizing transformation across batch strata in Scanpy [5]:

$$
\sigma^2_g = \mathrm{Var}_{\mathrm{norm}}(x_{\cdot g})
$$

To prevent technical artifacts and non-malignant lineage dynamics from confounding manifold construction, a confounder mask $\mathcal{M}_{\mathrm{conf}}$ excludes:
- Mitochondrial genes (`mt-`, `MT-`)
- Ribosomal protein genes (`Rps`, `Rpl`, `RPS`, `RPL`)
- Variable immunoglobulin chains (`Ighv`, `Igkv`, `Iglv`, `IGHV`, `IGKV`, `IGLV`)
- T-cell receptor chains (`Trav`, `Trbv`, `Trgv`, `Trdv`, `Trac`, `Trbc`, `Trgc`, `Trdc`)

The retained HVG set $\mathcal{H} = \mathcal{G}_{\mathrm{HVG}} \setminus \mathcal{M}_{\mathrm{conf}}$ provides the input feature space for latent representation learning.

---

## 3. Deep Generative Count Modeling (scVI)

### 3.1 Observation Likelihood

The single-cell Variational Inference (scVI) framework [2, 3] models count variation using a deep generative model parameterized by neural networks. For cell $c$ and gene $g \in \mathcal{H}$, the observed count $x_{cg}$ follows a Negative Binomial distribution:

$$
x_{cg} \mid \mathbf{z}_c, L_c, b_c, \mathbf{u}_c \sim \mathrm{NB}(\mu_{cg}, \theta_{g,b_c})
$$

The expected count $\mu_{cg}$ couples total library size $L_c$ with the normalized gene rate $\rho_{cg}$:

$$
\mu_{cg} = L_c \cdot \rho_g(\mathbf{z}_c, b_c, \mathbf{u}_c)
$$

with variance function:

$$
\mathrm{Var}(x_{cg} \mid \mathbf{z}_c, L_c, b_c, \mathbf{u}_c) = \mu_{cg} + \frac{\mu_{cg}^2}{\theta_{g,b_c}}
$$

where:
- $\mathbf{z}_c \in \mathbb{R}^{96}$ represents the low-dimensional continuous latent state of cell $c$.
- $L_c$ is the observed library depth factor.
- $b_c$ is the batch/library covariate (`sample_id`).
- $\mathbf{u}_c$ represents observed biological covariates, including experimental condition and continuous cell-cycle scores ($S$ and $G_2/M$).
- $\theta_{g,b_c} > 0$ denotes gene- and batch-specific inverse dispersion.

### 3.2 Variational Inference and Optimization

The posterior distribution over latent variables $p(\mathbf{z}_c \mid \mathbf{x}_c, b_c, \mathbf{u}_c)$ is approximated using a variational inference encoder network:

$$
q_\phi(\mathbf{z}_c \mid \mathbf{x}_c, b_c, \mathbf{u}_c) = \mathcal{N}\left(\boldsymbol{\mu}_\phi(\mathbf{x}_c), \mathrm{diag}\left(\boldsymbol{\sigma}_\phi^2(\mathbf{x}_c)\right)\right)
$$

Model parameters $(\theta, \phi)$ are trained jointly by maximizing the Evidence Lower Bound (ELBO):

$$
\mathcal{L}_{\mathrm{ELBO}}(\theta, \phi) = \sum_{c} \left( \mathbb{E}_{q_\phi}\left[ \log p_\theta(\mathbf{x}_c \mid \mathbf{z}_c, L_c, b_c, \mathbf{u}_c) \right] - D_{\mathrm{KL}}\left( q_\phi(\mathbf{z}_c \mid \mathbf{x}_c, b_c, \mathbf{u}_c) \,\|\, p(\mathbf{z}_c) \right) \right)
$$

with isotropic Gaussian prior $p(\mathbf{z}_c) = \mathcal{N}(\mathbf{0}, \mathbf{I}_{96})$.

**Architectural hyperparameters:**
- Latent dimensionality: $\dim(\mathbf{z}) = 96$
- Fully connected layers: 2 hidden layers with 128 units each
- Dropout rate: $0.10$
- Normalization: Layer normalization applied across both encoder and decoder networks
- Optimization: AdamW with learning rate $\eta = 10^{-3}$, plateau learning rate decay, and early stopping.

### 3.3 Graph Topology and Leiden Clustering

The low-dimensional latent coordinates $\mathbf{Z} = [\mathbf{z}_1, \ldots, \mathbf{z}_N]^{\mathsf{T}}$ define cell-cell similarity. A $k$-nearest-neighbor affinity graph ($k=30$) is constructed using Euclidean distances in $\mathbb{R}^{96}$:

$$
d_{\mathrm{scVI}}(c_i, c_j) = \|\mathbf{z}_{c_i} - \mathbf{z}_{c_j}\|_2
$$

Cell communities are partitioned using the Leiden algorithm [6] by maximizing modularity across resolution parameters $\gamma \in \{0.3, 0.5, 1.0\}$:

$$
\mathcal{Q}(\gamma) = \frac{1}{2m} \sum_{ij} \left[ A_{ij} - \gamma \frac{k_i k_j}{2m} \right] \delta(\sigma_i, \sigma_j)
$$

---

## 4. Cross-Species Manifold Alignment

### 4.1 Orthology Projection

To project murine lymphoma phenotypes onto human clinical lymphoma and tonsil differentiation hierarchies, mouse gene identifiers are mapped to human orthologs using reciprocal 1:1 Ensembl BioMart homology [13]:

$$
\mathcal{M}_{1:1}: \mathcal{G}_{\mathrm{mouse}} \longrightarrow \mathcal{G}_{\mathrm{human}}
$$

Where Ensembl queries encounter network latency, the MGI vertebrate homology report (`HOM_MouseHumanSequence.rpt`) serves as a deterministic 1:1 mapping fallback. Unmapped or non-orthologous features are removed from the shared cross-species feature matrix.

### 4.2 Joint Latent Integration

Shared highly variable genes are identified by intersecting genes exhibiting high variance across both murine and human datasets:

$$
\mathcal{H}_{\mathrm{shared}} = \left\lbrace g \in \mathcal{G}_{\mathrm{orth}} \mid \sum_{d=1}^D \mathbb{I}(g \in \mathrm{HVG}_d) \ge 2 \right\rbrace
$$

A joint scVI variational autoencoder integrates mouse lymphoma, human DLBCL (primary clinical biopsies), and healthy human tonsil B and plasma cells using a nested batch key:

$$
\mathbf{s}_c = [\mathbf{s}_c^{\mathrm{batch}}, \; \mathbf{s}_c^{\mathrm{species}}]
$$

This formulation models technical batch and species divergence while aligning shared biological cell states along the common 96-dimensional manifold $\mathbf{z}_c$.

---

## 5. Developmental Potency Modeling (CytoTRACE2)

### 5.1 Transcriptional Entropy and Potency Formulation

Developmental potency and differentiation commitment are scored using CytoTRACE2 [7, 8]. The underlying formulation exploits the relationship between open chromatin plasticity and transcriptional entropy, quantified by the active gene count signature:

$$
\mathrm{GCS}_c = \sum_{g} \mathbb{I}(x_{cg} > 0)
$$

### 5.2 Harmonic Graph Diffusion

Initial potency scores $\mathbf{S}^{(0)}$ derived from pre-trained gene-expression models are propagated and smoothed over the single-cell affinity graph $\mathbf{W}$ via regularized harmonic diffusion:

$$
\mathbf{S}^* = (1 - \alpha)\left(\mathbf{I} - \alpha \mathbf{P}\right)^{-1} \mathbf{S}^{(0)}
$$

where $\mathbf{P} = \mathbf{D}^{-1}\mathbf{W}$ represents the row-stochastic random-walk transition matrix with degree matrix $D_{ii} = \sum_j W_{ij}$, and $\alpha = 0.9$ is the Markov restart parameter.

### 5.3 Directionality and Differentiation Pseudotime

The standardized CytoTRACE2 score $q_c \in [0, 1]$ directly reflects developmental potency, where higher values indicate stem-like, plastic, or undifferentiated phenotypes. Differentiation pseudotime $\tau_c$ is defined as the reverse coordinate:

$$
\tau_c = 1 - q_c
$$

tracking progressive maturation from undifferentiated stem-like B-cell states ($\tau_c \to 0$, $q_c \to 1$) toward terminally differentiated effector states ($\tau_c \to 1$, $q_c \to 0$).

---

## 6. Context-Aware Foundation Model Transfer (Geneformer)

### 6.1 Rank-Value Transcriptome Tokenization

Geneformer [9] processes single-cell transcriptomes as rank-ordered token sequences. For cell $c$, non-zero gene expressions are normalized by total cell depth $L_c$:

$$
\tilde{x}_{cg} = 10^4 \cdot \frac{x_{cg}}{L_c}
$$

and sorted in descending order of relative abundance:

$$
\boldsymbol{\tau}_c = \left( t_{(1)}, t_{(2)}, \ldots, t_{(M)} \right), \quad \tilde{x}_{c, t_{(1)}} \ge \tilde{x}_{c, t_{(2)}} \ge \cdots \ge \tilde{x}_{c, t_{(M)}}
$$

This rank-value encoding provides mathematical scale invariance:

$$
\boldsymbol{\tau}_c(\kappa \cdot \mathbf{x}_c) = \boldsymbol{\tau}_c(\mathbf{x}_c) \quad \forall \; \kappa > 0
$$

making the input robust to variation in technical sequencing depth. The sequence is truncated to the model context length of $M = 2,048$ tokens.

### 6.2 Fine-Tuning on Tonsil B and Plasma Cell Reference States

A pre-trained 6-layer transformer is fine-tuned for cell state classification using annotated human tonsil B and plasma cells (dark zone centroblasts, light zone centrocytes, memory B cells, and plasma cells) [9]. The classification objective optimizes cross-entropy loss:

$$
\mathcal{L}_{\mathrm{CE}} = - \frac{1}{N} \sum_{i=1}^N \sum_{k=1}^K y_{ik} \log p_{ik}
$$

with AdamW optimizer ($\eta = 5 \times 10^{-5}$, weight decay $0.01$, mixed-precision FP16).

### 6.3 Projection and Uncertainty Quantification

Mouse lymphoma cells mapped to human orthologs are tokenized and passed through the fine-tuned classifier. Predicted probabilities for class $k$ are computed via softmax:

$$
p_{ck} = \frac{\exp(a_{ck})}{\sum_{j=1}^K \exp(a_{cj})}, \quad \hat{k}_c = \arg\max_k p_{ck}
$$

Prediction confidence is quantified using the maximum posterior probability:

$$
C_c = \max_k p_{ck}
$$

High-confidence predictions ($C_c \ge 0.70$) define cells exhibiting confident phenotypic alignment with canonical tonsil B and plasma cell maturation compartments.

---

## 7. Dynamic Trajectory Modeling & Gene Program Discovery

### 7.1 Condition-Discriminative Feature Selection

To isolate gene programs whose dynamics are modulated by *Crebbp* loss, features are selected based on both overall potency correlation and condition-specific interaction:

**1. Baseline Filter:**
Genes are retained with detection fraction $\ge 3\%$ and mean normalized expression $\ge 0.05$ CP10k.

**2. Global Potency Correlation:**
The correlation of gene expression along the CytoTRACE2 potency coordinate $q$ is computed as:

$$
r_g = \mathrm{Spearman}\left(\{y_{cg}\}_{c=1}^N, \{q_c\}_{c=1}^N\right)
$$

**3. Condition-Interaction Score:**
Within each condition $k \in \mathcal{K}$ containing at least 10 cells, the within-condition Spearman correlation $r_{gk}$ is computed. The interaction score $I_g$ measures divergence across conditions:

$$
I_g = \sqrt{\frac{1}{K_g} \sum_{k=1}^{K_g} (r_{gk} - \bar{r}_g)^2}, \quad \bar{r}_g = \frac{1}{K_g} \sum_{k=1}^{K_g} r_{gk}
$$

**4. Composite Ranking Score:**
Min-max normalized global correlation $\mathcal{N}(|r_g|)$ and interaction score $\mathcal{N}(I_g)$ are combined:

$$
A_g = 0.5 \cdot \mathcal{N}(|r_g|) + 0.5 \cdot \mathcal{N}(I_g)
$$

The top candidate genes ($N = 3,000$) proceed to trajectory profile construction.

### 7.2 Markov Affinity Graph Diffusion (MAGIC)

To address technical dropouts and expose continuous expression dynamics along the trajectory, normalized expression profiles undergo Markov affinity diffusion via MAGIC [10]:

$$
\widehat{\mathbf{X}} = \mathbf{M}^t \mathbf{X}
$$

where $\mathbf{M}$ is the row-normalized adaptive Gaussian kernel transition matrix computed in PCA space ($\dim = 100$, $k = 30$ neighbors, decay parameter $1.0$), and $t = 12$ diffusion time steps provide low-pass graph spectral filtering.

### 7.3 Dual-Metric Trajectory Distance

Cells are ordered along ascending CytoTRACE2 potency $q_c$. Expression matrices are smoothed across cell rank using a Gaussian kernel ($\sigma = 3.0$) and rolling-window averaging ($W = 100$ cells).

For genes $g_1$ and $g_2$ with smoothed expression profiles $\mathbf{v}_1$ and $\mathbf{v}_2$, trajectory dissimilarity combines temporal shape with absolute expression amplitude:

$$
D_{\mathrm{dual}}(g_1, g_2) = \alpha_{\mathrm{shape}} \cdot d_{\mathrm{shape}}(g_1, g_2) + (1 - \alpha_{\mathrm{shape}}) \cdot d_{\mathrm{amp}}(g_1, g_2)
$$

The parameter $\alpha_{\mathrm{shape}} = 0.5$ balances dynamic profile shape and absolute expression magnitude equally.

**Trajectory Shape Distance:**
Pearson correlation distance is evaluated on standardized $z$-score profiles:

$$
d_{\mathrm{shape}}(g_1, g_2) = 1 - \frac{\sum_{i=1}^N (z_{1i} - \bar{z}_1)(z_{2i} - \bar{z}_2)}{\sqrt{\sum_{i=1}^N (z_{1i} - \bar{z}_1)^2 \cdot \sum_{i=1}^N (z_{2i} - \bar{z}_2)^2}}
$$

yielding $d_{\mathrm{shape}} \in [0, 2]$.

**Trajectory Amplitude Distance:**
Euclidean distance is normalized by the empirical 95th percentile scale factor $S_{95} = \mathrm{Percentile}_{95}(d_{\mathrm{raw}})$:

$$
d_{\mathrm{amp}}(g_1, g_2) = \min\left(2.0, \quad 2.0 \cdot \frac{\|\mathbf{v}_1 - \mathbf{v}_2\|_2}{S_{95}}\right)
$$

yielding $d_{\mathrm{amp}} \in [0, 2]$.

### 7.4 Hierarchical Program Clustering

Agglomerative hierarchical clustering with Ward's linkage is performed on $D_{\mathrm{dual}}$:

$$
\Delta d(u, v) = \frac{n_u n_v}{n_u + n_v} \|\mathbf{m}_u - \mathbf{m}_v\|_2^2
$$

with optimal leaf ordering applied to minimize adjacent leaf dissimilarity. Discovered gene programs are ordered along the developmental axis by program center of mass:

$$
\mathrm{CoM}_g = \frac{\sum_{b=1}^B b \cdot \bar{y}_{gb}}{\sum_{b=1}^B \bar{y}_{gb}}
$$

---

## 8. Statistical Inference & Functional Enrichment

### 8.1 Non-Parametric Marker Testing

Cluster-defining marker genes are identified using the two-sided Wilcoxon rank-sum test with continuity correction. For gene $g$ between group 1 ($n_1$ cells) and group 2 ($n_2$ cells):

$$
U = \sum_{i=1}^{n_1} \sum_{j=1}^{n_2} \left[ \mathbb{I}(y_{1i} > y_{2j}) + \frac{1}{2} \mathbb{I}(y_{1i} = y_{2j}) \right]
$$

Under the null hypothesis, the asymptotic standardized test statistic is:

$$
Z = \frac{U - \frac{n_1 n_2}{2}}{\sqrt{\mathrm{Var}(U)}}
$$

with tie-adjusted variance:

$$
\mathrm{Var}(U) = \frac{n_1 n_2}{12}\left(N + 1 - \frac{\sum_{t} (t^3 - t)}{N(N - 1)}\right)
$$

### 8.2 Multiple Hypothesis Correction

To control the False Discovery Rate across $M$ tested genes, the Benjamini-Hochberg (BH) step-up procedure is applied:

$$
P_{(1)} \le P_{(2)} \le \cdots \le P_{(M)}
$$

$$
P_{\mathrm{adj},(i)} = \min_{j \ge i} \left( \min\left(1, \; \frac{M}{j} P_{(j)}\right) \right)
$$

Statistical significance is defined at an FDR threshold of $P_{\mathrm{adj}} < 0.05$.

### 8.3 Preranked Gene Set Enrichment Analysis (GSEA)

Functional pathway enrichment along the potency trajectory is evaluated using preranked GSEA [11] implemented via GSEApy [12]. Genes are ordered by their Spearman potency correlation score:

$$
s_g = r_g = \mathrm{Spearman}(y_{\cdot g}, q_\cdot)
$$

The enrichment score $\mathrm{ES}(\mathcal{S})$ for gene set $\mathcal{S}$ evaluates the maximum deviation of a weighted running-sum statistic:

$$
P_{\mathrm{hit}}(i) = \sum_{g_j \in \mathcal{S}, j \le i} \frac{|s_{g_j}|^p}{N_R}, \quad P_{\mathrm{miss}}(i) = \sum_{g_j \notin \mathcal{S}, j \le i} \frac{1}{N - N_H}
$$

$$
\mathrm{ES}(\mathcal{S}) = \max_{1 \le i \le N} \left( P_{\mathrm{hit}}(i) - P_{\mathrm{miss}}(i) \right)
$$

with weighting exponent $p = 1.0$. Significance is evaluated against empirical null distributions across permutation iterations, yielding Normalized Enrichment Scores (NES) and family-wise error rates across MSigDB Hallmark, Reactome, and KEGG pathway collections.

---

## 9. Mathematical Assumptions Summary Matrix

| Modeling Domain | Mathematical Formulation | Primary Parameterization | Objective / Role in Analysis |
| :--- | :--- | :--- | :--- |
| **Count Modeling** | Negative Binomial likelihood | $\mathbf{z} \in \mathbb{R}^{96}$, dispersion $\theta_{g,b_c}$, 2 hidden layers | Deconvolves technical library variation from biological manifold coordinates |
| **Quality Control** | CellBender + empirical filter | $N_{\mathrm{genes}} \ge 200$, $\text{MT} \le 10\%$, Scrublet $\le 0.06$ | Excludes broken droplets, cell-free ambient RNA, and heterotypic doublets |
| **Batch Integration** | Conditional ELBO maximization | `sample_id` batch key + condition/cell cycle covariates | Aligns disparate sequencing batches while preserving biological states |
| **Species Alignment** | 1:1 Reciprocal orthology | Ensembl BioMart + MGI homology catalog | Harmonizes mouse and human feature spaces for joint manifold learning |
| **Potency Scoring** | Harmonic graph diffusion | Random-walk restart $\alpha = 0.90$, continuous $q_c \in [0, 1]$ | Quantifies developmental plasticity and stem-like state transitions |
| **State Transfer** | Rank-value attention encoding | 2,048 tokens, 6-layer BERT, cross-entropy loss | Maps malignant query cells to tonsil B and plasma cell compartments |
| **Trajectory Denoising** | Markov affinity diffusion | MAGIC $t=12$, $k=30$ neighbors, decay $1.0$ | Mitigates single-cell dropouts along continuous developmental paths |
| **Program Dissimilarity**| Dual-metric trajectory distance | $\alpha_{\mathrm{shape}} = 0.50$ (Pearson shape + Euclidean amplitude) | Identifies co-regulated dynamic gene programs along potency gradients |
| **Statistical Testing** | Wilcoxon rank-sum + BH-FDR | Two-sided rank test, FDR $\alpha = 0.05$ | Discovers robust, non-parametric cluster-specific biomarker panels |
| **Pathway Analysis** | Weighted Kolmogorov-Smirnov | Preranked GSEA ($p=1.0$), MSigDB Hallmark/Reactome | Identifies coordinated metabolic and oncogenic pathway activation |

---

## 10. Pipeline Module Implementation Map

| Computational Module | Executable Implementation Script |
| :--- | :--- |
| Mouse cohort preprocessing, scVI latent modeling, and CytoTRACE2 | [`mouse_scvi_cytotrace2_cellbender.py`](scripts/mouse_scvi_cytotrace2_cellbender.py) |
| Cross-species orthology mapping and human–mouse scVI integration | [`scvi_human_dlbcl_mouse_malignant_integration.py`](scripts/scvi_human_dlbcl_mouse_malignant_integration.py) |
| Harmonized joint potency prediction | [`cytotrace2_dlbcl_mouse_tonsil.py`](scripts/cytotrace2_dlbcl_mouse_tonsil.py) |
| Geneformer reference model fine-tuning on tonsil atlas | [`train_geneformer_tonsil_multi.py`](scripts/train_geneformer_tonsil_multi.py) |
| Geneformer inference and high-confidence state prediction | [`geneformer_predict_and_plot_manuscript.py`](scripts/geneformer_predict_and_plot_manuscript.py) |
| Condition-discriminative gene selection, MAGIC, and trajectory clustering | [`ordering_cytotrace_2_mouse_geneformer.py`](scripts/ordering_cytotrace_2_mouse_geneformer.py) |
| Metabolic trajectory profiling and correlation modeling | [`pgc1_cytotrace2_trajectory.py`](scripts/pgc1_cytotrace2_trajectory.py) |
| Non-parametric Wilcoxon rank-sum differential expression | [`wilcoxon_rank_mouse_integrated.py`](scripts/wilcoxon_rank_mouse_integrated.py) |

---

## 11. References

1. **CellBender:** Fleming, S. J., Chaffin, M. D., Arduini, A., Akkad, A. D., Banks, E., Marioni, J. C., Philippakis, A. A., Ellinor, P. T., & Babadi, M. (2023). Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. *Nature Methods*, 20(9), 1323–1335. https://doi.org/10.1038/s41592-023-01943-7
2. **scVI:** Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15(12), 1053–1058. https://doi.org/10.1038/s41592-018-0229-2
3. **scvi-tools:** Gayoso, A., Lopez, R., Xing, G., Boyeau, P., Valiollah Pour Amiri, V., Hong, J., Chen, W., Wu, K., Jayasuriya, M., Mehlman, E., Lange, M., Yarats, D., Regier, J., & Yosef, N. (2022). A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*, 40(2), 163–166. https://doi.org/10.1038/s41587-021-01206-w
4. **Scrublet:** Wolock, S. L., Lopez, R., & Klein, A. M. (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. *Cell Systems*, 8(4), 281–291.e9. https://doi.org/10.1016/j.cels.2018.11.005
5. **Scanpy:** Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19(1), 15. https://doi.org/10.1186/s13059-017-1382-0
6. **Leiden Algorithm:** Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9(1), 5233. https://doi.org/10.1038/s41598-019-41695-z
7. **CytoTRACE:** Gulati, G. S., Sikandar, S. S., Wesche, D. J., Manjunath, A., Bharadwaj, A., Berger, M. J., Ilagan, F., Kuo, A. H., Hurlbut, N. K., Newman, A. M., & Clarke, M. F. (2020). Single-cell transcriptional diversity is a hallmark of developmental potential. *Science*, 367(6476), 405–411. https://doi.org/10.1126/science.aax0249
8. **CytoTRACE2:** Kang, M. T., Gulati, G. S., & Newman, A. M. (2024). CytoTRACE 2: Cellular Potency and Lineage Reconstruction from Single-Cell RNA Sequencing. *Nature*, in press / bioRxiv. https://github.com/digitalcytometry/cytotrace2
9. **Geneformer:** Theodoris, C. V., Xiao, L., Chopra, A., Chaffin, M. D., Al Sayed, Z. R., Hill, M. C., Mantineo, H., Brydon, E. M., Zeng, Z., Liu, X. S., & Ellinor, P. T. (2023). Transfer learning enables predictions in network biology. *Nature*, 618(7965), 616–624. https://doi.org/10.1038/s41586-023-06139-9
10. **MAGIC:** van Dijk, D., Sharma, R., Nainys, J., Yim, K., Kathail, P., Carr, A. J., Burdziak, C., Moon, K. R., Chaffer, C. L., Pattabiraman, D., Bierie, B., Mazutis, L., Wolf, G., Krishnaswamy, S., & Pe'er, D. (2018). Recovering Gene Interactions from Single-Cell Data Using Data Diffusion. *Cell*, 174(3), 716–729.e27. https://doi.org/10.1016/j.cell.2018.05.061
11. **GSEA:** Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
12. **GSEApy:** Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a Python package for gene set enrichment analysis. *Bioinformatics*, 39(1), btac757. https://doi.org/10.1093/bioinformatics/btac757
13. **Ensembl BioMart:** Kinsella, R. J., Kähäri, A., Haider, S., Zamora, J., Proctor, G., Spudich, G., Almeida-King, J., Staines, D., Derwent, P., Kerhornou, A., Kersey, P., & Flicek, P. (2011). Ensembl BioMarts: a hub for data retrieval across taxonomic space. *Database*, 2011, bar030. https://doi.org/10.1093/database/bar030
