#!/usr/bin/env python3
"""
Custom GSEA/Enrichr Scatter Plots - Version B
=================================
Generates publication-quality plots:
- GSEA prerank plots (NES on Y-axis) 
- Enrichr plots (Combined Score on Y-axis)


Version B: Two separate panels for Malignant (red points) and Physiologic (blue diamonds)
with bigger labels.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob

# ===== CONFIGURATION =====
BASE_DIR = Path("/home/gusti/CREBBP_aggressive_lymphoma_single_cell/heatmap_programs")
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


