#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the CREBBP aggressive lymphoma single-cell analysis pipeline.
Tests foundational statistical, mathematical, and data-transformation routines.
"""

import unittest
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from scipy.spatial.distance import pdist, squareform


class TestSingleCellPreprocessing(unittest.TestCase):
    """Test data filtering and quality control mathematical operations."""

    def setUp(self):
        np.random.seed(42)
        self.n_cells = 50
        self.n_genes = 100
        counts = np.random.negative_binomial(n=5, p=0.7, size=(self.n_cells, self.n_genes))
        counts[counts < 2] = 0
        self.sparse_counts = sp.csr_matrix(counts)
        self.gene_names = [f"Gene_{i}" for i in range(self.n_genes)]
        self.gene_names[0] = "mt-Nd1"
        self.gene_names[1] = "mt-Nd2"
        self.gene_names[2] = "Rps3"
        self.gene_names[3] = "Rpl5"
        self.gene_names[4] = "Ighv1-1"

    def test_confounder_mask(self):
        """Verify confounder mask detects mitochondrial, ribosomal, and V(D)J genes."""
        varnames = pd.Index(self.gene_names)
        mt_mask = varnames.str.startswith("mt-") | varnames.str.startswith("MT-")
        ribo_mask = varnames.str.startswith(("Rps", "Rpl", "RPS", "RPL"))
        ig_mask = varnames.str.startswith(("Ighv", "Igkv", "Iglv", "IGHV", "IGKV", "IGLV"))
        combined_mask = mt_mask | ribo_mask | ig_mask

        self.assertTrue(combined_mask[0])   # mt-Nd1
        self.assertTrue(combined_mask[1])   # mt-Nd2
        self.assertTrue(combined_mask[2])   # Rps3
        self.assertTrue(combined_mask[3])   # Rpl5
        self.assertTrue(combined_mask[4])   # Ighv1-1
        self.assertFalse(combined_mask[5])  # Gene_5

    def test_qc_mitochondrial_fraction(self):
        """Test calculation of mitochondrial percentage."""
        mt_indices = [0, 1]
        dense_counts = self.sparse_counts.toarray()
        total_counts = np.sum(dense_counts, axis=1)
        total_counts[total_counts == 0] = 1
        mt_counts = np.sum(dense_counts[:, mt_indices], axis=1)
        pct_mt = (mt_counts / total_counts) * 100.0

        self.assertEqual(len(pct_mt), self.n_cells)
        self.assertTrue(np.all(pct_mt >= 0.0))
        self.assertTrue(np.all(pct_mt <= 100.0))


class TestTrajectoryMathematics(unittest.TestCase):
    """Test pseudotime binning, trajectory correlation, and dual-metric clustering."""

    def setUp(self):
        np.random.seed(42)
        self.n_bins = 20
        t = np.linspace(0, 1, self.n_bins)
        g0 = np.exp(-((t - 0.2) ** 2) / (2 * 0.05 ** 2))
        g1 = np.exp(-((t - 0.8) ** 2) / (2 * 0.05 ** 2))
        g2 = t
        g3 = 1 - t
        matrix = np.vstack([g0, g1, g2, g3]) + np.random.normal(0, 0.01, (4, self.n_bins))
        self.mat_log = pd.DataFrame(matrix, index=["EarlyPeak", "LatePeak", "Increasing", "Decreasing"])

    def test_spearman_correlation(self):
        """Verify Spearman correlation correctly distinguishes monotone trends."""
        pseudotime = np.linspace(0, 1, self.n_bins)
        r_inc, _ = stats.spearmanr(self.mat_log.loc["Increasing"].values, pseudotime)
        r_dec, _ = stats.spearmanr(self.mat_log.loc["Decreasing"].values, pseudotime)
        self.assertGreater(r_inc, 0.95)
        self.assertLess(r_dec, -0.95)

    def test_dual_metric_distance(self):
        """Verify combined distance balances correlation shape and expression amplitude."""
        alpha_shape = 0.5
        z_scores = stats.zscore(self.mat_log.values, axis=1)
        d_shape = pdist(z_scores, metric="correlation")
        d_amp = pdist(self.mat_log.values, metric="euclidean")
        scale = np.percentile(d_amp, 95) + 1e-12
        d_amp_s = np.clip(2.0 * (d_amp / scale), 0.0, 2.0)
        d_combined = alpha_shape * d_shape + (1.0 - alpha_shape) * d_amp_s

        self.assertEqual(len(d_combined), 6)  # 4 choose 2 pairs
        self.assertTrue(np.all(d_combined >= 0.0))
        self.assertTrue(np.all(d_combined <= 2.0))


class TestStatisticalInference(unittest.TestCase):
    """Test non-parametric hypothesis testing and FDR control."""

    def test_mann_whitney_u(self):
        group_a = np.array([12, 14, 15, 18, 19, 21, 23, 24])
        group_b = np.array([2, 3, 5, 6, 7, 8, 9, 10])
        stat, pval = stats.mannwhitneyu(group_a, group_b, alternative='two-sided')
        self.assertEqual(stat, len(group_a) * len(group_b))
        self.assertLess(pval, 0.001)

    def test_benjamini_hochberg(self):
        p_values = np.array([0.001, 0.01, 0.02, 0.04, 0.15, 0.45, 0.80])
        m = len(p_values)
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        adj_p = np.zeros(m)
        cummin = 1.0
        for i in range(m - 1, -1, -1):
            rank = i + 1
            cummin = min(cummin, sorted_p[i] * m / rank)
            adj_p[i] = cummin
        fdr = np.zeros(m)
        fdr[sorted_indices] = np.clip(adj_p, 0.0, 1.0)
        self.assertTrue(np.all(fdr >= p_values))
        self.assertLess(fdr[0], 0.01)


if __name__ == "__main__":
    unittest.main()
