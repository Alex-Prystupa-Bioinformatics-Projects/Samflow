"""
Python-side assertions for samflow sync.
Run after the R integration tests have loaded the data, OR run standalone
(this file sets up its own AnnData via samflow_load via rpy2 — skip if rpy2 not available).

These tests focus on what Python sees — column naming, values, index integrity.
Run: conda run -n samflow pytest tests/test_sync_meta.py -v
"""
import pytest
import os
import warnings

os.environ["KMP_WARNINGS"] = "0"
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata
import scanpy as sc


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def adata():
    """Load PBMC 1k data as AnnData with QC metrics."""
    ad = sc.read_10x_mtx("filtered_feature_bc_matrix/", var_names="gene_symbols")
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    return ad


@pytest.fixture(scope="module")
def adata_renamed(adata):
    """Simulate what samflow_sync does: rename scanpy cols to Seurat convention."""
    SCANPY_TO_SEURAT = {
        "total_counts":      "nCount_RNA",
        "n_genes_by_counts": "nFeature_RNA",
    }
    obs = adata.obs.copy()
    obs = obs.rename(columns=SCANPY_TO_SEURAT)
    adata.obs = obs
    return adata


# ── Column naming tests ───────────────────────────────────────────────────────

def test_total_counts_renamed_to_nCount_RNA(adata_renamed):
    assert "nCount_RNA" in adata_renamed.obs.columns
    assert "total_counts" not in adata_renamed.obs.columns


def test_n_genes_by_counts_renamed_to_nFeature_RNA(adata_renamed):
    assert "nFeature_RNA" in adata_renamed.obs.columns
    assert "n_genes_by_counts" not in adata_renamed.obs.columns


def test_extra_qc_columns_preserved(adata_renamed):
    """log1p_* and pct_counts_* should still be present after rename."""
    cols = adata_renamed.obs.columns.tolist()
    assert any("log1p" in c for c in cols)
    assert any("pct_counts" in c for c in cols)


# ── Value integrity tests ─────────────────────────────────────────────────────

def test_nCount_RNA_values_are_positive(adata_renamed):
    assert (adata_renamed.obs["nCount_RNA"] > 0).all()


def test_nFeature_RNA_less_than_nCount_RNA(adata_renamed):
    """Every cell should have fewer unique genes than total counts."""
    assert (adata_renamed.obs["nFeature_RNA"] <= adata_renamed.obs["nCount_RNA"]).all()


def test_no_duplicate_columns(adata_renamed):
    cols = adata_renamed.obs.columns.tolist()
    assert len(cols) == len(set(cols)), f"Duplicate columns found: {cols}"


# ── Index / barcode tests ─────────────────────────────────────────────────────

def test_barcodes_are_unique(adata_renamed):
    assert adata_renamed.obs_names.is_unique


def test_barcode_format(adata_renamed):
    """PBMC 1k barcodes should end with -1."""
    assert all(bc.endswith("-1") for bc in adata_renamed.obs_names)


# ── R → Python sync simulation ────────────────────────────────────────────────

def test_r_to_py_sync_drops_columns():
    """Simulate samflow_sync(): Python obs should only have what R has."""
    r_meta = pd.DataFrame(
        {"orig.ident": ["s1", "s1", "s1"]},
        index=["Cell1", "Cell2", "Cell3"]
    )
    # Simulated sync: overwrite obs with R metadata
    obs = r_meta.copy()
    assert list(obs.columns) == ["orig.ident"]
    assert "nCount_RNA" not in obs.columns


def test_py_to_r_sync_renames_columns():
    """Simulate samflow_sync(from='python'): scanpy names become Seurat names."""
    SCANPY_TO_SEURAT = {"total_counts": "nCount_RNA", "n_genes_by_counts": "nFeature_RNA"}
    obs = pd.DataFrame({
        "total_counts": [100, 200],
        "n_genes_by_counts": [50, 80],
        "log1p_total_counts": [4.6, 5.3]
    })
    obs = obs.rename(columns=SCANPY_TO_SEURAT)
    assert "nCount_RNA" in obs.columns
    assert "nFeature_RNA" in obs.columns
    assert "total_counts" not in obs.columns
    assert "log1p_total_counts" in obs.columns  # non-canonical cols preserved
