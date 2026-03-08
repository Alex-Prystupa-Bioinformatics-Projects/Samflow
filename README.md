# SAMFLOW

A Quarto-based framework for truly equal bidirectional R/Python access to the same single-cell object — without lossy conversion.

SAMFLOW keeps a Seurat object (R) and an AnnData object (Python) in sync within a single Quarto notebook, using reticulate as the bridge. Each language works with its native object; SAMFLOW handles the translation.

## Why

Most single-cell workflows force a choice: Seurat in R or scanpy in Python. Converting between them is lossy, one-directional, and breaks the analysis narrative. SAMFLOW lets you use Seurat for clustering, Harmony for integration, scVI for deep learning, and CellRank for trajectory — all on the same cells, in the same notebook.

## Status

v1 in progress — preprocessing pipeline through HVG selection implemented and tested.

### What works
- `samflow_load()` — loads 10x data into both Seurat and AnnData simultaneously with union metadata merge, Seurat naming enforced
- `samflow_sync()` — R → Python directional metadata sync
- `samflow_sync(from = "python")` — Python → R directional metadata sync
- `samflow_normalize()` — log-normalizes both objects with matching scale factor (10k); auto-stashes `layers["counts"]` and `layers["lognorm"]` so raw and log-norm data are never lost
- `samflow_find_hvg()` — runs `FindVariableFeatures` in R, copies the exact HVG gene set to Python's `var["highly_variable"]`; one selection, both languages use it
- 59 tests passing (49 R / 10 Python)

### Design principles
- Operations where **native flexibility matters** (clustering, annotation, plotting) → work in either language, call `samflow_sync()` after
- Operations where **cross-language consistency is critical** (normalization, HVG selection, scaling) → dedicated `samflow_*()` wrappers run both sides together
- Data management is invisible — SAMFLOW stashes layers automatically, users never manually copy matrices

### v1 Constraints
- RNA assay only (no multi-assay)
- Single sample (no integration)
- NormalizeData only (no SCTransform)
- No spatial

## Usage

```r
# In a Quarto notebook (engine: knitr)
library(Seurat)
library(reticulate)
use_condaenv("samflow", required = TRUE)
```

```python
import scanpy as sc
```

```r
source("samflow_functions.R")

samflow_load("filtered_feature_bc_matrix/")
samflow_normalize()
samflow_find_hvg(nfeatures = 2000)
# samflow_scale()   ← coming next
# samflow_pca()     ← coming next

# Work in R with Seurat
samflow_obj@meta.data$my_cluster <- Idents(samflow_obj)
samflow_sync()  # push to Python
```

```python
# Work in Python with AnnData — same cells, same metadata, same HVGs
import scvi
scvi.model.SCVI.setup_anndata(samflow_obj)
```

## Setup

**R packages:** Seurat, reticulate, glue

**Python (conda env `samflow`):**
```bash
conda create -n samflow python=3.12
conda activate samflow
pip install scanpy anndata pytest
```

**Run tests:**
```bash
make test        # all tests
make test-unit   # fast Tier 1 only (no data needed)
make test-log    # save timestamped logs to tests/logs/
```

## Architecture

```
samflow_functions.R
├── samflow_load(path)            # load both objects + union metadata merge
├── samflow_normalize(scale=10k)  # log-norm both, stash counts + lognorm layers
├── samflow_find_hvg(nfeatures)   # R selects HVGs, copies exact set to Python
├── samflow_sync(from="r")        # directional metadata sync
├── .samflow_obj_r                # internal Seurat storage
├── samflow_obj                   # active binding (public-facing name)
└── SCANPY_TO_SEURAT              # canonical column name mapping
```

## What's next

- `samflow_scale()` — ScaleData equivalent (HVGs only, not synced — regenerated per language)
- `samflow_pca()` — RunPCA equivalent with cross-language embedding sync
- Count matrix sync (`@assays$RNA@counts` ↔ `layers["counts"]`)
- `var` / `meta.features` sync
