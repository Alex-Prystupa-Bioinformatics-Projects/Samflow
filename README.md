# SAMFLOW

A Quarto-based framework for truly equal bidirectional R/Python access to the same single-cell object — without lossy conversion.

SAMFLOW keeps a Seurat object (R) and an AnnData object (Python) in sync within a single Quarto notebook, using reticulate as the bridge. Each language works with its native object; SAMFLOW handles the translation.

## Why

Most single-cell workflows force a choice: Seurat in R or scanpy in Python. Converting between them is lossy, one-directional, and breaks the analysis narrative. SAMFLOW lets you use Seurat for clustering, Harmony for integration, scVI for deep learning, and CellRank for trajectory — all on the same cells, in the same notebook.

## Status

v1 — metadata sync implemented and tested. Count matrix sync coming next.

### What works
- `samflow_load()` — loads PBMC-style 10x data into both Seurat and AnnData simultaneously, runs a union merge so both objects start with identical metadata
- `samflow_sync()` — R → Python directional sync (R is truth)
- `samflow_sync(from = "python")` — Python → R directional sync (Python is truth)
- Seurat naming convention enforced (`total_counts` → `nCount_RNA`, `n_genes_by_counts` → `nFeature_RNA`)
- 43 tests passing (33 R / 10 Python)

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

# Work in R with Seurat
samflow_obj@meta.data$my_cluster <- Idents(samflow_obj)
samflow_sync()  # push to Python
```

```python
# Work in Python with AnnData — same cells, same metadata
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
make test
```

## Architecture

```
samflow_functions.R
├── samflow_load(path)          # load + union merge
├── samflow_sync(from="r")      # directional sync
├── .samflow_obj_r              # internal Seurat storage
├── samflow_obj                 # active binding (public name)
├── .samflow_sync_meta_union()  # load-time bidirectional merge
├── .samflow_sync_meta_r_to_py()
└── .samflow_sync_meta_py_to_r()
```

Seurat column naming is the canonical standard. Known scanpy QC column names are remapped on sync (`SCANPY_TO_SEURAT` mapping).
