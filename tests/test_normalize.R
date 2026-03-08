library(testthat)
library(Seurat)
library(reticulate)
use_condaenv("samflow", required = TRUE)

source("../samflow_functions.R")

reticulate::py_run_string("
import os, warnings
os.environ['KMP_WARNINGS'] = '0'
warnings.filterwarnings('ignore')
import scanpy as sc
import scipy.sparse
")

DATA_PATH <- "../filtered_feature_bc_matrix/"

samflow_load(DATA_PATH)
samflow_normalize()

# ── R side ───────────────────────────────────────────────────────────────────

test_that("normalize: R @data slot is populated", {
  data_mat <- GetAssayData(samflow_obj, layer = "data")
  expect_true(nrow(data_mat) > 0)
  expect_true(any(data_mat > 0))
})

test_that("normalize: R @data differs from raw counts", {
  counts <- GetAssayData(samflow_obj, layer = "counts")
  data   <- GetAssayData(samflow_obj, layer = "data")
  # Compare sums — zeros are identical before/after normalization so don't index blindly
  expect_false(isTRUE(all.equal(sum(counts), sum(data))))
})

# ── Python side ───────────────────────────────────────────────────────────────

test_that("normalize: Python layers['counts'] exists with raw counts", {
  reticulate::py_run_string("
_has_counts = 'counts' in samflow_obj.layers
_counts_nnz = int(samflow_obj.layers['counts'].nnz)
")
  expect_true(reticulate::py$`_has_counts`)
  expect_gt(reticulate::py$`_counts_nnz`, 0)
})

test_that("normalize: Python X differs from raw counts after normalization", {
  reticulate::py_run_string("
import numpy as np
# Compare matrix sums — scalar indexing on sparse matrices differs across scipy versions
_x_sum   = float(samflow_obj.X.sum())
_raw_sum = float(samflow_obj.layers['counts'].sum())
_x_differs = (_x_sum != _raw_sum)
")
  expect_true(reticulate::py$`_x_differs`)
})

test_that("normalize: Python layers['lognorm'] is stashed and matches X", {
  reticulate::py_run_string("
_has_lognorm  = 'lognorm' in samflow_obj.layers
_lognorm_sum  = float(samflow_obj.layers['lognorm'].sum())
_x_sum        = float(samflow_obj.X.sum())
_lognorm_eq_x = abs(_lognorm_sum - _x_sum) < 1e-3
")
  expect_true(reticulate::py$`_has_lognorm`)
  expect_true(reticulate::py$`_lognorm_eq_x`)
})

test_that("normalize: scale factor stored in Python uns", {
  reticulate::py_run_string("
_scale_factor = samflow_obj.uns['samflow']['scale_factor']
")
  expect_equal(reticulate::py$`_scale_factor`, 10000)
})

# ── Cross-language spot check ─────────────────────────────────────────────────

test_that("normalize: R and Python log-norm values match within tolerance", {
  # Pull first 5 genes x first 3 cells from R (genes x cells)
  r_mat <- as.matrix(GetAssayData(samflow_obj, layer = "data")[1:5, 1:3])

  # Pull same values from Python (cells x genes — transposed)
  reticulate::py_run_string("
import numpy as np
_py_mat = samflow_obj.X[:3, :5]
if scipy.sparse.issparse(_py_mat):
    _py_mat = _py_mat.toarray()
import numpy as np
_py_vals = _py_mat.flatten().tolist()
")
  py_vals <- reticulate::py$`_py_vals`

  # R is genes x cells, Python is cells x genes — transpose R to match
  r_vals <- as.vector(t(r_mat))

  expect_equal(r_vals, py_vals, tolerance = 1e-5)
})
