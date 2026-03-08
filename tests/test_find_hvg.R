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
")

DATA_PATH <- "../filtered_feature_bc_matrix/"

samflow_load(DATA_PATH)
samflow_normalize()
samflow_find_hvg(nfeatures = 2000)

# ── R side ────────────────────────────────────────────────────────────────────

test_that("find_hvg: R has correct number of variable features", {
  expect_equal(length(VariableFeatures(samflow_obj)), 2000)
})

test_that("find_hvg: R variable features are gene names (strings)", {
  expect_type(VariableFeatures(samflow_obj), "character")
})

# ── Python side ───────────────────────────────────────────────────────────────

test_that("find_hvg: Python var has highly_variable column", {
  reticulate::py_run_string("_has_hv = 'highly_variable' in samflow_obj.var.columns")
  expect_true(reticulate::py$`_has_hv`)
})

test_that("find_hvg: Python highly_variable is boolean", {
  reticulate::py_run_string("_is_bool = samflow_obj.var['highly_variable'].dtype == bool")
  expect_true(reticulate::py$`_is_bool`)
})

test_that("find_hvg: Python HVG count matches R", {
  reticulate::py_run_string("_n_hvg = int(samflow_obj.var['highly_variable'].sum())")
  expect_equal(reticulate::py$`_n_hvg`, length(VariableFeatures(samflow_obj)))
})

# ── Cross-language: same genes selected ──────────────────────────────────────

test_that("find_hvg: R and Python have identical HVG gene sets", {
  r_hvg <- sort(VariableFeatures(samflow_obj))
  reticulate::py_run_string("
_py_hvg = sorted(samflow_obj.var_names[samflow_obj.var['highly_variable']].tolist())
")
  py_hvg <- sort(unlist(reticulate::py$`_py_hvg`))
  expect_equal(r_hvg, py_hvg)
})
