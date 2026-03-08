library(testthat)
library(Seurat)
library(reticulate)
library(dplyr)
use_condaenv("samflow", required = TRUE)

source("../samflow_functions.R")

# Suppress scanpy output during tests
reticulate::py_run_string("
import os, warnings
os.environ['KMP_WARNINGS'] = '0'
warnings.filterwarnings('ignore')
import scanpy as sc
")

DATA_PATH <- "../filtered_feature_bc_matrix/"

# Helper: get Python obs as R data.frame
# reticulate auto-converts pandas DataFrame to R data.frame on access
py_obs <- function() reticulate::py$samflow_obj$obs

# ── After samflow_load() ──────────────────────────────────────────────────────

samflow_load(DATA_PATH)

test_that("after load: R and Python have identical column names in identical order", {
  expect_equal(colnames(samflow_obj@meta.data), colnames(py_obs()))
})

test_that("after load: nCount_RNA values match between R and Python", {
  expect_equal(samflow_obj@meta.data$nCount_RNA, py_obs()[["nCount_RNA"]])
})

test_that("after load: nFeature_RNA values match between R and Python", {
  expect_equal(samflow_obj@meta.data$nFeature_RNA, py_obs()[["nFeature_RNA"]])
})

test_that("after load: cell barcodes are identical in R and Python", {
  expect_equal(rownames(samflow_obj@meta.data), rownames(py_obs()))
})

test_that("after load: scanpy QC column names follow Seurat convention", {
  r_cols <- colnames(samflow_obj@meta.data)
  expect_false("total_counts" %in% r_cols)
  expect_false("n_genes_by_counts" %in% r_cols)
  expect_true("nCount_RNA" %in% r_cols)
  expect_true("nFeature_RNA" %in% r_cols)
})

# ── samflow_sync() — R → Python ───────────────────────────────────────────────

test_that("samflow_sync(): dropping R column removes it from Python", {
  # <<- required: test_that runs in a local env, <<- ensures the active binding fires
  samflow_obj@meta.data <<- samflow_obj@meta.data %>% select(orig.ident)
  samflow_sync()
  expect_equal(colnames(py_obs()), c("orig.ident"))
})

test_that("samflow_sync(): adding R column pushes it to Python with correct values", {
  samflow_obj@meta.data$test_col <<- "hello"
  samflow_sync()
  expect_true(all(py_obs()[["test_col"]] == "hello"))
})

# ── samflow_sync(from = "python") — Python → R ───────────────────────────────

test_that("samflow_sync(from='python'): new Python column appears in R", {
  reticulate::py_run_string("samflow_obj.obs['py_test_col'] = 'world'")
  samflow_sync(from = "python")
  expect_true("py_test_col" %in% colnames(samflow_obj@meta.data))
  expect_true(all(samflow_obj@meta.data$py_test_col == "world"))
})

test_that("samflow_sync(from='python'): scanpy-named col arrives in R with Seurat naming", {
  reticulate::py_run_string("samflow_obj.obs['total_counts'] = 999.0")
  samflow_sync(from = "python")
  expect_true("nCount_RNA" %in% colnames(samflow_obj@meta.data))
  expect_false("total_counts" %in% colnames(samflow_obj@meta.data))
})
