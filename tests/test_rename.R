library(testthat)

# The rename mapping (mirrors SCANPY_TO_SEURAT in samflow_functions.R)
SCANPY_TO_SEURAT <- c(
  "total_counts"      = "nCount_RNA",
  "n_genes_by_counts" = "nFeature_RNA"
)

# Pure helper — apply rename logic to a data.frame (no reticulate needed)
rename_scanpy_cols <- function(df, mapping) {
  known <- intersect(names(mapping), colnames(df))
  if (length(known) > 0) {
    colnames(df)[colnames(df) %in% known] <- mapping[known]
  }
  df
}

test_that("total_counts is renamed to nCount_RNA", {
  df <- data.frame(total_counts = 1:3)
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_true("nCount_RNA" %in% colnames(result))
  expect_false("total_counts" %in% colnames(result))
})

test_that("n_genes_by_counts is renamed to nFeature_RNA", {
  df <- data.frame(n_genes_by_counts = 1:3)
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_true("nFeature_RNA" %in% colnames(result))
  expect_false("n_genes_by_counts" %in% colnames(result))
})

test_that("unknown columns are left unchanged", {
  df <- data.frame(log1p_total_counts = 1:3, pct_counts_top_50 = 4:6)
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_true("log1p_total_counts" %in% colnames(result))
  expect_true("pct_counts_top_50" %in% colnames(result))
})

test_that("already-Seurat-named columns are not double-renamed", {
  df <- data.frame(nCount_RNA = 1:3, nFeature_RNA = 4:6)
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_equal(colnames(result), c("nCount_RNA", "nFeature_RNA"))
})

test_that("values are preserved after rename", {
  df <- data.frame(total_counts = c(100, 200, 300),
                   n_genes_by_counts = c(50, 80, 120))
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_equal(result$nCount_RNA, c(100, 200, 300))
  expect_equal(result$nFeature_RNA, c(50, 80, 120))
})

test_that("both canonical columns renamed in one pass", {
  df <- data.frame(total_counts = 1:3, n_genes_by_counts = 4:6,
                   log1p_total_counts = 7:9)
  result <- rename_scanpy_cols(df, SCANPY_TO_SEURAT)
  expect_setequal(colnames(result), c("nCount_RNA", "nFeature_RNA", "log1p_total_counts"))
})
