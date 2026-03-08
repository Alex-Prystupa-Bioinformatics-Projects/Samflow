library(testthat)

# Pure helper — union merge of two data.frames (no reticulate needed)
# Mirrors the logic in .samflow_sync_meta_union()
union_merge <- function(r_meta, py_obs) {
  py_only <- setdiff(colnames(py_obs), colnames(r_meta))
  if (length(py_only) > 0) {
    r_meta[py_only] <- py_obs[rownames(r_meta), py_only, drop = FALSE]
  }

  r_only <- setdiff(colnames(r_meta), colnames(py_obs))
  if (length(r_only) > 0) {
    py_obs[r_only] <- r_meta[rownames(py_obs), r_only, drop = FALSE]
  }

  py_obs <- py_obs[, colnames(r_meta), drop = FALSE]

  list(r_meta = r_meta, py_obs = py_obs)
}

# Shared fixtures
make_r_meta <- function() {
  data.frame(
    orig.ident   = c("s1", "s1", "s1"),
    nCount_RNA   = c(100L, 200L, 300L),
    nFeature_RNA = c(50L, 80L, 120L),
    row.names    = c("Cell1", "Cell2", "Cell3")
  )
}

make_py_obs <- function() {
  data.frame(
    nCount_RNA             = c(100, 200, 300),
    nFeature_RNA           = c(50, 80, 120),
    log1p_total_counts     = c(4.6, 5.3, 5.7),
    pct_counts_in_top_50   = c(30.0, 35.0, 40.0),
    row.names              = c("Cell1", "Cell2", "Cell3")
  )
}

test_that("Python-only columns are added to R", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_true("log1p_total_counts" %in% colnames(result$r_meta))
  expect_true("pct_counts_in_top_50" %in% colnames(result$r_meta))
})

test_that("R-only columns are added to Python", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_true("orig.ident" %in% colnames(result$py_obs))
})

test_that("shared columns are not duplicated", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_equal(sum(colnames(result$r_meta) == "nCount_RNA"), 1)
  expect_equal(sum(colnames(result$py_obs) == "nCount_RNA"), 1)
})

test_that("R and Python end up with identical column names", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_setequal(colnames(result$r_meta), colnames(result$py_obs))
})

test_that("Python obs column order: R cols first, Python-only appended", {
  result <- union_merge(make_r_meta(), make_py_obs())
  py_cols <- colnames(result$py_obs)
  r_cols  <- colnames(make_r_meta())
  expect_equal(py_cols[seq_along(r_cols)], r_cols)
})

test_that("cell values are correctly transferred from R to Python", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_equal(result$py_obs[["orig.ident"]], c("s1", "s1", "s1"))
})

test_that("cell values are correctly transferred from Python to R", {
  result <- union_merge(make_r_meta(), make_py_obs())
  expect_equal(result$r_meta[["log1p_total_counts"]], c(4.6, 5.3, 5.7))
})
