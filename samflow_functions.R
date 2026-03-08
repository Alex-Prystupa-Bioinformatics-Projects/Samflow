library(reticulate)

SCANPY_TO_SEURAT <- c(
  "total_counts"      = "nCount_RNA",
  "n_genes_by_counts" = "nFeature_RNA"
)

# Safe rename: looks up each column individually — order-independent
.rename_scanpy_cols <- function(df) {
  new_names <- colnames(df)
  for (i in seq_along(new_names)) {
    if (new_names[i] %in% names(SCANPY_TO_SEURAT)) {
      new_names[i] <- SCANPY_TO_SEURAT[[new_names[i]]]
    }
  }
  colnames(df) <- new_names
  df
}

# Push a data.frame to Python as adata.obs via column-level manipulation
.push_obs_to_python <- function(df) {
  py$samflow_obs <- df
  reticulate::py_run_string("
_target = samflow_obs.columns.tolist()
for col in list(samflow_obj.obs.columns):
    if col not in _target:
        del samflow_obj.obs[col]
for col in _target:
    samflow_obj.obs[col] = samflow_obs[col].values
samflow_obj.obs = samflow_obj.obs[_target]
del _target
")
}

# ── Sync functions ────────────────────────────────────────────────────────────

# Union merge — used only by samflow_load()
# Both objects end up with all columns, Seurat naming as standard
.samflow_sync_meta_union <- function() {
  r_meta <- .samflow_obj_r@meta.data
  py_obs <- .rename_scanpy_cols(as.data.frame(reticulate::py$samflow_obj$obs))

  # Push Python-only columns to R
  py_only <- setdiff(colnames(py_obs), colnames(r_meta))
  if (length(py_only) > 0) {
    r_meta[py_only] <- py_obs[rownames(r_meta), py_only, drop = FALSE]
  }

  # Push R-only columns to Python
  r_only <- setdiff(colnames(r_meta), colnames(py_obs))
  if (length(r_only) > 0) {
    py_obs[r_only] <- r_meta[rownames(py_obs), r_only, drop = FALSE]
  }

  # Reorder Python obs: R column order first, Python-only appended
  py_obs <- py_obs[, colnames(r_meta), drop = FALSE]

  .samflow_obj_r@meta.data <<- r_meta
  .push_obs_to_python(py_obs)

  invisible(NULL)
}

# Directional sync — R is truth, Python is overwritten to match R
.samflow_sync_meta_r_to_py <- function() {
  .push_obs_to_python(.samflow_obj_r@meta.data)
  invisible(NULL)
}

# Directional sync — Python is truth, R is overwritten to match Python
.samflow_sync_meta_py_to_r <- function() {
  py_obs <- .rename_scanpy_cols(as.data.frame(reticulate::py$samflow_obj$obs))
  .samflow_obj_r@meta.data <<- py_obs
  .push_obs_to_python(py_obs)
  invisible(NULL)
}

# Public sync — directional, explicit user intent
samflow_sync <- function(from = "r") {
  if (from == "r") {
    .samflow_sync_meta_r_to_py()
  } else if (from == "python") {
    .samflow_sync_meta_py_to_r()
  }
  invisible(NULL)
}

samflow_load <- function(path) {
  .samflow_obj_r <<- Read10X(path) |> CreateSeuratObject()

  reticulate::py_run_string(glue::glue(
    "samflow_obj = sc.read_10x_mtx('{path}', var_names='gene_symbols')\n",
    "sc.pp.calculate_qc_metrics(samflow_obj, inplace=True)"
  ))

  .samflow_sync_meta_union()

  makeActiveBinding("samflow_obj", function(v) {
    if (missing(v)) .samflow_obj_r else .samflow_obj_r <<- v
  }, globalenv())

  invisible(NULL)
}
