# SAMFLOW Operations Reference

> **Two patterns govern all SAMFLOW operations:**
> - **Wrapper needed** — anything that modifies the expression matrix, embedding spaces, or feature sets where cross-language consistency is critical
> - **Native + `samflow_sync()`** — anything that produces metadata (labels, scores, annotations) — run in the best language, push result via `samflow_sync()`

---

## Preprocessing

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Cell/gene filtering | `subset()` | `sc.pp.filter_cells/genes()` | `samflow_filter()` wrapper | Must remove identical cells/genes from both — any mismatch breaks all downstream |
| Normalization | `NormalizeData()` | `normalize_total()` + `log1p()` | `samflow_normalize()` ✓ | Scale factor must match exactly; auto-stashes `layers["counts"]` + `layers["lognorm"]` |
| Doublet detection | `DoubletFinder` | `scrublet` / `scDblFinder` | `samflow_find_doublets()` wrapper | Run once in best language, sync score + flag as metadata |
| Ambient RNA correction | `SoupX` | `CellBender` | `samflow_filter()` after | Correction changes counts — must reload both sides |
| Mitochondrial % | Auto on load | Auto on load | Handled by `samflow_load()` ✓ | Part of initial metadata union |
| Cell cycle scoring | `CellCycleScoring()` | `sc.tl.score_genes_cell_cycle()` | Native + `samflow_sync()` | Scores + phase are just metadata |

---

## Feature Selection & Scaling

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| HVG selection | `FindVariableFeatures()` | `sc.pp.highly_variable_genes()` | `samflow_find_hvg()` ✓ | R selects, Python receives exact gene set — one selection, both languages use it |
| Scaling | `ScaleData()` | `sc.pp.scale()` | `samflow_scale()` wrapper | HVG scope + clip params must match; NOT synced — regenerated per language before PCA |

---

## Dimensionality Reduction

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| PCA | `RunPCA()` | `sc.tl.pca()` | `samflow_pca()` wrapper | Run once in R, sync embeddings to Python — avoids sign flip divergence between implementations |
| UMAP | `RunUMAP()` | `sc.tl.umap()` | `samflow_umap()` wrapper | Stochastic — must run once and sync coordinates, never independently |
| tSNE | `RunTSNE()` | `sc.tl.tsne()` | Same pattern as UMAP | Same reason |
| Diffusion map | custom | `sc.tl.diffmap()` | Python-native + `samflow_sync()` | Sync embedding to R as custom reduction |
| PHATE | `phateR` | `phate` | Python-native + `samflow_sync()` | Sync embedding coordinates |

---

## Graph Construction & Clustering

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Neighbor graph | `FindNeighbors()` | `sc.pp.neighbors()` | Run independently | Graphs NOT synced by design — different algorithms (SNN vs kNN), regenerate per language |
| Leiden/Louvain clustering | `FindClusters()` | `sc.tl.leiden()` | Native + `samflow_sync()` | Labels are metadata — run in either language, sync via obs |
| Resolution sweep | `FindClusters(resolution=)` | `sc.tl.leiden(resolution=)` | Native + `samflow_sync()` | Results table + chosen labels → metadata sync |
| Hierarchical clustering | `BuildClusterTree()` | `sc.tl.dendrogram()` | Language-native | Visualization artifact, no sync needed |

---

## Cell Type Annotation

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Manual labeling | `RenameIdents()` | `adata.obs["cell_type"] =` | `samflow_sync()` after | Labels are metadata |
| Reference-based (automated) | `Azimuth` / `SingleR` | `celltypist` / `sctype` | Native + `samflow_sync()` | R ecosystem stronger for SingleR; Python for celltypist — run in best language |
| Gene module scoring | `AddModuleScore()` | `sc.tl.score_genes()` | Native + `samflow_sync()` | Scores → obs columns |
| Marker-based scoring | `UCell` (R) | `decoupleR` (Python) | Native + `samflow_sync()` | Activity scores → metadata |

---

## Differential Expression

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Pairwise DE | `FindMarkers()` | `sc.tl.rank_genes_groups()` | Language-native | Returns results table, not a slot — no sync needed |
| One-vs-all markers | `FindAllMarkers()` | `sc.tl.rank_genes_groups()` | Language-native | Same |
| Pseudobulk DE | `DESeq2` / `edgeR` via R | `pydeseq2` | R-native preferred | R ecosystem much stronger; results are tables |
| MAST hurdle model | `MAST` (R) | — | R-native | No Python equivalent |
| Mixed-effects DE | `NEBULA` / `glmmTMB` (R) | — | R-native | No equivalent in scanpy |

---

## Trajectory & RNA Velocity

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| RNA velocity | — | `scVelo` | Python-native + `samflow_sync()` | Sync latent time + velocity embedding to R |
| Trajectory inference | `Monocle3` (R) | `CellRank` (Python) | Language-native preference | Sync pseudotime as metadata via `samflow_sync()` |
| Fate probability | — | `CellRank` | Python-native + `samflow_sync()` | Fate probs → obs columns |
| Dynamical modeling | — | `scVelo` dynamical | Python-native | Computationally intensive — Python only |

---

## Batch Correction & Integration

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Harmony | `RunHarmony()` | `harmonypy` | `samflow_pca()` pattern | Run once, sync corrected embedding both ways |
| scVI / scANVI | — | `scvi-tools` | Python-native + `samflow_sync()` | Sync latent embedding + denoised expression to R |
| BBKNN | — | `bbknn` | Python-native | Graph-level correction — not synced, used for Python clustering |
| Scanorama | — | `scanorama` | Python-native + `samflow_sync()` | Sync corrected embedding |
| Label transfer | `Azimuth` / Seurat v5 | `scArches` / `scVI` | Language-native | Predicted labels → metadata sync |
| RPCA integration | `Seurat v5` | — | R-native + `samflow_sync()` | Sync corrected embedding to Python |

---

## Gene Programs & Regulons

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| NMF / cNMF | custom | `cNMF` | Python-native + `samflow_sync()` | Sync usage scores as obs columns |
| SCENIC (TF regulons) | `SCENIC` (R) | `pySCENIC` | Python-native preferred | `pySCENIC` faster; sync TF activity scores to R |
| Gene set enrichment | `fgsea` / `clusterProfiler` (R) | `decoupleR` / `GSEApy` | R-native preferred | Rich R ecosystem; results are tables |
| Pathway scoring | `AddModuleScore()` | `sc.tl.score_genes()` | Native + `samflow_sync()` | Scores → obs |

---

## Cell-Cell Communication

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| Ligand-receptor | `CellChat` / `NicheNet` (R) | `CellPhoneDB` / `LIANA` (Python) | Language-native | Results are interaction tables, not single-cell slots — no sync needed |
| Spatial CCC | `Seurat spatial` | `squidpy` | Out of scope v1 | Requires spatial coordinates |

---

## Multi-omics (Future)

| Operation | Seurat | scanpy | SAMFLOW approach | Notes |
|---|---|---|---|---|
| WNN (RNA+ATAC) | `FindMultiModalNeighbors()` | `muon` / `MultiVI` | Out of scope v1 | Requires multi-assay support |
| Chromatin accessibility | `Signac` (R) | `snapatac2` (Python) | Out of scope v1 | — |
| CITE-seq (protein) | `Seurat v5` | `muon` | Out of scope v1 | — |
| Spatial transcriptomics | `Seurat spatial` | `squidpy` / `scanpy` | Out of scope v1 | — |

---

## Summary

| Pattern | When to use | Examples |
|---|---|---|
| **`samflow_*()` wrapper** | Cross-language consistency critical; matrix/embedding/feature modifications | `samflow_normalize`, `samflow_find_hvg`, `samflow_pca`, `samflow_umap`, `samflow_filter` |
| **Native + `samflow_sync()`** | Produces metadata only; language has clear best tool | Clustering, annotation, DE, trajectory, integration labels |
| **Language-native only** | No equivalent in other language; results are tables not slots | Pseudobulk DE, cell-cell communication, SCENIC, RNA velocity |
