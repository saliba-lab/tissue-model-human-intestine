"
Analysis of the control dataset - embedding & clustering

Usage:
    control-analysis-core.R [options]
    
Options:
    -h --help         Show this screen.
    --name=<file>     Name of output file [default: data/control/seurat.Rds]
" -> doc

# Global variables -------------------------------------------------------------

# Input arguments
args <- docopt::docopt(doc)
file <- args[["--name"]]
seed <- 42

# Load data --------------------------------------------------------------------
ds <- readRDS(file)

# Compute embedding & clusters -------------------------------------------------

ds <- Seurat::NormalizeData(ds)
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 5000)

# PCA
hvg <- ds@assays$RNA@var.features
set.seed(seed)
pca <- irlba::irlba(
  Matrix::t(
    scale(Matrix::t(ds@assays$RNA@data[hvg, ]), center = TRUE, scale = TRUE)
  ),
  nv = 35
)
ds[["pca"]] <- Seurat::CreateDimReducObject(
  key = "PC_", embeddings = pca$v, loadings = pca$u
)
rownames(ds[["pca"]]@cell.embeddings) <- colnames(ds)
rownames(ds[["pca"]]@feature.loadings) <- hvg

# UMAP
set.seed(seed)
umap <- uwot::umap(
  X = ds@reductions$pca@cell.embeddings,
  init = ds@reductions$pca@cell.embeddings[, 1:2], local_connectivity = 3,
  metric = "cosine", n_neighbors = 100, min_dist = 0, spread = 0.22,
  ret_model = TRUE
)
ds[["umap"]] <- Seurat::CreateDimReducObject(
  key        = "UMAP_",
  embeddings = umap$embedding, misc = list(model = umap)
)
rownames(ds[["umap"]]@cell.embeddings) <- colnames(ds)
ds[["umap"]]@cell.embeddings <- ds[["umap"]]@cell.embeddings[, 2:1]

# Create SNN-graph for clustering
set.seed(seed)
ds@graphs$snn <- bluster::makeSNNGraph(
  x = ds@reductions$pca@cell.embeddings, k = 10, type = "rank"
)

# Compute clusters
ds$cluster <- factor(leiden::leiden(
  object = igraph::as_adjacency_matrix(ds@graphs$snn, sparse = FALSE),
  resolution_parameter = 0.8, seed = seed, n_iterations = 5
))

ds$cluster <- factor(leiden::leiden(
  object = ds@graphs$snn,
  resolution_parameter = 0.8, seed = seed, n_iterations = 5
))

# Save dataset -----------------------------------------------------------------
message(paste("Writing", file))
saveRDS(ds, file)
