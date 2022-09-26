"
Analysis of the control dataset - quality control

Usage:
    control-analysis-qc.R [options]
    
Options:
    -h --help         Show this screen.
    --plot-dir=<path> Directory for visualizations [default: analysis/control/qc]
    --name=<file>     Name of output file [default: data/control/seurat.Rds]
" -> doc

# Global variables -------------------------------------------------------------
args <- docopt::docopt(doc)
file <- args[["--name"]]
to <- args[["--plot-dir"]]
if (!stringr::str_ends(to, "/")) {
  to <- paste0(to, "/")
}
dir.create(to, recursive = TRUE)

# Load data --------------------------------------------------------------------
ds <- readRDS(file)

# De-multiplex replicates ------------------------------------------------------

# Plot relationship
ht.min <- c(Hashtag1 = 3.4, Hashtag2 = 3.9)
df <- as.data.frame.matrix(Matrix::t(ds@assays$HTO@data))
ggplot2::ggplot(df, ggplot2::aes(Hashtag1, Hashtag2)) +
  ggplot2::geom_point(size = .1, alpha = .1) +
  ggplot2::geom_vline(xintercept = ht.min[1], col = "indianred", size = 1) +
  ggplot2::geom_hline(yintercept = ht.min[2], col = "indianred", size = 1) +
  ggplot2::theme_classic(20) + ggplot2::coord_fixed()

# Create annotation based on manual thresholds
tbl <- matrix(
  c("Negative", "A", "B", "Doublet"), nrow = 2,
  dimnames = list(c(FALSE, TRUE), c(FALSE, TRUE))
)
classify <- function(x) {
  y <- as.character(x > ht.min)
  tbl[y[1], y[2]]
}
df$Replicate <- factor(apply(df, 1, classify), as.character(tbl)[c(2,3,4,1)])

# Review annotation
ggplot2::ggplot(df, ggplot2::aes(Hashtag1, Hashtag2, col = Replicate)) +
  ggplot2::geom_point(size = 2, shape = 21, alpha = .1) +
  ggplot2::geom_vline(xintercept = ht.min[1], col = "indianred", size = 1) +
  ggplot2::geom_hline(yintercept = ht.min[2], col = "indianred", size = 1) +
  ggplot2::theme_classic(20) +
  ggplot2::guides(color = ggplot2::guide_legend(
    override.aes = list(size = 5, shape = 20)
  )) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::labs(col = NULL, x = "Hashtag 1", y = "Hashtag 2") +
  ggplot2::guides(color = ggplot2::guide_legend(
    override.aes = list(size = 5, alpha = 1, shape = 20)
  )) +
  ggplot2::theme(legend.position = c(0.8, 0.8))
filename <- paste0(to, "hashtag-distribution.png")
ggplot2::ggsave(filename, height = 6, width = 6)

# Add replicate information to metadata
ds$Replicate <- df$Replicate

# Quality control --------------------------------------------------------------

# Calculate quality clusters by kmeans on quality metrics
df <- log1p(data.frame(
  libsize    = ds$nCount_RNA,
  features   = ds$nFeature_RNA,
  percent.mt = ds$percent.mt
))
set.seed(42)
ds$qc.group <- factor(kmeans(df, 9)$cluster)

# Plot quality metrics
dat <- ds@meta.data[, c(
  "nFeature_RNA", "nCount_RNA", "percent.mt", "Replicate", "qc.group"
)]
dat <- tidyr::gather(dat, "key", "value", -nCount_RNA, -qc.group, -Replicate)
ggplot2::ggplot(
  dat, ggplot2::aes(value, log10(nCount_RNA), col = qc.group)
) +
  ggplot2::geom_point(shape = 1, size = .5) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::facet_grid(~key, scales = "free") +
  ggplot2::guides(
    col = ggplot2::guide_legend(override.aes = list(
      size = 8, shape = 20
    ))
  ) +
  ggplot2::labs(x = NULL, y = "log10(Library size)")
filename <- paste0(to, "quality.png")
ggplot2::ggsave(filename, width = 12, height = 6, bg = "white")

# Compute embedding & clusters
# Normalization & Feature Selection
ds <- Seurat::NormalizeData(ds)
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 5000)
# PCA
hvg <- ds@assays$RNA@var.features
set.seed(1993)
pca <- irlba::irlba(ds@assays$RNA@data[hvg, ], nv = 25)
ds@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  key = "PC_", embeddings = pca$v, loadings = pca$u
)
# UMAP
set.seed(1993)
ds@reductions[["umap"]] <- Seurat::CreateDimReducObject(
  key        = "UMAP_",
  embeddings = uwot::umap(
    X = ds@reductions$pca@cell.embeddings,
    metric = "cosine", n_neighbors = 20
  )
)
# Create SNN-graph for clustering
set.seed(1993)
ds@graphs$snn <- bluster::makeSNNGraph(ds@reductions$pca@cell.embeddings)

# --- Show quality groups ---
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds@meta.data$qc.group
)
ggplot2::ggplot(
  data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = c(0.2,0.8)) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 3, override.aes = list(size = 5))
  ) +
  ggplot2::labs(col = "QC group")

filename <- paste0(to, "umap-qc-group.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Show number of detected genes per cell ---
data$col <- ds$nFeature_RNA
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.2,0.8)
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      direction = "vertical", barwidth = 1, frame.colour = "black",
      ticks = FALSE, title.position = "top"
    )
  ) +
  ggplot2::labs(col = "Genes") +
  viridis::scale_color_viridis(option = "C", direction = -1)
filename <- paste0(to, "umap-nGene.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Show percentage of mitochondrial counts ---
data$col <- ds$percent.mt
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.2,0.8)
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      direction = "vertical", barwidth = 1, frame.colour = "black",
      ticks = FALSE
    )
  ) +
  ggplot2::labs(col = "% Mt-RNAs") +
  viridis::scale_color_viridis(option = "C", direction = -1)
filename <- paste0(to, "umap-pMT.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Show replicates ---
data$col <- ds$Replicate
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = c(0.2,0.8)) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, override.aes = list(
    size = 5
  ))
  ) +
  ggplot2::labs(col = "Replicate") +
  ggplot2::scale_color_brewer(palette = "Set1")
filename <- paste0(to, "umap-hashtag.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Compute and show quality clusters ---
ds$cluster <- factor(leiden::leiden(
  object = igraph::as_adjacency_matrix(ds@graphs$snn, sparse = FALSE),
  resolution_parameter = 0.2, seed = 1993,
  partition_type = "RBConfigurationVertexPartition"
))

data$col <- ds$cluster
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = c(0.2,0.8)) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, override.aes = list(
    size = 5
  ))
  ) +
  ggplot2::labs(col = "Cluster") +
  ggplot2::scale_color_brewer(palette = "Set2")
filename <- paste0(to, "umap-qc-cluster.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Remove low quality clusters ---
quality <- c("1"  = "High", "2"  = "Low", "3"  = "Low")
ds$quality <- quality[ds$cluster]
data$col <- ds$quality
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = c(0.2,0.8)) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, override.aes = list(
    size = 5
  ))
  ) +
  ggplot2::labs(col = "Quality") +
  ggplot2::scale_color_brewer(palette = "Set1")
filename <- paste0(to, "umap-quality.png")
ggplot2::ggsave(filename, width = 8, height = 6, bg = "white")

# --- Subset Seurat object ---
cells <- names(which(
  ds$quality == "High" & ds$Replicate %in% c("A", "B")
))
ds <- subset(ds, cells = cells)

# Re-level replicates
ds$Replicate <- factor(ds$Replicate)

# Save dataset -----------------------------------------------------------------
message(paste("Writing", file))
saveRDS(ds, file)