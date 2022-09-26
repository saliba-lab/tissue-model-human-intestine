"
Analysis of the control dataset - cell types

Usage:
    control-analysis-celltypes.R [options]
    
Options:
    -h --help         Show this screen.
    --plot-dir=<path> Directory for visualizations [default: analysis/control/celltypes]
    --name=<file>     Name of output file [default: data/control/seurat.Rds]
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

# Global variables -------------------------------------------------------------

# Input arguments
args <- docopt::docopt(doc)
file <- args[["--name"]]
seed <- 42

# Colors
cols <- list(
  Replicate = c("A" = "darkorange", "B" = "steelblue"),
  Phase = c("G1" = "grey", "G2M" = "seagreen", "S" = "purple")
)

# Output arguments
to <- args[["--plot-dir"]]
if (!stringr::str_ends(to, "/")) {
  to <- paste0(to, "/")
}
dir.create(to, recursive = TRUE)

# Load data --------------------------------------------------------------------
ds <- readRDS(file)

# Show cell metadata on embeddings ---------------------------------------------

# --- Replicates ---
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$Replicate
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.5) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = c(0.25, 0.88)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    nrow = 1, override.aes = list(size = 5)
  )) +
  ggplot2::labs(col = "Replicate") +
  ggplot2::scale_color_manual(values = cols$Replicate)

filename <- paste0(to, "01_umap-replicate.png")
ggplot2::ggsave(filename, width = 4.5, height = 6, bg = "white")

# --- Expressed features ---
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$nFeature_RNA
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col/1000)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = "left", legend.margin = ggplot2::margin(0, -140)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 23, barwidth = 1, ticks = FALSE)
  ) +
  ggplot2::labs(col = "Genes (thsd)") +
  viridis::scale_color_viridis(option = "A", direction = -1)

filename <- paste0(to, "S01_umap-nfeatures.png")
ggplot2::ggsave(filename, width = 5, height = 6, bg = "white")

# --- Library size ---
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$nCount_RNA
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col/1000)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = "left", legend.margin = ggplot2::margin(0, -145)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 23, barwidth = 1, ticks = FALSE)
  ) +
  ggplot2::labs(col = "Counts (thsd)") +
  viridis::scale_color_viridis(option = "A", direction = -1)

filename <- paste0(to, "S02_umap-libsize.png")
ggplot2::ggsave(filename, width = 5, height = 6, bg = "white")

# --- Mitochondrial percentage ---
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$percent.mt
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = "left", legend.margin = ggplot2::margin(0, -160)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 23, barwidth = 1, ticks = FALSE)
  ) +
  ggplot2::labs(col = "Mt. counts (%)") +
  viridis::scale_color_viridis(option = "A", direction = -1)

filename <- paste0(to, "S03_umap-percent-mt.png")
ggplot2::ggsave(filename, width = 5, height = 6, bg = "white")

# Annotate cell cycle ----------------------------------------------------------

n2i <- ds@misc$features$id
names(n2i) <- ds@misc$features$name

genes <- Seurat::cc.genes.updated.2019

genes$s.genes <- n2i[genes$s.genes]
genes$g2m.genes <- n2i[genes$g2m.genes]

cells_rankings <- AUCell::AUCell_buildRankings(
  ds@assays$RNA@data, plotStats=TRUE
)
cells_AUC <- AUCell::AUCell_calcAUC(genes, cells_rankings)
data <- as.data.frame(t(cells_AUC@assays@data$AUC))
colnames(data) <- stringr::str_remove(colnames(data), ".genes")
set.seed(seed)
data$group <- factor(kmeans(data, 3)$cluster)

ggplot2::ggplot(data, ggplot2::aes(g2m, s, col = group)) +
  ggplot2::geom_point()
ds$phase <- factor(c("1" = "G1", "2" = "G2M", "3" = "S")[data$group])

data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$phase
)
ann <- dplyr::summarise(
  dplyr::group_by(data, col), x = median(x), y = median(y)
)
ann$x <- rep(-0.65, 3)
ann$y <- c(0.8,-1.2,-0.8)

data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$phase
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::geom_text(data = ann, ggplot2::aes(label = col), size = 10) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(legend.position = "") +
  ggplot2::guides() +
  ggplot2::labs(subtitle = "Cell cycle phase") +
  ggplot2::scale_color_manual(values = cols$Phase)

filename <- paste0(to, "S04_umap-phase.png")
ggplot2::ggsave(filename, width = 5, height = 6, bg = "white")

# Clustering -------------------------------------------------------------------
set.seed(seed)
ds@graphs$snn <- bluster::makeSNNGraph(
  x = ds@reductions$pca@cell.embeddings, k = 20, type = "rank"
)
ds$cluster <- factor(
  igraph::cluster_louvain(ds@graphs$snn, resolution = 1.3)$membership
)

data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$cluster
)
ann <- dplyr::summarise(
  dplyr::group_by(data, col), x = median(x), y = median(y)
)
ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::geom_label(data = ann, ggplot2::aes(label = col), size = 7) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = "",
    plot.subtitle = ggplot2::element_text(
      margin = ggplot2::margin(b = -60, t = 30)
    )
  ) +
  ggplot2::labs(col = NULL, subtitle = "Clusters")
filename <- paste0(to, "02_umap-clusters.png")
ggplot2::ggsave(filename, width = 4.5, height = 6, bg = "white")

# Annotate clusters ------------------------------------------------------------

c2l <- c(
  "1" = "Progenitor",
  "2" = "Stem",
  "3" = "Secretory",
  "4" = "Progenitor",
  "5" = "Enterocytes",
  "6" = "M-like",
  "7" = "Enterocytes",
  "8" = "Dead/empty",
  "9" = "HLA-G+",
  "10" = "Dead/empty"
)

ds$type <- c2l[ds$cluster]
ds$type[ds$phase %in% c("S", "G2M")] <- "Prolif."

cols$type <- c(
  "Prolif." = "dodgerblue",
  "Stem" = "cyan3",
  "Progenitor" = "slategrey",
  "Enterocytes" = "darkorange",
  "M-like" = "purple",
  "Secretory"   = "seagreen",
  "HLA-G+" = "deeppink",
  "Dead/empty" = "grey"
)
ds$type <- factor(ds$type, names(cols$type))

data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2],
  col = ds$type
)
ann <- dplyr::summarise(
  dplyr::group_by(data, col), x = median(x), y = median(y)
)

ggplot2::ggplot(data, ggplot2::aes(x, y, col = col)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::scale_color_manual(values = cols$type) +
  ggplot2::geom_label(data = ann, ggplot2::aes(label = col), size = 7) +
  ggplot2::coord_fixed(clip = "off") +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = "",
    plot.subtitle = ggplot2::element_text(
      margin = ggplot2::margin(b = -60, t = 30)
    )
  ) +
  ggplot2::labs(col = NULL, subtitle = "Celltype")

filename <- paste0(to, "03_umap-type.png")
ggplot2::ggsave(filename, width = 5, height = 6, bg = "white")

filename <- paste0(to, "03_umap-type.eps")
ggplot2::ggsave(filename, width = 5, height = 6)

# --- Barplot of celltype proportions ---

# Focus on high quality cells
cells <- colnames(ds)[which(ds$type != "Dead/empty")]

# Extract data & summarize proportions
prop <- data.frame(type = ds$type[cells])
prop$freq <- 1
prop <- dplyr::group_by(prop, type)
prop <- dplyr::summarise(prop, freq = sum(freq))
prop <- dplyr::mutate(prop, prop = round(freq / sum(freq) * 100), 1)

# Re-level celltype factor
prop$Type <- factor(prop$type, levels = levels(ds$type))

# Plot
ggplot2::ggplot(
  data = prop,
  mapping = ggplot2::aes(
    x = type,
    y = prop,
    fill = type,
    label = freq
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(nudge_y = 2, size = 5) +
  ggplot2::scale_fill_manual(values = cols$type) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::labs(y = "Proportion of cells in %") +
  ggplot2::theme(
    legend.position = "",
    axis.title.x    = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    panel.grid.major.y = ggplot2::element_line(size = .5, color = "grey80"),
    panel.grid.minor.y = ggplot2::element_line(size = .5)
  )

fn <- paste0(to, "04_barplot-type.png")
ggplot2::ggsave(fn, height = 6, width = 6)

# --- Barplot of celltype proportions by replicate ---

# Extract data & summarize proportions
prop <- data.frame(
  type = ds$type,
  rep = ds$Replicate
)
prop$freq <- 1
prop <- dplyr::group_by(prop, rep, type)
prop <- dplyr::summarise(prop, freq = sum(freq))
prop <- dplyr::mutate(prop, prop = round(freq / sum(freq) * 100), 1)

# Re-level celltype factor
prop$Type <- factor(prop$type, levels = levels(ds$type))

# Plot
ggplot2::ggplot(
  data = prop,
  mapping = ggplot2::aes(
    x = type,
    y = prop,
    fill = type,
    label = freq
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(nudge_y = 2, size = 5) +
  ggplot2::facet_wrap(~rep) +
  ggplot2::scale_fill_manual(values = cols$type) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::labs(y = "Proportion of cells in %") +
  ggplot2::theme(
    legend.position = "",
    axis.title.x    = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    panel.grid.major.y = ggplot2::element_line(size = .5, color = "grey80"),
    panel.grid.minor.y = ggplot2::element_line(size = .5)
  )

filename <- paste0(to, "04_barplot-type_replicate.png")
ggplot2::ggsave(filename, height = 6, width = 8)

# Differential expression across groups ----------------------------------------
mf <- model.frame(~ type + Replicate, ds@meta.data)
markers <- scran::findMarkers(
  ds@assays$RNA@data, mf[,1], block = mf[,2], pval.type = "some"
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$group <- i
  markers[[i]]$id <- rownames(markers[[i]])
  markers[[i]]$name <- names(n2i[match(markers[[i]]$id, n2i)])
}
markers <- dplyr::bind_rows(as.list(markers))
markers$group <- factor(markers$group, levels(mf[,1]))

# Overview by volcano plot
markers <- markers[markers$FDR != 1, ]
p <- ggplot2::ggplot(
  markers, ggplot2::aes(summary.logFC, -log10(FDR), col = group, label = name)
) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~group) +
  ggplot2::theme_classic(15)
p <- plotly::ggplotly(p)
fn <- paste0(to, "volcano-group", ".html")
htmlwidgets::saveWidget(widget = p, file = fn)

# Select genes by FDR
deg <- markers[markers$FDR < 1e-15 & markers$summary.logFC > 0, ]
deg <- deg[!duplicated(deg$name), ]
dim(deg)

# Limit proliferation signal to 20 genes
index <- which(deg$group == "Prolif.")
deg <- deg[-tail(index, length(index) - 40), ]

mtx <- ds@assays$RNA@data[deg$id, do.call(order, mf)]
mtx <- t(scale(Matrix::t(mtx)))
rownames(mtx) <- deg$name

# Set colorscale limits
limits <- c(-2, 2)
breaks <- seq(from = min(limits), to = max(limits), length.out = 100)

p <- pheatmap::pheatmap(
  mtx, breaks = breaks, annotation_col = mf,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_colnames = FALSE, show_rownames = FALSE,
  annotation_colors = cols,
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
  )(length(breaks))
)
fn <- paste0(to, "heatmap-cluster_scran", ".png")
ggplot2::ggsave(fn, p, width = 10, height = 5)

# Summarized
cann <- unique_combinations(mf)
groups <- factor(
  apply(mf, 1, function(x) paste0(x, collapse = "-")),
  apply(cann, 1, function(x) paste0(x, collapse = "-"))
)
mtx <- summarize_groups(ds@assays$RNA@data[deg$id, ], groups)
rownames(mtx) <- deg$name

mtx <- t(scale(t(mtx)))

p <- pheatmap::pheatmap(
  mtx, breaks = breaks, annotation_col = cann, annotation_colors = cols,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_colnames = FALSE, show_rownames = FALSE,
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
  )(length(breaks))
)
fn <- paste0(to, "heatmap-type_scran_summarized", ".png")
ggplot2::ggsave(fn, p, width = 10, height = 5)

# Export table
markers$DE <- rownames(markers) %in% rownames(deg)
index <- c("id", "name", "group", "DE", "p.value", "FDR", "summary.logFC")
index <- c(
  names(markers)[match(index, names(markers))],
  names(markers)[!names(markers) %in% index]
)
markers <- markers[, index]
fn <- paste0(to, "markers-type", ".csv")
write.table(markers, fn, sep = ",", row.names = FALSE)

# Dotplot of cell type markers -------------------------------------------------

# Select genes
genes <- list(
  "Prolif."     = c("MKI67", "PCNA", "MCM6", "BRCA2"),
  "Stem"        = c("LGR5", "SOX4", "OLFM4", "GDF15"),
  "M-like"      = c("BIRC3", "SDC4", "CCL20", "LAMC2", "NFKBIA"),
  "Progenitor"  = c("CCND2", "SYPL1", "SLC19A3"),
  "Enterocytes" = c("ZG16", "MS4A10", "SLC4A7", "RBP2"),
  "Secretory"   = c("LEMD1", "TNFSF15", "SLC5A3", "MRPS6"),
  "HLA-G+"      = c("HLA-G", "APOL4", "TRIP6", "SULT1C2"),
  "Low quality" = c("MTRNR2L12")
)

genes <- list(
  "immature" = c("ALPI", "TMEM37", "FABP2", "SI", "ANPEP"),
  "mature" = c("RBP2", "CYP3A4", "APOA4", "APOB"),
  "other" = c("ATOH1", "HES1", "DLL1", "DLL4", "BEST4", "TNFRSF11A"),
  "epi" = c("KRT8", "EPCAM", "CDH1")
)

# Select genes
genes <- list(
  "Prolif."     = c("MKI67", "PCNA", "NUSAP1"),
  "Stem"        = c("LGR5", "OLFM4", "SOX4"),
  "Progenitor"  = c(),
  "Enterocytes" = c("ANPEP", "FABP2", "APOA4", "RBP2", "APOC3", "CYP3A4"),
  "M-like"      = c("CXCL3", "CXCL8", "CCL20", "TNFAIP2"),
  "Secretory"   = c("LEMD1", "TFF1", "RAB3B", "MMP1", "DHRS9", "F3"),
  "HLA-G+"      = c("HLA-G", "APOL4", "TRIP6"),
  "Dead/empty" = c("MTRNR2L12", "KRT8", "EPCAM", "CDH1")
)

# Tidy data & compute summaries
df <- data.frame(
  group = ds$type
)
for (gene in unlist(genes)) {
  id <- n2i[gene]
  df[[gene]] <- ds@assays$RNA@data[id, ]
}
df <- tidyr::gather(df, "gene", "count", -group)
df$gene <- factor(df$gene, unique(df$gene))
df$N <- 1
df$expr <- df$count > 0
df <- dplyr::group_by(df, group, gene)
df <- dplyr::summarise(
  df, count = mean(count), N = sum(N), expr = sum(expr)
)
df$pct <- round(df$expr / df$N, 3) * 100

# Scale data & set limits
df <- dplyr::group_by(df, gene)
df <- dplyr::mutate(
  df, zscore = scale(count)[,1]
)
limits <- c(-2, 2)
df$zscore[df$zscore > max(limits)] <- max(limits)
df$zscore[df$zscore < min(limits)] <- min(limits)

# Re-order group
df$group <- factor(as.character(df$group), names(genes))

# Plot
ggplot2::ggplot(df, ggplot2::aes(gene, group, col = zscore, size = pct)) +
  ggplot2::geom_point() +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, face = "italic"
    )
  ) +
  ggplot2::scale_color_distiller(palette = "RdBu") +
  ggplot2::labs(x = NULL, y = NULL, col = "z-score", size = "% expr.") +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 15, barheight = .75, frame.colour = "black", ticks = FALSE,
      title.position = "top", title.hjust = 0.5
    ),
    size = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  )

filename <- paste0(to, "05_dotplot-type.png")
ggplot2::ggsave(filename, height = 6, width = 12)

# UMAP of celltype markers -----------------------------------------------------

# Get data
data <- tidyr::tibble(
  x = ds@reductions$umap@cell.embeddings[, 1],
  y = ds@reductions$umap@cell.embeddings[, 2]
)
for (i in unlist(genes)) {
  j <- n2i[i]
  data[[i]] <- ds[["RNA"]]@data[j, ]
}
data <- tidyr::gather(data, "gene", "norm", -x, -y)
data <- dplyr::mutate(dplyr::group_by(data, gene), zscore = scale(norm)[, 1])

# Re-order genes
data$gene <- factor(data$gene, unique(data$gene))

# Trim values
limits <- c(-2, 2)
data$zscore[data$zscore > max(limits)] <- max(limits)
data$zscore[data$zscore < min(limits)] <- min(limits)

# Plot
ggplot2::ggplot(
  data, ggplot2::aes(x, y, col = zscore)
) +
  ggplot2::geom_point(size = .1) +
  ggplot2::facet_wrap(~gene, nrow = 2) +
  ggplot2::scale_color_distiller(palette = "RdBu") +
  ggplot2::theme_void(20) +
  ggplot2::coord_fixed() +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 10, barwidth = 1, ticks = FALSE
    )
  ) +
  ggplot2::labs(color = NULL)

filename <- paste0(to, "05_umap-markers.png")
ggplot2::ggsave(filename, height = 5, width = 25, bg = "white")

# Dotplot of cytokines ---------------------------------------------------------

genes <- list(
  "IL1A", "IL1B", "IL1RN", "IL6", "IL7", "CXCL8", "IL15", "IL17C", "IL17D",
  "IL18", "IL23A", "IL32", "IL33", "IL37", "TNF"
)

# Tidy data & compute summaries
df <- data.frame(
  group = ds$type
)
for (gene in unlist(genes)) {
  id <- n2i[gene]
  df[[gene]] <- ds@assays$RNA@data[id, ]
}
df <- tidyr::gather(df, "gene", "count", -group)
df$gene <- factor(df$gene, unique(df$gene))
df$N <- 1
df$expr <- df$count > 0
df <- dplyr::group_by(df, group, gene)
df <- dplyr::summarise(
  df, count = mean(count), N = sum(N), expr = sum(expr)
)
df$pct <- round(df$expr / df$N, 3) * 100

# Scale data & set limits
df <- dplyr::group_by(df, gene)
df <- dplyr::mutate(
  df, zscore = scale(count)[,1]
)
limits <- c(-2, 2)
df$zscore[df$zscore > max(limits)] <- max(limits)
df$zscore[df$zscore < min(limits)] <- min(limits)

# Plot
ggplot2::ggplot(df, ggplot2::aes(gene, group, col = zscore, size = pct)) +
  ggplot2::geom_point() +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, face = "italic"
    ),
    legend.text = ggplot2::element_text(size = 15),
    legend.title = ggplot2::element_text(size = 15)
  ) +
  ggplot2::scale_color_distiller(palette = "RdBu") +
  ggplot2::scale_size_area(max_size = 8) +
  ggplot2::labs(x = NULL, y = NULL, col = "z-score", size = "% expr.") +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 15, barheight = .75, frame.colour = "black", ticks = FALSE,
      title.position = "top", title.hjust = 0.5
    ),
    size = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  )

filename <- paste0(to, "05_dotplot-cytokine.png")
ggplot2::ggsave(filename, height = 5, width = 10)

# Gene set overlaps ------------------------------------------------------------

# Fetch reference gene sets
ref <- "tables/Elmentaite2021_signatures.xls"
dict <- readxl::read_excel(ref)[, c("genes", "cluster")]
index <- unique(dict$cluster)[1:22]
dict <- dict[which(dict$cluster %in% index), ]
dict <- dict[dict$genes %in% ds@misc$features$name, ]
dict$id <- n2i[dict$genes]

# Fetch control gene sets
de <- markers[markers$FDR < 1e-15 & markers$summary.logFC > 0, ]

# Compare overlap
gs.ctrl <- split(de$name, de$group)
gs.ref <- split(dict$genes, dict$cluster)

gs.comp <- matrix(
  ncol = 4, nrow = (length(names(gs.ctrl)) * length(names(gs.ref)))
)
colnames(gs.comp) <- c("from", "to", "overlap", "total")
gs.comp <- as.data.frame(gs.comp)
gs.genes <- list()
n <- 1
for (i in names(gs.ctrl)) {
  for (j in names(gs.ref)) {
    
    key <- paste0(i, "_vs_", j)
    ctrl <- gs.ctrl[[i]]
    ref <- gs.ref[[j]]
    
    gs.comp[n, "from"] <- i
    gs.comp[n, "to"] <- j
    gs.comp[n, "overlap"] <- length(which(ctrl %in% ref))
    gs.comp[n, "total"] <- length(c(ctrl, ref))
    gs.genes[[key]] <- ctrl[which(ctrl %in% ref)]
    
    n <- n+1
    
  }
}

# Remove non-overlapping comparisons
index <- which(gs.comp$overlap > 0)
gs.comp <- gs.comp[index, ]
gs.genes <- gs.genes[index]

# Plot
ggplot2::ggplot(
  gs.comp, ggplot2::aes(from, overlap/total, fill = to, label = overlap)
) +
  ggplot2::geom_col() +
  ggplot2::geom_label() +
  ggplot2::facet_wrap(~ to) +
  ggplot2::theme_classic(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
  )

filename <- paste0(to, "barplot_geneset-overlap.png")
ggplot2::ggsave(filename, height = 6, width = 12)

# Assess correlation to in vivo reference cell types ---------------------------

ref <- "tables/Elmentaite2021_signatures.xls"
dict <- readxl::read_excel(ref)[, c("genes", "cluster")]
index <- unique(dict$cluster)[1:22]
dict <- dict[which(dict$cluster %in% index), ]
dict <- dict[dict$genes %in% ds@misc$features$name, ]
dict$id <- n2i[dict$genes]

de <- markers[markers$FDR < 1e-15 & markers$summary.logFC > 0, ]
genes <- split(de$id, de$group)
genesets <- split(dict$id, dict$cluster)

# Contingency table
matrix(
  c("A", "B", "C", "D"), nrow = 2,
  dimnames = list(
    c("DE", "Not.DE"), c("In.gene.set", "Not.in.gene.set"))
)

# Calculate fisher test for each cluster-gene set pair
result <- list()
for (i in names(genes)) {
  if (length(genes[[i]]) > 0) {
    
    print(i)
    de <- as.character(genes[[i]])
    bg <- as.character(markers$gene[markers$cluster == i])
    bg <- bg[which(!bg %in% de)]
    
    result[[i]] <- data.frame(
      Geneset = names(genesets), Cluster = i, p.value = NA, overlap = NA
    )
    
    for (ii in names(genesets)) {
      print(ii)
      
      A <- length(de[de %in% genesets[[ii]]])
      B <- length(bg[bg %in% genesets[[ii]]])
      
      C <- length(de[!de %in% genesets[[ii]]])
      D <- length(bg[!bg %in% genesets[[ii]]])
      
      deTable <- matrix(
        c(A, B, C, D), nrow = 2,
        dimnames = list(
          c("DE", "Not.DE"), c("In.gene.set", "Not.in.gene.set"))
      )
      
      result[[i]]$p.value[result[[i]]$Geneset == ii] <- fisher.test(
        deTable, alternative = "greater"
      )$p.value
      
      result[[i]]$overlap[result[[i]]$Geneset == ii] <- round(
        A / (A + B + C), 3
      ) * 100
      
    }
    # End of ii
    
    result[[i]]$p.adjust <- p.adjust(result[[i]]$p.value, method = "BH")
    
  }
  
}
# End of i

data <- dplyr::bind_rows(result)
data$Cluster <- factor(data$Cluster, levels(ds$type))

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = Geneset,
    y = Cluster,
    fill = overlap
  )
) +
  ggplot2::geom_tile(col = "white", size = 1) +
  ggplot2::scale_fill_gradient2(low = "white", high = "red") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    )
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      barwidth = 39, barheight = 1, ticks = FALSE, frame.colour = "black",
      title.position = "top", title.hjust = 0.5
    )
  )
fn <- paste0(to, "geneset-overlap", ".eps")
ggplot2::ggsave(fn, width = 10, height = 6, bg = "white")

# Calculate activity of gene sets with AUCell  ---------------------------------

# Build gene expression rankings
cells_rankings <- AUCell::AUCell_buildRankings(
  ds@assays$RNA@data, plotStats=TRUE
)

# Calculate enrichment for the gene signatures (AUC)
de <- markers[markers$FDR < 1e-15 & markers$summary.logFC > 0, ]
modules <- split(de$id, de$group)
modules$`Dead/empty` <- NULL
cells_AUC <- AUCell::AUCell_calcAUC(modules, cells_rankings)
df <- as.data.frame(t(cells_AUC@assays@data$AUC))

# Extract AUC and turn into long format
df$celltype <- ds$type
df$x <- ds@reductions$umap@cell.embeddings[, 1]
df$y <- ds@reductions$umap@cell.embeddings[, 2]

# Wrangle data into long format
df <- tidyr::gather(df, "module", "score", -celltype, -x, -y)

# Re-order
df <- df[order(df$score, decreasing = TRUE), ]
df$module <- factor(df$module, names(modules))

df <- dplyr::mutate(
  dplyr::group_by(df, module), scaled = scale(score)[,1]
)

limits <- c(-2, 2)
df$scaled[df$scaled > max(limits)] <- max(limits)
df$scaled[df$scaled < min(limits)] <- min(limits)

df <- df[order(df$scaled), ]

# Plot
ggplot2::ggplot(df, ggplot2::aes(x, y, col = scaled)) +
  ggplot2::geom_point(size = .1, alpha = .25) +
  ggplot2::scale_color_distiller(palette = "RdBu") +
  ggplot2::facet_wrap(~module, nrow = 1) +
  ggplot2::theme_void(15) +
  ggplot2::theme(legend.position = "right") +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 8, barwidth = 1,
      frame.colour = "black", ticks = FALSE, title = NULL
    )
  ) +
  ggplot2::coord_fixed()
fn <- paste0(to, "umap-type-AUC", ".png")
ggplot2::ggsave(fn, width = 12, height = 2.5, bg = "white")

ggplot2::ggplot(
  df, ggplot2::aes(celltype, score, col = scaled, fill = celltype)
) +
  ggplot2::geom_point(position = "jitter") +
  ggplot2::scale_color_distiller(palette = "RdBu") +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::scale_fill_manual(values = cols$type) +
  ggplot2::facet_wrap(~module, scales = "free", nrow = 1) +
  ggplot2::theme_classic(10) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::guides(
    col = ggplot2::guide_none(), fill = ggplot2::guide_none()
  )
fn <- paste0(to, "violin-type-AUC.png")
ggplot2::ggsave(fn, height = 3, width = 12, bg = "white")

# Enrichment of gene ontology (GO) terms ---------------------------------------

# Change ENSEMBL IDs to ENTREZ
input <- list()
for (i in names(modules)) {
  input[[i]] <- modules[[i]]
  input[[i]] <- clusterProfiler::bitr(
    geneID   = input[[i]],
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )[[2]]
}

# Calculate GO enrichment
result <- clusterProfiler::compareCluster(
  geneCluster = input,
  fun         = "enrichGO",
  OrgDb       = "org.Hs.eg.db",
  ont         = "ALL"
)

# Extract data
data <- result@compareClusterResult

# Add ratio of between overlap and input
gr <- as.data.frame(
  stringr::str_split(data$GeneRatio, "/", simplify = TRUE)
)
gr$V1 <- as.numeric(as.character(gr$V1))
gr$V2 <- as.numeric(as.character(gr$V2))
data$GeneRatio <- gr$V1 / gr$V2

# Select top categories per condition
data <- dplyr::group_by(data, Cluster)
categories <- dplyr::top_n(data, 7, -p.adjust)$Description
length(categories)
data <- data[data$Description %in% categories, ]

# Sort GO terms
data$Description <- factor(
  x = data$Description,
  levels = rev(unique(categories))
)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Description,
    col  = p.adjust,
    size = GeneRatio
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_distiller(palette = "YlOrRd", direction = -1) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.title  = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12
    )
  ) +
  ggplot2::scale_size_area(max_size = 5) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      reverse = TRUE, barheight = 10, barwidth = 1, ticks = FALSE
    )
  )

# Save plot & table
fn <- paste0(to, "dotplot-marker-go", ".png")
ggplot2::ggsave(fn, width = 9, height = 7)

write.table(
  x         = result@compareClusterResult,
  file      = stringr::str_replace(fn, "png", "csv"),
  sep       = ",",
  row.names = FALSE
)

# ChIP-seq enrichment analysis (ChEA3) based on DE genes -----------------------

genes <- modules

# ChEA3 query for each cluster
result <- list()
for (i in names(genes)) {
  print(i)
  # Skip empty categories
  if(length(genes[[i]]) > 0) {
    # Prepare input as list
    gs <- names(n2i[match(genes[[i]], n2i)])
    payload = list(query_name = "myQuery", gene_set = gs)
    
    # POST to ChEA3 server
    response <- httr::POST(
      url    = "https://amp.pharm.mssm.edu/chea3/api/enrich/",
      body   = payload,
      encode = "json"
    )
    json <- httr::content(response, "text")
    
    # Convert JSON to data.frame and store in list
    result[[i]] <- jsonlite::fromJSON(json)
  }
}

# Select integrated result (meanRank) and add cluster column
for (j in names(result)) {
  print(j)
  result[[j]] <- result[[j]]$`Integrated--meanRank`
  result[[j]][["cluster"]] <- j
  result[[j]][["bg"]] <- length(genes[[j]])
}
result <- dplyr::bind_rows(result)

# Adjust columns / add metrics
result$Rank <- as.numeric(result$Rank)
result$Score <- as.numeric(result$Score)
result$map <- unlist(
  lapply(stringr::str_split(result$Overlapping_Genes, ","), length)
)
result$geneRatio <- round(result$map / result$bg, 3)

# Order result by Score
result$cluster <- factor(result$cluster, names(genes))
result <- result[order(result$cluster, result$Score), ]

# Select TFs by mean rank
tfs <- result$TF[result$Score < 35]
data <- result[which(result$TF %in% tfs), ]
data$TF <- factor(data$TF, levels = unique(tfs))
data$cluster <- factor(data$cluster, names(genes))

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = TF,
    y = cluster,
    fill = Score
  )
) +
  ggplot2::geom_tile(col = "white") +
  viridis::scale_fill_viridis(
    option = "D", direction = 1, trans = "log10"
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    )
  ) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Mean Rank") +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      barheight = 7, barwidth = 1, frame.colour = "black", ticks = FALSE,
      reverse = TRUE
    )
  )

# Save plot
fn <- paste0(to, "heatmap-marker-TFs.png")
ggplot2::ggsave(fn, width = 12, height = 3)
write.csv(
  result[, c(7,3,2,4,8,9,10,5,6,1)],
  stringr::str_replace(fn, "png", "csv"), row.names = FALSE
)

# Save dataset -----------------------------------------------------------------
message(paste("Writing", file))
saveRDS(ds, file)
