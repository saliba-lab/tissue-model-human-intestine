"
Download the control dataset

Usage:
    dataset-control.R [options]
Options:
    -h --help         Show this screen.
    --name=<file>     Name of output file [default: data/control/seurat.Rds]
" -> doc

# Global variables -------------------------------------------------------------
args <- docopt::docopt(doc)
outfile <- args[["--name"]]
dir.create(dirname(outfile), recursive = TRUE)

# Download data ----------------------------------------------------------------
url <- "https://nubes.helmholtz-berlin.de/s/EabnT8WSrR4Bmqn/download"
dir <- tempdir()
file <- paste0(dir, "/matrix.h5")
download.file(url, file)

# Extract data components ------------------------------------------------------
rhdf5::h5ls(file)
var <- dplyr::bind_rows(
  rhdf5::h5read(file, "matrix/features")[c("id", "name", "feature_type")]
)
ds <- Matrix::sparseMatrix(
  i = rhdf5::h5read(file, "matrix/indices")+1,
  p = rhdf5::h5read(file, "matrix/indptr"),
  x = rhdf5::h5read(file, "matrix/data"),
  dimnames = list(
    var$id,
    rhdf5::h5read(file, "matrix/barcodes")
  )
)

# Split the data into transcripts and hashtags
hto <- ds[split(var$id, var$feature_type)[["Antibody Capture"]], ]
ds <- ds[split(var$id, var$feature_type)[["Gene Expression"]], ]

# Create seurat object ---------------------------------------------------------
ds <- Seurat::CreateSeuratObject(
  counts = ds, project = "hITM_control", assay = "RNA"
)
ds@misc$features <- var

# Compute percentage of mitochondrial counts
index <- match(var$id[grepl("^MT-", var$name)], rownames(ds))
ds$percent.mt <- round(
  Matrix::colSums(ds@assays$RNA@counts[index, ]) / ds$nCount_RNA, 3
) * 100

# Add hashtag counts & transform
ds[["HTO"]] <- Seurat::CreateAssayObject(counts = hto) ; rm(hto)
ds@assays$HTO@data <- log1p(ds@assays$HTO@counts)

# Clean-up ---------------------------------------------------------------------
file.remove(file)

# Save dataset -----------------------------------------------------------------
saveRDS(ds, outfile)