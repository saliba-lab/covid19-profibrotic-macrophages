################################################################################
# scRNA-seq analysis of BAL macrophages

# ------------------------------------------------------------------------------
# Download data from public repository

# Barcodes
con <- "https://nubes.helmholtz-berlin.de/s/8DYsAwEMD3CniNj/download"
file <- "~/Downloads/BAL-macrophages_barcodes.tsv"
download.file(url = con, destfile = file)
barcodes <- readLines(file)

# Features
con <- "https://nubes.helmholtz-berlin.de/s/rDAk9Sa7n8rGGG6/download"
file <- "~/Downloads/BAL-macrophages_features.tsv"
download.file(url = con, destfile = file)
features <- readr::read_tsv(file, col_names = c("ENSEMBL", "SYMBOL", "TYPE"))

# Matrix
con <- "https://nubes.helmholtz-berlin.de/s/dS39WttqQz7qedp/download"
file <- "~/Downloads/BAL-macrophages_matrix.mtx"
download.file(url = con, destfile = file)
matrix <- Matrix::readMM(gzfile(file))
colnames(matrix) <- barcodes
rownames(matrix) <- features$ENSEMBL

# ------------------------------------------------------------------------------
# Re-shape data into Seurat object

# Check features
unique(features$TYPE)
table(stringr::str_detect(features$ENSEMBL, "ENSG"))
table(stringr::str_detect(features$ENSEMBL, "SCoV2"))

scov2 <- matrix[which(stringr::str_detect(rownames(matrix), "SCoV2")), ]
matrix <- matrix[which(stringr::str_detect(rownames(matrix), "ENSG")), ]

# Compare matrix size
lapply(
  c("Gene expression"       = matrix, 
    "SARS-CoV-2 mRNAs"      = scov2
  ), 
  dim
)

# Create Seurat object and remove raw matrix
object <- Seurat::CreateSeuratObject(
  counts  = matrix, 
  assay   = "RNA",
  project = "COVID-19 BAL macrophages"
)
rm(matrix)

# Add viral counts and hashtags as assays
object@assays[["SCoV2"]] <- Seurat::CreateAssayObject(counts = scov2)
rm(scov2)

# Store features in the object
object@misc$features <- features

# Add patient metadata
metadata <- dplyr::tibble(
  dataset = c("C3", "C4", "C5", "C6", "C7", "D9", "E3", "E4"),
  patient = c("Neph_X1", "C19-62", "C19-83", "C19-82", 
              "C19-85", "C19-98", "C19-136", "C19-120"),
  dpso    = c("day 14", "day 14", "day 7", "day 14", 
              "day 25", "day 7", "day 25", "day 10"),
  age     = c("51", "72", "56", "72", "63", "22", "62", "76"),
  sex     = c("Female", "Male", "Female", "Male", 
              "Male", "Male", "Female", "Male")
)
object@meta.data <- cbind(
  object@meta.data, 
  metadata[stringr::str_split(colnames(object), "-", simplify = TRUE)[, 2], ]
  )
object@meta.data$dataset <- factor(object@meta.data$dataset)
object@meta.data$dpso <- factor(
  object@meta.data$dpso, c("day 7", "day 10", "day 14", "day 25")
  )
object@meta.data$patient <- factor(
  object@meta.data$patient,
  unique(object@meta.data$patient[order(object@meta.data$dpso)])
)
object@meta.data$age <- factor(object@meta.data$age)
object@meta.data$sex <- factor(object@meta.data$sex)

# ------------------------------------------------------------------------------
# Quality control & filtering

# Calculate proportion of mitochondrial gene counts
mt_counts <- Matrix::colSums(object@assays$RNA@counts[
  features$ENSEMBL[which(stringr::str_detect(features$SYMBOL, "^MT-"))]
  , ]) 
object@meta.data$percent.mt <- round(
  mt_counts / object@meta.data$nCount_RNA * 100, 1
)

# Set cell (column) filter
# Based on cutoffs for QC metrics
cells <- rownames(object@meta.data)[which(
  object@meta.data$nFeature_RNA > 0 &
    object@meta.data$nFeature_RNA < Inf &
    object@meta.data$nCount_RNA > 0 &
    object@meta.data$nCount_RNA < Inf &
    object@meta.data$percent.mt > 0 & 
    object@meta.data$percent.mt < Inf
)]

# Calculate proportion of cells retained
round(prop.table(table(rownames(object@meta.data) %in% cells)) * 100, 2)

# Set gene (row) filter
genes <- list()
for (assay in names(object@assays)) {
  print(assay)
  genes[[assay]] <- names(which(
    Matrix::rowSums(slot(object@assays[[assay]], "counts")) > 0
  ))
}
genes <- as.character(unlist(genes))

# Apply filter
object <- subset(object, cells = cells, features = genes)

# ------------------------------------------------------------------------------
# Normalize mRNA counts

object <- Seurat::NormalizeData(
  object               = object,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
object <- Seurat::FindVariableFeatures(
  object           = object,
  selection.method = "vst", 
  nfeatures        = 3000
)

# ------------------------------------------------------------------------------
# Calculate low-dimensional embedding & clustering

# MNN-corrected PCA
set.seed(1993)
object@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    object@assays$RNA@data[object@assays$RNA@var.features, ],
    batch = object@meta.data$patient,
    d     = 40
  )@int_colData$reducedDims$corrected,
  key        = "PC_",
  assay      = "RNA"
)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:40, seed.use = 1993)

# Clustering
set.seed(1993)
object <- Seurat::FindNeighbors(object = object, dims = 1:40)
object@meta.data$Cluster <- factor(leiden::leiden(
  object@graphs$RNA_snn, resolution_parameter = 0.9, seed = 1993
))

# ------------------------------------------------------------------------------
# Annotate clusters & select high-quality cells

# Annotate
cluster.annotation <- c(
  "1"  = "Neutrophils",
  "2"  = "Low quality",
  "3"  = "Neutrophils",
  "4"  = "T cells",
  "5"  = "T cells",
  "6"  = "Neutrophils",
  "7"  = "T cells",
  "8"  = "Macrophages",
  "9"  = "Low quality",
  "10" = "Macrophages",
  "11" = "Low quality",
  "12" = "Macrophages",
  "13" = "NK cells",
  "14" = "T cells",
  "15" = "T cells",
  "16" = "Plasma cells",
  "17" = "Plasma cells",
  "18" = "Low quality",
  "19" = "Low quality",
  "20" = "Erythrocytes",
  "21" = "Low quality",
  "22" = "Plasma cells",
  "23" = "Macrophages",
  "24" = "Macrophages",
  "25" = "B and DC",
  "26" = "Plasma cells",
  "27" = "Low quality",
  "28" = "Plasma cells",
  "29" = "Neutrophils"
)
object@meta.data$Celltype <- factor(
  cluster.annotation[as.character(object$Cluster)],
  unique(cluster.annotation)[c(4,1,8,6,3,5,7,2)]
  )

# Plot
Seurat::DimPlot(object, group.by = "Celltype") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::coord_fixed()

# Select
cells <- colnames(object)[which(
  object@meta.data$patient != "Neph_X1" &
    object@meta.data$Celltype != "Low quality"
)]
genes <- list()
for (assay in names(object@assays)) {
  print(assay)
  genes[[assay]] <- names(which(
    Matrix::rowSums(
      slot(subset(object, cells = cells)@assays[[assay]], "counts")
      ) > 0
  ))
}
genes <- as.character(unlist(genes))

# Filter
object <- subset(object, cells = cells, features = genes)

# ------------------------------------------------------------------------------
# Normalize mRNA counts

object <- Seurat::NormalizeData(
  object               = object,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
object <- Seurat::FindVariableFeatures(
  object           = object,
  selection.method = "vst", 
  nfeatures        = 3000
)

# ------------------------------------------------------------------------------
# Calculate low-dimensional embedding & clustering

# MNN-corrected PCA
set.seed(1993)
object@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    object@assays$RNA@data[object@assays$RNA@var.features, ],
    batch = object@meta.data$patient,
    d     = 40
  )@int_colData$reducedDims$corrected,
  key        = "PC_",
  assay      = "RNA"
)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:40, seed.use = 1993)

# Clustering
set.seed(1993)
object <- Seurat::FindNeighbors(object = object, dims = 1:40)
object@meta.data$Cluster <- factor(leiden::leiden(
  object@graphs$RNA_snn, resolution_parameter = 0.9, seed = 1993
))

# ------------------------------------------------------------------------------
# Show celltypes

color.celltype <- c(
  "Neutrophils"  = "darkorange2",
  "Macrophages"  = "indianred",
  "T cells"      = "slateblue",
  "NK"           = "purple",
  "Plasma"       = "lightseagreen",
  "Ery"          = "pink2",
  "Epithelium"   = "goldenrod",
  "B and DC"     = "steelblue"
)

# end of document
################################################################################