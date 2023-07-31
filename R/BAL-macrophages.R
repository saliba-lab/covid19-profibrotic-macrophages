# Libraries --------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
library(batchelor)

# Variables --------------------------------------------------------------------
url <- list(
  "counts"="https://syncandshare.desy.de/index.php/s/pZ7aWcNYjMnWC2E",
  "metadata"="https://syncandshare.desy.de/index.php/s/ZqzjMZxLgNDcQSR"
)
url <- lapply(url, paste0, "/download")

data_path <- "data/bal/"
dir.create(data_path, recursive = TRUE)

files <- list(
  "counts"="data/bal/counts.h5",
  "metadata"="data/bal/metadata.xlsx"
)

plot_dir <- "results/bal/"
dir.create(plot_dir, recursive = TRUE)

# R datasets
out_bal <- "data/BAL.Rds"
out_bal_macrophages <- "data/BAL-macrophages.Rds"

# Colors -----------------------------------------------------------------------

# Timepoint
color.timepoint <- c(
  "day 7"  = "darkorange",
  "day 10" = "indianred",
  "day 14" = "purple",
  "day 25" = "navy",
  "" = "black"
)
object@misc$colors$Timepoint <- color.timepoint

# Celltype
color.celltype <- c(
  "Neutrophils"      = "darkorange2",
  "Macrophages"      = "indianred",
  "T cells"          = "deepskyblue",
  "NK cells"         = "purple",
  "Plasma cells"     = "lightseagreen",
  "Erythrocytes"     = "pink2",
  "Epithelial cells" = "goldenrod",
  "B and DC"         = "steelblue"
)
object@misc$colors <- list(Celltype = color.celltype)

# Macrophage subtypes
color.cluster <- c(
  "Monocytes"     = "#009E73",
  "Mono/Mφ"       = "#E69F00",
  "CD163/LGMN-Mφ" = "#D55E00",
  "AMφ-1"         = "#A020F0",
  "AMφ-2"         = "#0072B2",
  "Prolif. AMφ"   = "#56B4E9"
)
object@misc$colors$Subtype <- color.cluster

# Download data ----------------------------------------------------------------
options(timeout = Inf)

for (name in names(files)) {
  file <- files[[name]]
  if (file.exists(file)) {
    print(paste("File", file, "already exists and will be skipped..."))
    next
  } else {
    print(paste("Downloading", name))
    download.file(url[[name]], file)
  }
}

# Read data --------------------------------------------------------------------
counts <- Read10X_h5(files$counts)
metadata <- readxl::read_excel(files$metadata)

# Create Seurat object ---------------------------------------------------------

# Viral genes
index <- which(stringr::str_detect(rownames(counts), "SCoV2"))
viral <- counts[index, ]

# Human genes
index <- which(!stringr::str_detect(rownames(counts),"SCoV2"))
counts <- counts[index, ]

# Create Seurat object
object <- Seurat::CreateSeuratObject(
  counts  = counts, 
  assay   = "RNA",
  project = "COVID-19 BAL macrophages"
)
rm(counts) # remove raw data

# Add viral counts as assay
object[["SCoV2"]] <- Seurat::CreateAssayObject(counts = viral)
rm(viral)

# Add patient metadata
object@meta.data <- cbind(
  object@meta.data, 
  metadata[stringr::str_split(colnames(object), "-", simplify = TRUE)[, 2], ]
  )
object$dataset <- factor(object@meta.data$dataset)
object$dpso <- paste("day", object$dpso)
object$dpso <- factor(
  object$dpso, c("day 7", "day 10", "day 14", "day 25")
  )
object$patient <- factor(
  object$patient,
  unique(object$patient[order(object@meta.data$dpso)])
)
object$age <- factor(object$age)
object$sex <- factor(object$sex)

# Quality control --------------------------------------------------------------

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

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

# Normalization ----------------------------------------------------------------

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

# Integration ------------------------------------------------------------------

# MNN-corrected PCA
set.seed(1993)
object@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    object@assays$RNA@data[object@assays$RNA@var.features, ],
    batch = object$patient,
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
object <- Seurat::FindClusters(object, algorithm = 3, resolution = 0.9)
object$Cluster <- object$seurat_clusters

# Annotate cell types ----------------------------------------------------------

# Plot
Seurat::DimPlot(object, group.by = "Cluster", label = TRUE) +
  ggplot2::coord_fixed()

# Annotate
c2l <- list(
  "Neutrophils" = c(0,2,5),
  "Macrophages" = c(6,11,13,25),
  "T cells" = c(12,14,4,3,7),
  "NK cells" = c(17),
  "Plasma cells" = c(15,27,22,23,19),
  "Erythrocytes" = c(20),
  "Epithelial cells" = c(21),
  "B and DC" = c(24)
)
v <- c()
for (i in names(lapply(c2l, length))) {v <- c(v, rep(i, length(c2l[[i]])))}
names(v) <- unlist(c2l)

object$Celltype <- "Low quality"
index <- object$Cluster %in% names(v)
object$Celltype[index] <- v[as.character(object$Cluster[index])]
object$Celltype <- factor(object$Celltype, names(c2l))

# Plot
Seurat::DimPlot(object, group.by = "Celltype") +
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

# Remove low quality cells -----------------------------------------------------
object <- subset(object, cells = cells, features = genes)

# Normalization ----------------------------------------------------------------

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

# Embedding --------------------------------------------------------------------

# Re-level factor
object$patient <- factor(object$patient)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:40, seed.use = 1993)

# Plot: Celltypes --------------------------------------------------------------

# Re-level factor (remove 'Low quality')
object$Celltype <- factor(object$Celltype)

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = object@meta.data$Celltype
)

# Set annotation coordinates
ann <- data.frame(
  col = levels(data$col),
  x   = c(-3, -2, 9,5,  5,  4,-1,-7),
  y   = c( 8,-14,-3,8,-11,-17,-2,-5)
)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col,
    label = col
  )
) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::geom_text(
    data = ann, size = 8
  ) +
  ggplot2::scale_color_manual(values = color.celltype) +
  ggplot2::theme_void(base_size = 20) +
  ggplot2::theme(
    legend.position = ""
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL_Celltypes.png")
ggplot2::ggsave(fn, width = 5, height = 5.5, bg="white")

# Create DE table between cell types
markers <- scran::findMarkers(
  x         = object@assays$RNA@data,
  groups    = object@meta.data$Celltype,
  pval.type = "some",
  test.type = "wilcox",
  direction = "up",
  block     = object@meta.data$patient
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$cluster <- i
  markers[[i]]$gene <- row.names(markers[[i]])
}
markers <- dplyr::bind_rows(as.list(markers))
markers <- markers[, c(11,12,1,2,3:10,13)]
fn <- paste0(plot_dir, "BAL_Celltype-markers.csv")
write.csv(markers, fn, row.names = FALSE)

# Plot: Patients ---------------------------------------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = object@meta.data$patient
)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col
  )
) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::theme_void(base_size = 25) +
  ggplot2::theme(
    legend.position = "right"
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::labs(title = "Patient", col = NULL) +
  ggplot2::guides(
    col = ggplot2::guide_legend(
      override.aes = list(size = 10)
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL_Patients.png")
ggplot2::ggsave(fn, width = 7, height = 6, bg="white")

# Create DE table between patients
markers <- scran::findMarkers(
  x         = object@assays$RNA@data,
  groups    = object@meta.data$patient,
  pval.type = "some",
  test.type = "wilcox",
  direction = "up"
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$cluster <- i
  markers[[i]]$gene <- row.names(markers[[i]])
}
markers <- dplyr::bind_rows(as.list(markers))
markers <- markers[, c(10,11,1,2,3:9,12)]
fn <- paste0(plot_dir, "BAL_Patient-markers.csv")
write.csv(markers, fn, row.names = FALSE)

# Plot: Timepoints (DPSO) ------------------------------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = object@meta.data$dpso
)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col,
    shape = col
  )
) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::scale_color_manual(values = color.timepoint) +
  ggplot2::theme_void(base_size = 25) +
  ggplot2::theme(
    legend.position = "right"
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::labs(title = "Timepoint", shape = NULL, color = NULL) +
  ggplot2::guides(
    col = ggplot2::guide_legend(
      override.aes = list(size = 10), keyheight = grid::unit(12, "mm")
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL_Timepoints.png")
ggplot2::ggsave(fn, width = 7, height = 6, bg="white")

# Create DE table between timepoints
markers <- scran::findMarkers(
  x         = object@assays$RNA@data,
  groups    = object@meta.data$dpso,
  pval.type = "some",
  test.type = "wilcox",
  direction = "up"
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$cluster <- i
  markers[[i]]$gene <- row.names(markers[[i]])
}
markers <- dplyr::bind_rows(as.list(markers))
markers <- markers[, c(7,8,1,2,3:6,9)]
fn <- paste0(plot_dir, "BAL_Timepoint-markers.csv")
write.csv(markers, fn, row.names = FALSE)

# Plot: SCoV2 counts -----------------------------------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = Matrix::colSums(object@assays$SCoV2@counts)
)
data <- data[order(data$col), ]
data$col[data$col == 0] <- NA # Hide zero-values

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col
  )
) +
  ggplot2::geom_point(size = 0.5) +
  viridis::scale_color_viridis(
    option = "A", direction = -1, trans = "log10"
  ) +
  ggplot2::theme_void(base_size = 25) +
  ggplot2::theme(
    legend.position = "right"
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::labs(color = NULL, title = "SARS-CoV-2 mRNAs") +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barheight = 20, barwidth = 1, frame.colour = "black", ticks = FALSE
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL_viral-mRNA.png")
ggplot2::ggsave(fn, width = 6, height = 6, bg="white")

# Dotplot: Cell type markers ---------------------------------------------------

# Select marker genes
genes <- list(
  "Macrophages"      = c("CD68", "CD163"),
  "Neutrophils"      = c("CSF3R", "CXCL8", "FCGR3B"),
  "NK cells"         = c("NKG7", "KLRC1", "KLRB1"),
  "T cells"          = c("CD3D", "CD8A", "CD8B", "CD4"),
  "B and DC"         = c("IRF8", "MS4A1", "CD19", "HLA-DRB5", "CLEC9A"),
  "Plasma cells"     = c("MZB1", "JCHAIN", "IGHM"),
  "Epithelial cells" = c("KRT18", "KRT8", "EPCAM"),
  "Erythrocytes"     = c("HBA1", "HBA2", "HBB")
)

genes <- unlist(genes[levels(object$Celltype)])
data <- tidyr::as_tibble(t(as.matrix(object@assays$RNA@data[genes, ])))

# Add metadata
data$Celltype <- object$Celltype

# Tidy data & compute summaries
data <- tidyr::gather(data, "Gene", "Expression", -Celltype)
data$Gene <- factor(data$Gene, levels = unique(data$Gene))
data$N <- rep(1, length(data$Celltype))
data$Expressed <- data$Expression > 0
data <- dplyr::group_by(data, Celltype, Gene)
data <- dplyr::summarize(
  .data        = data,
  "Expression" = mean(Expression),
  "Cells"      = sum(N),
  "Expressed"  = sum(Expressed)
)
data$PCT <- data$Expressed / data$Cells * 100

# Scale data & set limits
data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(
  data, zscore = scale(Expression)
)
limits <- c(-2, 2)
data$zscore[data$zscore > max(limits)] <- max(limits)
data$zscore[data$zscore < min(limits)] <- min(limits)

# Plot
ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x    = Gene,
    y    = Celltype,
    col  = zscore,
    size = PCT
  )
) +
  ggplot2::geom_point() +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 18, barheight = 1,
      ticks = FALSE, frame.colour = "black",
      title.position = "top", title.hjust = 0.5
    ),
    size  = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  ) +
  ggplot2::labs(
    col  = "z-score", size = "% expressed"
  ) +
  ggplot2::theme_classic(base_size = 25) +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x     = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, face = "italic"
    ),
    axis.title      = ggplot2::element_blank()
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdBu", limits = limits
  ) +
  ggplot2::scale_size_area(max_size = 10)

# Save plot
fn <- paste0(plot_dir, "BAL_Celltype-markers_dotplot.png")
ggplot2::ggsave(fn, width = 12, height = 6)

# Save DE table of marker genes
fn_in <- paste0(plot_dir, "BAL_Celltype-markers.csv")
fn_out <- paste0(plot_dir, "BAL_Celltype-markers_dotplot.csv")
markers <- read.csv(fn_in)
write.csv(markers[markers$gene %in% unlist(genes), ], fn_out)

# Barplot: Cell type x Patient -------------------------------------------------

# Fetch metadata
data <- tidyr::tibble(
  Patient   = object$patient,
  Celltype  = object$Celltype,
  Timepoint = object$dpso,
  Cells     = 1
)

# Count cells by patient and celltype
data <- dplyr::summarise(
  dplyr::group_by(data, Patient, Celltype, Timepoint), 
  Cells = sum(Cells)
)

# Count cells by patient
data <- dplyr::mutate(dplyr::group_by(data, Patient), Pop = sum(Cells))

# Calculate ratio
data <- dplyr::mutate(data, Percent = round(Cells/Pop * 100, 1))

# Store timepoints
tp <- data$Timepoint[!duplicated(data$Patient)]

# Create summary 
smry <- dplyr::summarise(dplyr::group_by(data, Celltype), Cells = sum(Cells))
smry$Pop <- sum(smry$Cells)
smry$Percent <- smry$Cells / smry$Pop * 100
smry$Timepoint <- NA
smry$Patient <- "Summary"

data <- dplyr::bind_rows(data, smry)
data$Patient <- factor(data$Patient, unique(data$Patient))

# Plot
ggplot2::ggplot(
  data    = data[order(data$Cells), ],
  mapping = ggplot2::aes(
    x        = Patient,
    y        = Percent,
    fill     = Celltype,
    col      = Celltype,
    label    = Cells
  )
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_label(
    position      = ggplot2::position_stack(vjust = 0.5),
    size          = 4, fill          = "white", 
    label.size    = 0.5, show.legend   = FALSE
  ) +
  ggplot2::scale_color_manual(values = color.celltype[levels(data$Celltype)]) +
  ggplot2::scale_fill_manual(values = color.celltype[levels(data$Celltype)]) +
  ggplot2::theme_classic(base_size = 25) +
  ggplot2::theme(
    legend.position = "right",
    aspect.ratio    = 0.6,
    axis.text.x     = ggplot2::element_text(
      # color = c(color.timepoint[tp], "black"), # BROKEN!!!
      angle = 45, hjust = 1, vjust = 1
    )
  ) +
  ggplot2::labs(
    y = "Proportion of cells [%]", x = NULL
  ) +
  ggplot2::annotate(
    geom = "text",
    label = c(as.character(tp), ""),
    x     = 1:8,
    y     = -10,
    size  = 5,
    color = c(color.timepoint[tp], "black")
  ) +
  ggplot2::guides(
    fill  = ggplot2::guide_legend(
      direction = "horizontal", nrow = 10, title.position = "top"
    ),
    color = ggplot2::guide_none()
  ) +
  ggplot2::expand_limits(y = -15)

# Save plot
fn <- paste0(plot_dir, "BAL_Celltype-proportions_barplot.png")
ggplot2::ggsave(fn, width = 12, height = 6)

# Save dataset -----------------------------------------------------------------
saveRDS(object, out_bal)

# Select the macrophages & re-normalize ----------------------------------------

# Select
cells <- colnames(object)[which(
  object@meta.data$Celltype == "Macrophages"
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

# Normalize mRNA counts
object <- Seurat::NormalizeData(
  object               = object,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

# Find HVGs
object <- Seurat::FindVariableFeatures(
  object           = object,
  selection.method = "vst", 
  nfeatures        = 3000
)

# Calculate low-dimensional embedding & clustering -----------------------------

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
object <- Seurat::FindClusters(object, algorithm=3, resolution = 0.7)
object$Cluster <- object$seurat_clusters

# Select high quality cells & re-normalize -------------------------------------

# Plot
DimPlot(object, label = TRUE) + coord_fixed()
FeaturePlot(object, "CD3D") + coord_fixed() +
  viridis::scale_color_viridis(option = "A", direction = -1)

# Select
cells <- colnames(object)[which(
  object@meta.data$Cluster != "9"
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

# Normalize mRNA counts
object <- Seurat::NormalizeData(
  object               = object,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

# Find HVGs
object <- Seurat::FindVariableFeatures(
  object           = object,
  selection.method = "vst", 
  nfeatures        = 3000
)

# Calculate low-dimensional embedding & clustering -----------------------------

# MNN-corrected PCA
set.seed(1993)
object@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    object@assays$RNA@data[object@assays$RNA@var.features, ],
    batch = object@meta.data$patient,
    d     = 30
  )@int_colData$reducedDims$corrected,
  key        = "PC_",
  assay      = "RNA"
)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:30, seed.use = 1993)

# Clustering
set.seed(1993)
object <- Seurat::FindNeighbors(object = object, dims = 1:30)
object <- Seurat::FindClusters(object, algorithm=3, resolution = 0.7)
object$Cluster <- object$seurat_clusters

# Annotate clusters ------------------------------------------------------------

# Plot
DimPlot(object, label = TRUE) + coord_fixed()

# Annotate
c2l <- list(
  "Monocytes" = c(6,1),
  "Mono/Mφ" = c(2,5),
  "CD163/LGMN-Mφ" = c(0,7),
  "AMφ-1" = c(3),
  "AMφ-2" = c(4),
  "Prolif. AMφ" = c(8)
)
v <- c()
for (i in names(lapply(c2l, length))) {v <- c(v, rep(i, length(c2l[[i]])))}
names(v) <- unlist(c2l)

object$Subtype <- "Low quality"
index <- object$Cluster %in% names(v)
object$Subtype[index] <- v[as.character(object$Cluster[index])]
object$Subtype <- factor(object$Subtype, names(c2l))

# Show clusters ----------------------------------------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = object$Subtype
)

# Select annotation coordinates
ann <- dplyr::summarize(dplyr::group_by(data, col), x = mean(x), y = mean(y))


# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col,
    label = col
  )
) +
  ggplot2::geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = ann, size = 8, label.padding = grid::unit(2, "mm")
  ) +
  ggplot2::scale_color_manual(values = color.cluster) +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = ""
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_Clusters.png")
ggplot2::ggsave(fn, width = 10, height = 6, bg="white")

# Show clusters on 3D UMAP -----------------------------------------------------

# Create 3D umap and fetch clusters
set.seed(1993)
object@reductions$umap_3d <- Seurat::CreateDimReducObject(
  embeddings = uwot::umap(
    X            = object@reductions$pca@cell.embeddings, 
    n_components = 3
  ), 
  key = "UMAP_"
)
data <- tidyr::as_tibble(object@reductions$umap_3d@cell.embeddings)
names(data) <- c("x", "y", "z")
data$col <- object$Subtype

# Plot
rgl::par3d(windowRect = c(50, 50, 800, 800))
rgl::plot3d(
  data, 
  col    = color.cluster[as.character(data$col)], 
  type   = "s", 
  radius = .1,
  axes   = FALSE,
  box    = FALSE,
  xlab   = "",
  ylab   = "",
  zlab   = ""
)
rgl::rgl.lines(c(min(data$x), max(data$x)), c(0, 0), c(0, 0), color = "black", lwd = 4)
rgl::rgl.lines(c(0, 0), c(min(data$y), max(data$y)), c(0, 0), color = "firebrick", lwd = 4)
rgl::rgl.lines(c(0, 0), c(0, 0), c(min(data$z), max(data$z)), color = "seagreen", lwd = 4)

# Save plot
rgl::rgl.snapshot(paste0(plot_dir, "BAL-macrophages_3D-umap-clusters.png"))
rgl::writeWebGL(
  filename = paste0(plot_dir, "BAL-macrophages_3D-umap-clusters.png")
  )
rgl::rgl.close()

# Show differentiation trajectory by slingshot ---------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = object$Subtype
)

# Calculate the slingshot trajectory
set.seed(1993)
trajectory <- slingshot::slingshot(
  data          = data[, c("x", "y")],
  clusterLabels = data$col, 
  allow.breaks  = TRUE
)

data$pst <- trajectory@assays@data$pseudotime[,1]
# Re-shape trajectory data
curves <- list()
for (curve in names(trajectory@metadata$curves)) {
  print(curve)
  curves[[curve]] <- tidyr::as_tibble(trajectory@metadata$curves[[curve]]$s)
  curves[[curve]]$name <- curve
}
curves <- dplyr::bind_rows(curves)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = pst
  )
) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_point(
    data = curves,
    mapping = ggplot2::aes(
      x = x,
      y = y,
      fill = name,
      col = NULL
    ),
    shape  = 21,
    size   = 4, 
    stroke = NA,
    fill = "black"
  ) +
  viridis::scale_color_viridis(direction = -1) +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.18,0.8)
  ) +
  ggplot2::guides(
    #fill = ggplot2::guide_none(),
    color  = ggplot2::guide_colorbar(
      barheight = 1, barwidth = 13, frame.colour = "black", ticks = FALSE, 
      direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = "Pseudotime") +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_Cluster-trajectory.png")
ggplot2::ggsave(fn, width = 10, height = 6, bg="white")

# Show viral mRNA counts -------------------------------------------------------

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col = Matrix::colSums(object@assays$SCoV2@counts)
)
data <- data[order(data$col), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col
  )
) +
  ggplot2::geom_point(size = 3) +
  viridis::scale_color_viridis(direction = -1, option = "A", trans = "log10") +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.2,0.87)
  ) +
  ggplot2::guides(
    col  = ggplot2::guide_colorbar(
      barheight = 1, barwidth = 15, frame.colour = "black", ticks = FALSE,
      direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = "SARS-CoV-2 mRNAs") +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_viral-mRNA.png")
ggplot2::ggsave(fn, width = 10, height = 6, bg="white")

# Show marker genes ------------------------------------------------------------

# 1. Normalized expression on UMAP
marker.genes <- c(
  "S100A8", "MRC1", "FCN1", "CCR2", "TGFB1", "TGFBI", "CCL2", 
  "OLFML2B", "MARCKS", "NRP1", "CD163", "MERTK", "LGMN", "SPP1", "TREM2", 
  "CD9", "APOE", "FBP1", "FABP4", "INHBA", "MME", "MKI67"
)

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2]
  )
for (gene in marker.genes) {
  data[[gene]] <- object@assays$RNA@data[gene, ]
}

# Re-shape and order data for plotting
data <- tidyr::gather(data, "Gene", "Expression", -x, -y)
data$Gene <- factor(data$Gene, levels = unique(data$Gene))
data <- data[order(data$Expression), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = Expression
  )
) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::facet_wrap(~Gene, nrow = 4) +
  viridis::scale_color_viridis(option = "A", direction = -1) +
  ggplot2::theme_void(base_size = 25) +
  ggplot2::theme(
    legend.position = c(0.82, 0.1),
    strip.text      = ggplot2::element_text(face = "italic")
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barheight = 1, barwidth = 15, frame.colour = "black", ticks = FALSE,
      direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = "Norm. Expression") +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    size  = 0.25,
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    size  = 0.25,
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_Cluster-markers_umap.png")
ggplot2::ggsave(fn, width = 12, height = 6, bg="white")

# 2. Scaled expression in dotplot
marker.genes <- list(
  "1 FCN1-Mono"     = c(
    "CD14", "FCGR3A", "S100A8", "S100A12", "FCN1", "VCAN", "FPR1", "CD93", 
    "SELL", "CCR2", "DDX58", "NFKBIA", "SERPINB9", "CXCL8", "NAMPT", 
    "SLC25A37", "APOBEC3A", "TGFB1", "ITGAM", "CLEC5A"
  ),
  "2 Mono/Mp"       = c(
    "MYC", "DDX21", "DPYD", "TYMP", "GNG5", "CCT5", "LSM7", "CALM1"
  ),
  "3 SPP1/LGMN-Mp"  = c(
    "CCL2", "TGFBI", "CTSB", "OLFML2B", "CMKLR1", 
    "CCL13", "CD163", "CD84", "MERTK", "LGMN", "SDC3", "SPP1", "PLA2G7", 
    "NRP1", 
    "TGFBI", "NPL", "CMKLR1", "MARCKS", "MAF", "A2M", "ADAP2",
    "IDH1", "FMN1", "CALR", "PMP22", "HSPA5", "SLAMF8", "SPRED1", "HS3ST1"
  ),
  "4 SPP1/TREM2-Mp" = c(
    "GPNMB", "LIPA", "MATK", "HTRA4", "COL8A2", "ABCC3", 
    "MRC1", "TREM2", "APOE", "ACP5", "GCHFR"
  ),
  "5 INHBA-AMp"     = c(
    "CD9", "EVL", "GSN", "CD52", "APOC1", "MARCO", "FBP1", 
    "GCHFR", "HLA-DRA", "ALOX5AP", "LPL", 
    "ITGB8", "CES1", "INHBA", "GPD1", "MME", "VMO1", "NUPR1", 
    "FN1", "PPIC",
    "S100A13", "PLA2G16", "FABP4"
  ),
  "6 Prolif. AMp"   = c(
    "MKI67", "TOP2A", "CENPA", "BIRC5"
  )
)

marker.genes <- unlist(marker.genes)

# Fetch data
data <- tidyr::tibble(
  Cluster = object$Subtype
)
for (gene in marker.genes) {
  data[[gene]] <- object@assays$RNA@data[gene, ]
}

# Re-shape and order data for plotting
data <- tidyr::gather(data, "Gene", "Expression", -Cluster)
data$Gene <- factor(data$Gene, levels = unique(data$Gene))
data$Expressed <- data$Expression > 0
data$Cells <- 1

data <- dplyr::group_by(data, Cluster, Gene)
data <- dplyr::summarize(
  .data        = data,
  "Expression" = mean(Expression),
  "Expressed"  = sum(Expressed),
  "Cells"      = sum(Cells)
)
data$Percent <- round(data$Expressed / data$Cells * 100, 1)
data <- dplyr::mutate(
  dplyr::group_by(data, Gene), zscore = scale(Expression)[, 1]
)

# Set colorscale limits
limits <- c(-2, 2)
data$zscore[data$zscore > max(limits)] <- max(limits)
data$zscore[data$zscore < min(limits)] <- min(limits)

# Plot
ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x    = Gene,
    y    = Cluster,
    col  = zscore,
    size = Percent
  )
) +
  ggplot2::geom_point() +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barwidth  = 18, barheight = .3,
      ticks = FALSE, frame.colour = "black",
      title.position = "top", title.hjust = 0.5
    ),
    size  = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  ) +
  ggplot2::labs(
    col  = "z-score", size = "% expressed"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    legend.position = "top",
    aspect.ratio    = 1/8,
    axis.text.x     = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, size = 7
    ),
    axis.title      = ggplot2::element_blank()
  ) +
  ggplot2::scale_color_distiller(
    palette = "RdBu", limits = limits
  ) +
  ggplot2::scale_size_area(max_size = 2.5)

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_Cluster-markers_dotplot.png")
ggplot2::ggsave(fn, width = 12, height = 3, bg="white")

# Show clusters across patients ------------------------------------------------

# 1. Barplot of proportions
data <- tidyr::tibble(
  Cluster = object@meta.data$Subtype,
  Patient = object@meta.data$patient,
  Cells   = 1
)
data <- dplyr::summarise(
  dplyr::group_by(data, Cluster, Patient), Cells = sum(Cells)
)
data <- dplyr::mutate(dplyr::group_by(data, Patient), Pop = sum(Cells))
data$Freq <- data$Cells / data$Pop
data$Percent <- round(data$Freq * 100, 1)

smry <- dplyr::summarise(dplyr::group_by(data, Cluster), Cells = sum(Cells))
smry$Pop <- sum(smry$Cells)
smry$Percent <- smry$Cells / smry$Pop * 100
smry$Timepoint <- NA
smry$Patient <- "Summary"

data <- dplyr::bind_rows(data, smry)
data$Patient <- factor(
  x = data$Patient, 
  levels = c(
    as.character(unique(
      object@meta.data$patient[order(object@meta.data$dpso)]
      )), "Summary"
  )
)

tp <- unique(object@meta.data[, c("patient", "dpso")])
tp <- rbind(tp, data.frame(patient="Summary", dpso=""))
index <- match(levels(data$Patient), tp$patient)
tp <- tp[index, ]

ggplot2::ggplot(
  data    = data[order(data$Cells), ],
  mapping = ggplot2::aes(
    x        = Patient,
    y        = Percent,
    fill     = Cluster,
    col      = Cluster,
    label    = Cells
  )
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_label(
    position      = ggplot2::position_stack(vjust = 0.5),
    size          = 4, fill          = "white", 
    label.size    = 0.5, show.legend   = FALSE
  ) +
  ggplot2::scale_color_manual(values = color.cluster[levels(data$Cluster)]) +
  ggplot2::scale_fill_manual(values = color.cluster[levels(data$Cluster)]) +
  ggplot2::theme_classic(base_size = 25) +
  ggplot2::theme(
    legend.position = "right",
    aspect.ratio    = 0.7,
    axis.text.x     = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    )
  ) +
  ggplot2::labs(
    y = "Proportion of cells [%]", x = NULL
  ) +
  ggplot2::annotate(
    geom = "text",
    label = tp$dpso,
    x     = 1:8,
    y     = -11,
    size  = 5,
    color = color.timepoint[tp$dpso]
  ) +
  ggplot2::expand_limits(y = -15) +
  ggplot2::guides(
    fill  = ggplot2::guide_legend(
      direction = "horizontal", nrow = 10, title.position = "top"
    ),
    color = ggplot2::guide_none()
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_barplot-cluster-patient.png")
ggplot2::ggsave(fn, width = 11, height = 6, bg="white")

# 2. UMAP embedding split by patients

# Fetch data
data <- tidyr::tibble(
  x = object@reductions$umap@cell.embeddings[, 1],
  y = object@reductions$umap@cell.embeddings[, 2],
  col   = object$Subtype,
  group = object$patient
)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col,
    label = col
  )
) +
  ggplot2::geom_point(size = .5) +
  ggplot2::facet_wrap(~group, nrow = 2) +
  ggplot2::scale_color_manual(values = color.cluster) +
  ggplot2::theme_void(base_size = 25) +
  ggplot2::theme(
    legend.position = ""
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.075, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.075, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_umap-cluster-patient.png")
ggplot2::ggsave(fn, width = 12, height = 4, bg="white")

# ------------------------------------------------------------------------------
# Differential expression between clusters

# DE test
markers <- scran::findMarkers(
  x         = object@assays$RNA@data,
  groups    = object@meta.data$Subtype,
  pval.type = "some",
  test.type = "wilcox",
  direction = "up", 
  block     = object@meta.data$patient
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$cluster <- i
  markers[[i]]$gene <- row.names(markers[[i]])
}
markers <- dplyr::bind_rows(as.list(markers))
markers <- markers[,c(12,13,1:11,14)]

# Save results of DE test
fn <- paste0(plot_dir, "BAL-macrophages_Cluster-markers.csv")
write.csv(markers, fn, row.names = FALSE)
fn <- paste0(plot_dir, "BAL-macrophages_Cluster-markers_dotplot.csv")
write.csv(markers[markers$gene %in% marker.genes, ], fn, row.names = FALSE)

# Select genes based on FDR cutoff
cutoff <- 1e-15

# Select unique DE genes
de <- markers[markers$FDR < cutoff, ]
de <- de[!duplicated(de$gene), ]

# Limit proliferation DE genes
v <- de$gene[de$cluster == "Prolif. AMφ"]
de <- de[which(!de$gene %in% v[c(100:length(v))]), ]

# Re-order factor levels
de$cluster <- factor(de$cluster, levels(object$Subtype))

# Fetch normalized count data
data <- object@assays$RNA@data[
  de$gene, order(object$Subtype, object@meta.data$patient)
]

# Scale (z-scores)
data <- t(scale(t(as.matrix(data))))

# Create column annotations/gaps
cann <- data.frame(
  Cluster = object@meta.data$Subtype,
  Patient = object@meta.data$patient
)
row.names(cann) <- row.names(object@meta.data)

color.cann <- list(
  Cluster = color.cluster,
  Patient = RColorBrewer::brewer.pal(
    name = "Dark2", n = length(unique(object@meta.data$patient))
  )
)
names(color.cann$Patient) <- intersect(
  levels(object@meta.data$patient), unique(object@meta.data$patient)
)

# Create row annotations/gaps
rann <- data.frame(
  Cluster = de$cluster
)

# Set colorscale limits
limits <- c(-2, 2)
breaks <- seq(from = min(limits), to = max(limits), length.out = 100)

# Plot
plot <- pheatmap::pheatmap(
  mat    = data, 
  breaks = breaks,
  color  = colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
  )(length(breaks)),
  annotation_col    = cann,
  annotation_colors = color.cann,
  cluster_rows         = FALSE,
  cluster_cols         = FALSE,
  show_rownames        = FALSE,
  show_colnames        = FALSE, 
  annotation_names_col = FALSE,
  gaps_col             = head(as.numeric(cumsum(table(cann$Cluster))), -1),
  gaps_row             = head(as.numeric(cumsum(table(rann$Cluster))), -1), 
  silent = TRUE
)

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_DE-heatmap.png")
ggplot2::ggsave(fn, plot, width = 12, height = 6)

# ChIP-seq enrichment analysis (ChEA3) based on DE genes ----------------------

# Select non-unique DE genes
de <- markers[markers$FDR < cutoff, ]
de$cluster <- factor(de$cluster, levels(object$Subtype))
genes <- split(de$gene, de$cluster)

# ChEA3 query for each cluster
result <- list()
for (i in names(genes)) {
  print(i)
  # Skip empty categories
  if(length(genes[[i]]) > 0) {
    # Prepare input as list
    payload = list(query_name = "myQuery", gene_set = genes[[i]])
    
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
result$cluster <- factor(result$cluster, levels(object$Subtype))
result <- result[order(result$cluster, result$Score), ]

# Select TFs by mean rank
tfs <- result$TF[result$Score < 30]
data <- result[which(result$TF %in% tfs), ]
data$TF <- factor(data$TF, levels = unique(tfs))
data$cluster <- factor(data$cluster, levels(object$Subtype))

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
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    aspect.ratio = 0.25,
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
fn <- paste0(plot_dir, "BAL-macrophages_TF-enrichment.png")
ggplot2::ggsave(fn, width = 12, height = 3, bg="white")
fn <- paste0(plot_dir, "BAL-macrophages_TF-enrichment.csv")
write.csv(result[, c(7,3,2,4,8,9,10,5,6,1)], fn, row.names = FALSE)

# Enrichment of COVID-19 gene sets ---------------------------------------------

# Retrieve gene set dictionary
file <- tempfile()
download.file(
  "https://syncandshare.desy.de/index.php/s/ADoDP7imDpAm3Dw/download", file
)
dictionary <- read.csv(file)

# Select disease and gene set size
key <- "COVID-19"
n <- 50

dict <- list()
genesets <- list()
for (i in unique(dictionary$term[dictionary$disease == key])) {
  dict[[i]] <- head(dictionary[dictionary$term == i, ], n)
  genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
}
dict <- dplyr::bind_rows(dict)

# Overrepresentation analysis

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

# Order factors
term2ref <- unique(dict[, c("term", "ref")])$ref
names(term2ref) <- unique(dict[, c("term", "ref")])$term
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(11,1,8,12,5,6,2,13,14,7,3,15,4,16,9,10)]
)

# Select colors
cols <- c(
  "Liao et al."   = "seagreen",
  "Grant et al." = "navy",
  "Bharat et al." = "orange"
)
ref.color <- cols[term2ref[levels(data$Geneset)]]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = Geneset,
    y = Cluster,
    fill = -log10(p.adjust)
  )
) +
  ggplot2::geom_tile(col = "white", size = 1) +
  ggplot2::scale_fill_gradient2(low = "white", high = "red") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    legend.position = "top",
    aspect.ratio = 0.4,
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.ticks.x = ggplot2::element_line(
      color = ref.color, size = 7, arrow = grid::arrow(
        length = grid::unit(8, "pt"), angle = 90
      )
    )
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::expand_limits(y = -1.5) +
  ggplot2::annotate(
    geom = "text",
    y = rep(-0.4, 3),
    x = c(3,8.5,14),
    size  = 6,
    col   = cols,
    label = names(cols)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      barwidth = 39, barheight = 1, ticks = FALSE, frame.colour = "black", 
      title.position = "top", title.hjust = 0.5
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_COVID-19-geneset-enrichment.png")
ggplot2::ggsave(fn, width = 11, height = 6, bg="white")

# Module scores

# 1. UMAP embedding with cell-based scores
object@meta.data <- object@meta.data[, which(
  !stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
)]

object <- Seurat::AddModuleScore(
  object   = object, 
  features = genesets, 
  name     = "geneset",
  assay    = "RNA", 
  seed     = 1993
)

scores <- object@meta.data[, which(
  stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
)]
names(scores) <- names(genesets)

scores$x <- object@reductions$umap@cell.embeddings[, 1]
scores$y <- object@reductions$umap@cell.embeddings[, 2]

data <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(11,1,8,12,5,6,2,13,14,7,3,15,4,16,9,10)]
)

# Set colorscale limits
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)

data <- data[order(data$Score), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(size = 0.25) +
  ggplot2::facet_wrap(~Geneset, nrow = 4)+
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::theme_void(base_size = 20) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, 
      frame.colour = "black", ticks.colour = "black"
    )
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    size  = 0.25,
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment",
    size  = 0.25,
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_COVID-19-geneset-score-umap.png")
ggplot2::ggsave(fn, width = 9, height = 6, bg="white")

# 2. Violin plots as cluster-based population summaries

# Re-shape data
data <- scores
data$Cluster <- object@meta.data$Subtype
data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(11,1,8,12,5,6,2,13,14,7,3,15,4,16,9,10)]
)
data$col <- data$Score

# Set colorscale limits
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)

data <- dplyr::mutate(dplyr::group_by(data, Geneset), ms = mean(Score))

# Calculate p-values
ps <- unique(data[, c("Geneset", "Cluster")])
ps$p.value <- NA
for (i in levels(data$Geneset)) {
  print(i)
  for (ii in levels(data$Cluster)) {
    print(ii)
    ps$p.value[ps$Geneset == i & ps$Cluster == ii] <- wilcox.test(
      x = data$Score[data$Geneset == i & data$Cluster == ii],
      y = data$Score[data$Geneset == i], 
      alternative = "greater"
    )$p.value
  }
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))

ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter"
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) +
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) +
  ggplot2::facet_wrap(~Geneset, nrow = 4) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::scale_fill_manual(values = color.cluster) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 8, barheight = 1, frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_COVID-19-geneset-score-violins.png")
ggplot2::ggsave(fn, width = 12, height = 8)

# Enrichment of fibrosis gene sets ---------------------------------------------

# Select disease and gene set size
key <- "IPF"
n <- 50

dict <- list()
genesets <- list()
for (i in unique(dictionary$term[dictionary$disease == key])) {
  dict[[i]] <- head(dictionary[dictionary$term == i, ], n)
  genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
}
dict <- dplyr::bind_rows(dict)

# Overrepresentation analysis

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

# Order factors
data$Geneset <- factor(data$Geneset, unique(data$Geneset))
term2ref <- unique(dict[, c("term", "ref")])$ref
names(term2ref) <- unique(dict[, c("term", "ref")])$term
data$Geneset <- factor(data$Geneset, unique(data$Geneset))
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,5,8,4,6,12,9)]
)

# Select colors
cols <- c(
  "Morse et al."   = "seagreen",
  "Adams et al."   = "navy",
  "Reyfman et al." = "dodgerblue",
  "Ayaub et al."   = "purple"
)
ref.color <- cols[term2ref[levels(data$Geneset)]]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = Geneset,
    y = Cluster,
    fill = -log10(p.adjust)
  )
) +
  ggplot2::geom_tile(col = "white", size = 1) +
  ggplot2::scale_fill_gradient2(low = "white", high = "red") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    legend.position = "top",
    aspect.ratio = 1/3,
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.ticks.x = ggplot2::element_line(
      color = ref.color, size = 7, arrow = grid::arrow(
        length = grid::unit(8, "pt"), angle = 90
      )
    )
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::expand_limits(y = -1.5) +
  ggplot2::annotate(
    geom = "text",
    y = rep(-0.4, 4),
    x = c(2,5,8,11),
    size  = 6,
    col   = unique(ref.color),
    label = unique(names(ref.color))
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      barwidth = 40, barheight = 1, ticks = FALSE, frame.colour = "black", 
      title.position = "top", title.hjust = 0.5
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_IPF-geneset-enrichment.png")
ggplot2::ggsave(fn, width = 11, height = 6)

# Module scores

# 1. UMAP embedding with cell-based scores
object@meta.data <- object@meta.data[, which(
  !stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
)]

object <- Seurat::AddModuleScore(
  object   = object, 
  features = genesets, 
  name     = "geneset",
  assay    = "RNA", 
  seed     = 1993
)

scores <- object@meta.data[, which(
  stringr::str_detect(string = names(object@meta.data), pattern = "geneset")
)]
names(scores) <- names(genesets)

scores$x <- object@reductions$umap@cell.embeddings[, 1]
scores$y <- object@reductions$umap@cell.embeddings[, 2]

data <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
data$Geneset <- factor(data$Geneset, levels = unique(data$Geneset))
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
)

# Set colorscale limits
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)

data <- data[order(data$Score), ]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(size = 0.25) +
  ggplot2::facet_wrap(~Geneset, nrow = 3)+
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::theme_void(base_size = 20) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, 
      frame.colour = "black", ticks.colour = "black"
    )
  ) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.15, "inches"), type   = "closed"
    )
  )

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_IPF-geneset-score-umap.png")
ggplot2::ggsave(fn, width = 12, height = 6, bg="white")

# 2. Violin plots as cluster-based population summaries

data <- scores
data$Cluster <- object@meta.data$Subtype

data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)

data$Geneset <- factor(data$Geneset, unique(data$Geneset))
data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(1,11,10,2,3,7,8,12,4,6,5,9)]
)

data$col <- data$Score
data$col[data$col > max(limits)] <- max(limits)
data$col[data$col < min(limits)] <- min(limits)

data <- dplyr::mutate(dplyr::group_by(data, Geneset), ms = mean(Score))

# Calculate p-values
ps <- unique(data[, c("Geneset", "Cluster")])
ps$p.value <- NA
for (i in levels(data$Geneset)) {
  print(i)
  for (ii in levels(data$Cluster)) {
    print(ii)
    ps$p.value[ps$Geneset == i & ps$Cluster == ii] <- wilcox.test(
      x = data$Score[data$Geneset == i & data$Cluster == ii],
      y = data$Score[data$Geneset == i], 
      alternative = "greater"
    )$p.value
  }
}
ps$p.adjust <- p.adjust(ps$p.value)
ps$label <- round(-log10(ps$p.adjust))

ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Cluster,
    y    = Score,
    col  = col,
    fill = Cluster
  )
) +
  ggplot2::geom_point(
    position = "jitter"
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "darkgrey", size = 3
  ) +
  ggplot2::geom_hline(
    mapping = ggplot2::aes(yintercept = ms), col = "black", size = 0.1
  ) +
  ggplot2::geom_violin(scale = "width", draw_quantiles = 0.5) +
  ggplot2::facet_wrap(~Geneset, nrow = 3) +
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::scale_fill_manual(values = color.cluster) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    #strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 8, barheight = 1, frame.colour = "black", ticks = FALSE,
      direction = "horizontal"
    ),
    fill = ggplot2::guide_legend(
      order = 1, direction = "horizontal", title.position = "top"
    )
  ) +
  ggplot2::labs(col = NULL, x = NULL) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_label(
    data = ps,
    mapping = ggplot2::aes(
      label = label, y = 2.3, fill = label, col = NULL
    ),
    size = 4,
    color = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# Save plot
fn <- paste0(plot_dir, "BAL-macrophages_IPF-geneset-score-violins.png")
ggplot2::ggsave(fn, width = 12, height = 6)

# Save dataset -----------------------------------------------------------------
saveRDS(object, out_bal_macrophages)
