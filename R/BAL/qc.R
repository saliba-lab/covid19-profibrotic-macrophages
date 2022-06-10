# Global varibles --------------------------------------------------------------

file <- "data/BAL.Rds"
dest <- "analysis/BAL/qc/"
dir.create(dest, recursive = TRUE)

# Load data --------------------------------------------------------------------

ds <- readRDS(file)

# Filtering based on quality metrics -------------------------------------------

# Add percentage of mitochondrial genes
index <- grep("^MT-", rownames(ds[["RNA"]]@counts))
ds$percent.mito <- round(
  Matrix::colSums(ds[["RNA"]]@counts[index, ]) / ds$nCount_RNA, 3
) * 100

# Plot
data <- data.frame(
  libsize = ds$nCount_RNA,
  genes = ds$nFeature_RNA,
  percent.mito = ds$percent.mito,
  donor = ds$patient
)
data <- tidyr::gather(data, "key", "value", -libsize, -donor)
ggplot2::ggplot(data, ggplot2::aes(value, libsize, col = donor)) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::facet_wrap(~key, scales = "free") +
  ggplot2::theme_classic(20) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 5))
  ) +
  ggplot2::labs(x = NULL, y = "Library size")

fn <- paste0(dest, "libsize-genes-pMt", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Set gene filter
index <- which(Matrix::rowSums(ds[["RNA"]]@counts) != 0)
ds[["RNA"]] <- subset(ds[["RNA"]], features = index)

# Set cell filter
index <- colnames(ds)[which(
  dplyr::between(ds$nFeature_RNA, 0, Inf) &
  dplyr::between(ds$nCount_RNA, 0, Inf) &
  dplyr::between(ds$percent.mito, 0, Inf)
)]
ds <- subset(ds, cells = index)

# Plot
data <- data.frame(
  libsize = ds@meta.data$nCount_RNA,
  genes = ds$nFeature_RNA,
  percent.mito = ds$percent.mito,
  donor = ds$patient
)
data <- tidyr::gather(data, "key", "value", -libsize, -donor)
ggplot2::ggplot(data, ggplot2::aes(value, libsize, col = donor)) +
  ggplot2::geom_point(size = 0.1) +
  ggplot2::facet_wrap(~key, scales = "free") +
  ggplot2::theme_classic(20) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 5))
  ) +
  ggplot2::labs(x = NULL, y = "Library size")

fn <- paste0(dest, "libsize-genes-pMt_filtered", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Normalize data ---------------------------------------------------------------

ds <- Seurat::NormalizeData(ds, assay = "RNA")

# HVG selection & scaling ------------------------------------------------------

ds <- Seurat::FindVariableFeatures(ds, nfeatures = 3000)
ds <- Seurat::ScaleData(ds)

# ------------------------------------------------------------------------------
# Calculate low-dimensional embedding & clustering

npcs <- 40

# PCA (MNN-corrected)
hvg <- ds[["RNA"]]@var.features
ds[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds[["RNA"]]@data[hvg, ], batch = ds$patient, d = npcs
  )@int_colData$reducedDims$corrected,
  key = "PC_", assay = "RNA"
)

# UMAP
ds <- Seurat::RunUMAP(ds, dims = 1:npcs)

# Clustering
ds <- Seurat::FindNeighbors(ds, dims = 1:npcs)
ds <- Seurat::FindClusters(ds, algorithm = 1, resolution = 0.3)

# ------------------------------------------------------------------------------
# Select high quality classical monocytes

# Plot clusters & percent mitochondrial counts
cowplot::plot_grid(
  Seurat::DimPlot(
    ds, reduction = "umap", group.by = "seurat_clusters", 
    label = TRUE, label.box = TRUE, pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::guides(col = ggplot2::guide_none()),
  Seurat::DimPlot(
    ds, reduction = "umap", group.by = "patient", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15), rel_widths = c(.41, .55)
)

fn <- paste0(dest, "umap", "_", "cluster-patient", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 4, bg = "white")

# Plot count depths
cowplot::plot_grid(
  Seurat::FeaturePlot(
    ds, reduction = "umap", features = "nCount_RNA", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1),
  Seurat::FeaturePlot(
    ds, reduction = "umap", features = "percent.mito", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)

fn <- paste0(dest, "umap", "_", "libsize-pMt", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5, bg = "white")

# Plot marker genes
cowplot::plot_grid(
  Seurat::FeaturePlot(
    ds, reduction = "umap", 
    features = "S100A8", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1),
  Seurat::FeaturePlot(
    ds, reduction = "umap", 
    features = "CD3D", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)

fn <- paste0(dest, "umap", "_", "marker-genes", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5, bg = "white")

# Plot marker genes
cowplot::plot_grid(
  Seurat::FeaturePlot(
    ds, reduction = "umap", 
    features = "MKI67", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1),
  Seurat::FeaturePlot(
    ds, reduction = "umap", 
    features = "LGMN", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)
