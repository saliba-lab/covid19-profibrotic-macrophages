# Global varibles --------------------------------------------------------------

file <- "data/monocytes.Rds"
dest <- "analysis/monocytes/qc/"
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
  donor = ds$donor
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
  dplyr::between(ds$nFeature_RNA, 0, 5000),
  dplyr::between(ds$nCount_RNA, 0, 30000),
  dplyr::between(ds$percent.mito, 0, Inf)
)]
ds <- subset(ds, cells = index)

# Plot
data <- data.frame(
  libsize = ds@meta.data$nCount_RNA,
  genes = ds$nFeature_RNA,
  percent.mito = ds$percent.mito,
  donor = ds$donor
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

# mRNA counts
ds <- Seurat::NormalizeData(ds, assay = "RNA")

# Antibody counts
ds <- Seurat::NormalizeData(ds, assay = "HTO", method = "CLR")
ds <- Seurat::NormalizeData(ds, assay = "ADT", method = "CLR")

# Demultiplex conditions from HTO counts ---------------------------------------

# Retrieve HTO counts
data <- as.data.frame(t(as.matrix(ds[["HTO"]]@counts)))
colnames(data) <- stringr::str_split(colnames(data), "-", simplify = TRUE)[, 1]
set.seed(42)
data$id <- factor(kmeans(log1p(data), 5)$cluster)
data$bc <- colnames(ds)
data <- tidyr::gather(data, "Hashtag", "count", -bc, -id)

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = Hashtag,
    y = count,
    col = id
  )
) +
  ggplot2::geom_point(position = "jitter", size = .5) +
  ggplot2::geom_violin(scale = "width", ggplot2::aes(col = NULL), fill = NA) +
  ggplot2::theme_classic(20) +
  ggplot2::labs(x = NULL) +
  ggplot2::scale_y_continuous(trans = "log10")

fn <- paste0(dest, "hashtag-counts", "_", "kmeans", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Select count thresholds
thresh <- c(
  Hashtag5 = Inf,
  Hashtag6 = 45,
  Hashtag7 = 50,
  Hashtag8 = 60,
  Hashtag9 = 30
)

# Apply threshold
data$thresh <- thresh[data$Hashtag]
data$label <- "Negative"
index <- data$count > data$thresh
data$label[index] <- data$Hashtag[index]

# Classification function
classify <- function(x) {
  if (sum(x) < 1) {return("Negative")} else {
    if (sum(x) > 1) {return("Doublet")}  else {
      return(names(sort(x, decreasing = TRUE)[1]))
    }
  }
}

# Classification
cft <- as.matrix(table(data$bc, data$label))
index <- which(colnames(cft) == "Negative")
cft <- cft[colnames(ds), -index]
ds@meta.data$hash_id <- factor(apply(cft, 1, classify))

# Plot
data$label <- ds$hash_id[data$bc]
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = Hashtag,
    y = count,
    col = label
  )
) +
  ggplot2::geom_point(position = "jitter", size = .5) +
  ggplot2::geom_violin(scale = "width", ggplot2::aes(col = NULL), fill = NA) +
  ggplot2::theme_classic(20) +
  ggplot2::labs(x = NULL) +
  ggplot2::scale_y_continuous(trans = "log10") +
  ggplot2::geom_point(
    ggplot2::aes(y = thresh, group = Hashtag), 
    col = "black", shape = 95, size = 5, stroke = 15
    )

fn <- paste0(dest, "hashtag-counts", "_", "threshold", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Convert hashtags to stimulation conditions
label <- c(
  Hashtag5 = "SARS-CoV-2 MOI15",
  Hashtag6 = "Control",
  Hashtag7 = "SARS-CoV-2",
  Hashtag8 = "R848",
  Hashtag9 = "3p-hpRNA",
  Doublet  = "Doublet",
  Negative = "Negative"
)

ds@meta.data$condition <- factor(
  label[as.character(ds@meta.data$hash_id)], 
  levels = c("Control", "SARS-CoV-2", "3p-hpRNA", "R848", "Negative", "Doublet")
)

# HVG selection & scaling ------------------------------------------------------

ds <- Seurat::FindVariableFeatures(ds, nfeatures = 3000)
ds <- Seurat::ScaleData(ds)

# ------------------------------------------------------------------------------
# Calculate low-dimensional embedding & clustering

# PCA
ds <- Seurat::RunPCA(ds, npcs = 30)

# UMAP
ds <- Seurat::RunUMAP(ds, dims = 1:30)

# Clustering
ds <- Seurat::FindNeighbors(ds, dims = 1:30)
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
    ds, reduction = "umap", group.by = "condition", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15), rel_widths = c(.41, .55)
)

fn <- paste0(dest, "umap", "_", "cluster-condition", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5, bg = "white")

# Plot count depths
cowplot::plot_grid(
  Seurat::FeaturePlot(
    ds, reduction = "umap", features = "nCount_RNA", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1) +
    ggplot2::theme(legend.position = c(0.8, 0.8)),
  Seurat::FeaturePlot(
    ds, reduction = "umap", features = "percent.mito", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    viridis::scale_color_viridis(option = "A", direction = -1) +
    ggplot2::theme(legend.position = c(0.8, 0.8))
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
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1) +
    ggplot2::theme(legend.position = c(0.8, 0.8)),
  Seurat::FeaturePlot(
    ds, reduction = "umap", 
    features = "CD14", pt.size = .1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1) +
    ggplot2::theme(legend.position = c(0.8, 0.8))
)

fn <- paste0(dest, "umap", "_", "marker-genes", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5, bg = "white")
