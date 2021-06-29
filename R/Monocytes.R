################################################################################
# scRNA-seq analysis of stimulated monocytes

# ------------------------------------------------------------------------------
# Download data from public repository

# Barcodes
con <- "https://nubes.helmholtz-berlin.de/s/kr3pDEnTs8Zig59/download"
file <- "~/Downloads/monocytes-barcodes.tsv"
download.file(url = con, destfile = file)
barcodes <- readLines(file)

# Features
con <- "https://nubes.helmholtz-berlin.de/s/kDDPYYTjj7CMfEB/download"
file <- "~/Downloads/monocytes-features.tsv"
download.file(url = con, destfile = file)
features <- readr::read_tsv(file, col_names = c("ENSEMBL", "SYMBOL", "TYPE"))

# Matrix
con <- "https://nubes.helmholtz-berlin.de/s/ecHWyJy3eXxgYYc/download"
file <- "~/Downloads/monocytes-matrix.mtx"
download.file(url = con, destfile = file)
matrix <- Matrix::readMM(gzfile(file))
colnames(matrix) <- barcodes
rownames(matrix) <- features$ENSEMBL

# Donor information
con <- "https://nubes.helmholtz-berlin.de/s/xsYCzt6wpMDZafA/download"
file <- "~/Downloads/monocytes-donors-1.tsv"
download.file(url = con, destfile = file)
donors <- list(
  "1" = readr::read_tsv(file)
)

con <- "https://nubes.helmholtz-berlin.de/s/FoGG7wAJnDaZAXe/download"
file <- "~/Downloads/monocytes-donors-2.tsv"
download.file(url = con, destfile = file)
donors[["2"]] <- readr::read_tsv(file)
donors[["2"]]$barcode <- stringr::str_replace(
  donors[["2"]]$barcode, "1", "2"
  )

donors <- dplyr::bind_rows(donors)

# Check overlap between donors classification & matrix
table(donors$barcode == barcodes)

# ------------------------------------------------------------------------------
# Re-shape data into Seurat object

# Separate matrices of monocyte gene, SARS-CoV-2 gene and HTO-derived counts
adt <- matrix[
  features$ENSEMBL[which(features$TYPE == "Antibody Capture" & 
                           !stringr::str_detect(features$ENSEMBL, "Hashtag"))]
  , ]
hashtags <- matrix[which(stringr::str_detect(rownames(matrix), "Hashtag")), ]
scov2 <- matrix[which(stringr::str_detect(rownames(matrix), "SCoV2")), ]
matrix <- matrix[which(stringr::str_detect(rownames(matrix), "ENSG")), ]

# Compare matrix size
lapply(
  c("Gene expression"       = matrix, 
    "SARS-CoV-2 mRNAs"      = scov2, 
    "Hashtags"              = hashtags,
    "Antibody-derived tags" = adt
    ), 
  dim
  )

# Create Seurat object and remove raw matrix
object <- Seurat::CreateSeuratObject(
  counts  = matrix, 
  assay   = "RNA",
  project = "Stimulated monocytes"
)
rm(matrix)

# Add viral counts and hashtags as assays
object@assays[["SCoV2"]] <- Seurat::CreateAssayObject(counts = scov2)
rm(scov2)

object@assays[["HTO"]] <- Seurat::CreateAssayObject(counts = hashtags)
rm(hashtags)

object@assays[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)
rm(adt)

# Add donors as meta data
object@meta.data$donor <- as.character(c(
  "0" = "A", "1" = "B", "1/0" = "Doublet", "0/1" = "Doublet"
  )[donors$assignment])

# Store features in the object
object@misc$features <- features

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
cells <- rownames(object@meta.data)[which(
  object@meta.data$nFeature_RNA > 0 &
    object@meta.data$nFeature_RNA < 5000 &
    object@meta.data$nCount_RNA > 0 &
    object@meta.data$nCount_RNA < 30000 &
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
# Normalize

# Hashtags
object <- Seurat::NormalizeData(
  object = object, assay = "HTO", method = "CLR"
  )

# Protein antibodies
object <- Seurat::NormalizeData(
  object = object, assay = "ADT", normalization = "CLR", margin = 2
)

# Normalize mRNA counts
object <- Seurat::SCTransform(
  object = object, assay = "RNA", new.assay.name = "SCT", 
  variable.features.n = 2000
  )

# ------------------------------------------------------------------------------
# Demultiplex stimulation conditions from hashtags

# Retrieve normalized hashtag counts
set.seed(1993)
df <- tidyr::as_tibble(t(as.matrix(object@assays$HTO@data)))
df$bc <- rownames(object@meta.data)
df <- tidyr::gather(df, "Hashtag", "CLR", -bc)

# K-means clustering for each hashtag
df <- dplyr::mutate(
  dplyr::group_by(df, Hashtag), id = factor(kmeans(CLR, 3)$cluster)
)

# Plot
ggplot2::ggplot(
  data = df,
  mapping = ggplot2::aes(
    x = Hashtag,
    y = CLR,
    col = id
  )
) +
  ggplot2::geom_point(position = "jitter") +
  ggplot2::geom_violin(scale = "width", ggplot2::aes(col = NULL), fill = NA)

# Select true labelling events
df$HTO_max <- "Negative"
df$HTO_max[df$Hashtag == "Hashtag6" & df$id == "3"] <- "6"
df$HTO_max[df$Hashtag == "Hashtag7" & df$id == "3"] <- "7"
df$HTO_max[df$Hashtag == "Hashtag8" & df$id == "1"] <- "8"
df$HTO_max[df$Hashtag == "Hashtag9" & df$id == "2"] <- "9"

# Classification function
classify <- function(x) {
  if (sum(x) < 1) {return("Negative")} else {
    if (sum(x) > 1) {return("Doublet")}  else {
      return(names(sort(x, decreasing = TRUE)[1]))
    }
  }
}

# Classification
dat <- as.matrix(table(df$bc, df$HTO_max))
dat <- dat[colnames(object), -5]
object@meta.data$HTO_class <- factor(apply(dat, 1, classify))

# Convert hashtags to stimulation conditions
stimulation_hashtags <- c(
  "5"        = "SARS-CoV-2 MOI15",
  "6"        = "Control",
  "7"        = "SARS-CoV-2",
  "8"        = "R848",
  "9"        = "3p-hpRNA",
  "Doublet"  = "Doublet",
  "Negative" = "Negative"
)
object@meta.data$condition <- factor(
  stimulation_hashtags[as.character(object@meta.data$HTO_class)], 
  levels = c("Control", "SARS-CoV-2", "3p-hpRNA", "R848", "Negative", "Doublet")
  )

# ------------------------------------------------------------------------------
# Calculate low-dimensional embedding & clustering

# PCA
object <- Seurat::RunPCA(object = object, npcs = 30, seed.use = 1993)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:30, seed.use = 1993)

# Clustering
set.seed(1993)
object <- Seurat::FindNeighbors(object = object, dims = 1:30)
object <- Seurat::FindClusters(
  object = object, algorithm = 4, resolution = 0.9
  )

# ------------------------------------------------------------------------------
# Select high quality classical monocytes

# Plot count depths
cowplot::plot_grid(
  Seurat::FeaturePlot(
    object, reduction = "umap", features = "nFeature_RNA", pt.size = 1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1),
  Seurat::FeaturePlot(
    object, reduction = "umap", features = "nCount_RNA", pt.size = 1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)

# Plot clusters & percent mitochondrial counts
cowplot::plot_grid(
  Seurat::DimPlot(object, reduction = "umap", label = TRUE, pt.size = 1) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::guides(col = ggplot2::guide_none()),
  Seurat::FeaturePlot(
    object, reduction = "umap", features = "percent.mt", pt.size = 1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)

# Plot marker genes
cowplot::plot_grid(
  Seurat::FeaturePlot(
    object, reduction = "umap", 
    features = features$ENSEMBL[features$SYMBOL == "S100A8"], pt.size = 1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1),
  Seurat::FeaturePlot(
    object, reduction = "umap", 
    features = features$ENSEMBL[features$SYMBOL == "CD14"], pt.size = 1
  ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    viridis::scale_color_viridis(option = "A", direction = -1)
)

# Select clusters
monocyte_clusters <- c("5", "10", "18", "19")
object@meta.data$monocytes <- as.logical(
  object@meta.data$seurat_clusters %in% monocyte_clusters
)

# Subset object based on selection
cells <- row.names(object@meta.data[which(
  object@meta.data$monocytes & 
    !object@meta.data$condition %in% c("Doublet", "Negative")
  ), ])
object <- subset(object, cells = cells)

# ------------------------------------------------------------------------------
# Filter genes & re-normalize

# Set gene filter 
genes <- list()
for (assay in names(object@assays)) {
  print(assay)
  genes[[assay]] <- names(which(
    Matrix::rowSums(slot(object@assays[[assay]], "counts")) > 0
  ))
}
genes <- as.character(unlist(genes))

# Apply filter
object <- subset(object, features = genes)

# Protein antibodies
object <- Seurat::NormalizeData(
  object = object, assay = "ADT", normalization = "CLR", margin = 2
)

# Normalize mRNA counts
object <- Seurat::SCTransform(
  object = object, assay = "RNA", new.assay.name = "SCT", 
  variable.features.n = 2000
)

# ------------------------------------------------------------------------------
# Compute low dimensional embedding

# PCA
object <- Seurat::RunPCA(object = object, npcs = 30, seed.use = 1993)

# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:30, seed.use = 1993)

# ------------------------------------------------------------------------------
# Show experimental conditions

# Re-level factor
object@meta.data$condition <- factor(
  object@meta.data$condition, unique(object@meta.data$condition)[c(1,2,4,3)]
  
)

# Choose colors
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7"
)
color.condition <- c(
  "Control"    = cbPalette[1],
  "SARS-CoV-2" = cbPalette[7],
  "3p-hpRNA"   = cbPalette[4],
  "R848"       = cbPalette[6]
)

# Select data
data <- tidyr::as_tibble(object@reductions$umap@cell.embeddings)
names(data) <- c("x", "y")
data$col <- object@meta.data$condition

# Set text location
ann <- data.frame(
  x = c( 4,  3,  3, 10),
  y = c( 4, -1, -4, -2),
  col = levels(data$col)
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
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = ann, size = 15
  ) +
  ggplot2::scale_color_manual(values = color.condition) +
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

# ------------------------------------------------------------------------------
# Show SARS-CoV-2 counts

# Select data
data <- tidyr::as_tibble(object@reductions$umap@cell.embeddings)
names(data) <- c("x", "y")
data$col <- Matrix::colSums(object@assays$SCoV2@counts)

# Plot
ggplot2::ggplot(
  data = data[order(data$col), ],
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col
  )
) +
  ggplot2::geom_point(size = 2) +
  viridis::scale_color_viridis(option = "A", direction = -1, trans = "log10") +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.7, 0.2)
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      direction = "horizontal", barwidth = 15, barheight = 1,
      frame.colour = "black", ticks = FALSE, title.position = "top"
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

# ------------------------------------------------------------------------------
# Show donors

# Select colors
color.rep <- c(
  "A"       = "skyblue",
  "B"       = "indianred",
  "Doublet" = "grey"
)

# Select data
data <- tidyr::as_tibble(object@reductions$umap@cell.embeddings)
names(data) <- c("x", "y")
data$col <- object@meta.data$donor
data$shape <- stringr::str_split(colnames(object), "-", simplify = TRUE)[, 2]

# Plot
ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = col,
    shape = shape,
    label = col
  )
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = color.rep) +
  ggplot2::theme_void(base_size = 30) +
  ggplot2::theme(
    legend.position = c(0.6, 0.4)
  ) +
  ggplot2::labs(col = "Donor", shape = "Replicate") +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 5), order = 1),
    shape = ggplot2::guide_legend(override.aes = list(size = 5))
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

# ------------------------------------------------------------------------------
# Marker gene dotplot

genes <- list(
  "Control"    = c(
    "LYZ", "LCP1"
  ),
  "SARS-CoV-2" = c(
    "TGFBI", "GLUL", "MMP9", "CD84", "S100A6", "CAPG", "MRC1",
    "CD163", "LGMN", "CD9", "MERTK", "SLAMF8", "CMKLR1", "NRP1", "ADAP2",
    "SPRED1", "CALR"
  ),
  "3p-hpRNA"   = c(
    "IFIT1", "IFITM3", "IFIT3", "ISG15", "TNFSF10", "IFI6", "IFI44L"
  ),
  "R848"       = c(
    "CCL3", "IL1B", "CCL3L1", "CCL4", "IL6", "POU2F2", "CCL4L2"
  )
)
genes <- as.character(unlist(genes))

ids <- features$ENSEMBL[match(genes, features$SYMBOL)]
ids <- ids[ids %in% rownames(object@assays$RNA@data) & !duplicated(ids)]

# Tidy data for plotting
data <- tidyr::as_tibble(t(as.matrix(object@assays$RNA@data[ids, ])))
colnames(data) <- features$SYMBOL[match(ids, features$ENSEMBL)]

data$Cluster <- object@meta.data$condition

data <- tidyr::gather(data, "Gene", "Expression", -Cluster)
data$Gene <- factor(data$Gene, levels = unique(data$Gene))
data$N <- rep(1, length(data$Cluster))
data$Expressed <- data$Expression > 0

data <- dplyr::group_by(data, Cluster, Gene)
data <- dplyr::summarize(
  .data        = data,
  "Expression" = mean(Expression),
  "Cells"      = sum(N),
  "Expressed"  = sum(Expressed)
)
data$PCT <- data$Expressed / data$Cells * 100

data <- dplyr::ungroup(data)
data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(
  data, zscore = scale(Expression)
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
  ggplot2::theme_classic(base_size = 20) +
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

# ------------------------------------------------------------------------------
# Show normalized gene expression

# Select data
data <- tidyr::as_tibble(object@reductions$umap@cell.embeddings)
names(data) <- c("x", "y")
for (gene in genes) {
  id <- features$ENSEMBL[match(gene, features$SYMBOL)]
  if (id %in% rownames(object)) {
    data[[gene]] <- object@assays$SCT@data[id, ]
  }
}

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
  ggplot2::theme_void(base_size = 20) +
  ggplot2::theme(
    legend.position = "top"
  ) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barheight = 1, barwidth = 20, frame.colour = "black", ticks = FALSE
    )
  ) +
  ggplot2::labs(col = NULL) +
  ggplot2::coord_fixed() + 
  ggplot2::annotate(
    geom  = "segment", 
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
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.05, "inches"), type   = "closed"
    )
  )

# ------------------------------------------------------------------------------
# Differential expression between stimulation conditions

# Calculate DE genes
markers <- scran::findMarkers(
  x         = object@assays$SCT@data,
  groups    = object@meta.data$condition,
  pval.type = "some",
  test.type = "wilcox",
  direction = "up", 
  block     = object@meta.data$donor
)
for (i in names(markers)) {
  markers[[i]] <- as.data.frame(markers[[i]])
  markers[[i]]$cluster <- i
  markers[[i]]$ENSEMBL <- row.names(markers[[i]])
}
markers <- dplyr::bind_rows(as.list(markers))
markers$gene <- features$SYMBOL[match(markers$ENSEMBL, features$ENSEMBL)]
markers <- markers[, order(names(markers), decreasing = TRUE)]

cutoff <- 1e-15
de <- markers[markers$FDR < cutoff, ]

ids <- features$ENSEMBL[match(de$gene, features$SYMBOL)]
data <- object@assays$SCT@data[
  ids, order(object@meta.data$condition, object@meta.data$donor)
  ]
data <- t(scale(t(as.matrix(data))))

cann <- data.frame(
  Condition = object@meta.data$condition,
  Donor     = object@meta.data$donor,
  row.names = colnames(object)
)
color.cann <- list(
  Condition = color.condition,
  Donor     = color.rep
)

rann <- data.frame(
  Cluster = de$cluster
)
limits <- c(-2, 2)
breaks <- seq(from = min(limits), to = max(limits), length.out = 100)

pheatmap::pheatmap(
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
  gaps_col             = head(as.numeric(cumsum(table(cann$Condition))), -1),
  gaps_row             = head(as.numeric(
    cumsum(table(rann$Cluster)[unique(rann$Cluster)])
    ), -1), 
  cellwidth            = 0.4
)

# ------------------------------------------------------------------------------
# Transcription factor enrichment by ChEA3

genes <- split(
  de$gene, 
  factor(de$cluster, levels = c("Control", "SARS-CoV-2", "3p-hpRNA", "R848"))
  )

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

result$Rank <- as.numeric(result$Rank)
result$Score <- as.numeric(result$Score)
result$map <- unlist(
  lapply(stringr::str_split(result$Overlapping_Genes, ","), length)
)
result$geneRatio <- result$map / result$bg
result$cluster <- factor(result$cluster, unique(result$cluster))

# Order result by Score
result <- result[order(result$cluster, result$Score), ]
tfs <- character()

tfs <- result$TF[result$Score < 35]

data <- result[which(result$TF %in% tfs), ]
data$TF <- factor(data$TF, levels = unique(tfs))

ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x = TF,
    y = cluster,
    fill = Score
  )
) +
  ggplot2::geom_tile(col = "white") +
  viridis::scale_fill_viridis(option = "D", direction = 1, trans = "log10") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    aspect.ratio = .2,
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

# ------------------------------------------------------------------------------
# Enrichment of fibrosis gene sets

# Retrieve gene set dictionary
con <- "https://nubes.helmholtz-berlin.de/s/25ZGCXAHdBffZCP/download"
file <- "~/Downloads/genesets.xlsx"
download.file(url = con, destfile = file)
dictionary <- readxl::read_excel(file)
names(dictionary) <- c("ref", "term", "disease", "gene")

key <- "IPF"
n <- 50

dict <- list()
genesets <- list()
for (i in unique(dictionary$term[dictionary$disease == key])) {
  dict[[i]] <- head(dictionary[dictionary$term == i, ], n)
  genesets[[i]] <- as.character(head(dictionary$gene[dictionary$term == i], n))
}
dict <- dplyr::bind_rows(dict)

# Overrepresentation analysis ==================================================
# Contingency table
matrix(
  c("A", "B", "C", "D"), nrow = 2,
  dimnames = list(
    c("DE", "Not.DE"), c("In.gene.set", "Not.in.gene.set"))
)

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

term2ref <- unique(dict[, c("term", "ref")])$ref
names(term2ref) <- unique(dict[, c("term", "ref")])$term

data$Geneset <- factor(
  data$Geneset, unique(data$Geneset)[c(4,2,10,8,7,3,11,1,5,6,9,12)]
)

cols <- c(
  "Morse et al."   = "seagreen",
  "Adams et al."   = "navy",
  "Reyfman et al." = "dodgerblue",
  "Ayaub et al."   = "purple"
)
ref.color <- cols[term2ref[levels(data$Geneset)]]

ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x    = Geneset,
    y    = Cluster,
    fill = -log10(p.adjust)
  )
) +
  ggplot2::geom_tile(color = "white", size = 1) +
  ggplot2::scale_fill_gradient2(low = "white", high = "red") +
  ggplot2::scale_size_area(max_size = 12) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    aspect.ratio = 0.4,
    legend.position = "top",
    axis.text.x     = ggplot2::element_text(
      angle = 40, hjust = 1, vjust = 0.9
    ),
    axis.ticks.x = ggplot2::element_line(
      size = 7, color = ref.color, arrow = grid::arrow(
        length = grid::unit(8, "pt"), angle = 90
      )
    ),
    axis.title.y = ggplot2::element_text(angle = 0, color = "white"),
    legend.text  = ggplot2::element_text(size = 15),
    legend.title = ggplot2::element_text(size = 15)
  ) +
  ggplot2::labs(x = NULL, y = NULL, size = "% overlap") +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      barwidth = 34, barheight = 1, order = 1, title.vjust = 1,
      frame.colour = "black", ticks = FALSE,
      title.position = "top", title.hjust = 0.5
    ),
    size = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5
    )
  ) +
  ggplot2::annotate(
    geom = "text",
    x = c(2,5,8,11),
    y = rep(0.1, 4),
    label = unique(names(ref.color)),
    color = unique(ref.color),
    size  = 6
  ) +
  ggplot2::expand_limits(y = -0.5)

# Module scores ================================================================
idsets <- list()
for (i in names(genesets)) {
  print(i)
  ids <- as.character(na.omit(
    features$ENSEMBL[match(genesets[[i]], features$SYMBOL)]
  ))
  ids <- ids[ids %in% rownames(object)]
  
  object <- Seurat::AddModuleScore(
    object   = object,
    features = list(ids),
    name     = "geneset",
    seed     = 1993,
    nbin     = round(length(ids)/2)
  )
  names(object@meta.data)[length(names(object@meta.data))] <- i
}

scores <- object@meta.data[, names(object@meta.data) %in% dict$term]

scores$x <- object@reductions$umap@cell.embeddings[, 1]
scores$y <- object@reductions$umap@cell.embeddings[, 2]

data <- tidyr::gather(scores, "Geneset", "Score", -x, -y)
data$Geneset <- factor(data$Geneset, levels = unique(data$Geneset))

# Set colorscale limits
limits <- c(-0.5, 1)
data$Score[data$Score > max(limits)] <- max(limits)
data$Score[data$Score < min(limits)] <- min(limits)

ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y,
    col = Score
  )
) +
  ggplot2::geom_point(size = 0.75) +
  ggplot2::facet_wrap(~Geneset, nrow = 3)+
  ggplot2::scale_color_distiller(
    palette = "RdYlBu", direction = -1, limits = limits
  ) +
  ggplot2::theme_void(base_size = 20) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 15)
  ) +
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

# Violin Plots
data <- scores
data$Cluster <- object@meta.data$condition

data <- tidyr::gather(data, "Geneset", "Score", -Cluster, -x, -y)

data$Geneset <- factor(data$Geneset, unique(data$Geneset))
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

# Plot
ps$label[ps$label == Inf] <- max(ps$label[ps$label < Inf])

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
  ggplot2::scale_fill_manual(values = color.condition) +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(
    strip.text   = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) + 
  ggplot2::scale_y_continuous(limits = c(-0.5, 2.4)) +
  ggplot2::guides(
    col = ggplot2::guide_colorbar(
      barwidth = 15, barheight = 1, frame.colour = "black", ticks = FALSE,
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
      label = label, y = 2.2, fill = label, col = NULL
    ),
    size = 5,
    col = "white"
  ) +
  viridis::scale_fill_viridis(direction = -1, option = "A")

# end of document
################################################################################