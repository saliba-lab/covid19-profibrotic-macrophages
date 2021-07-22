################################################################################
# Shiny dashboard - Helper functions

#' Plot low dimensional embedding
#' 
#' @object    The Seurat object to pull data from
#' @embedding Name of the low-dimensional embedding
#' @coldata   Gene name or metadata slot
#' @brush     Brushed points (xmin, xmax, ymin, ymax)
#' @pt.size   Pointsize
#' 
plot.embedding <- function(
  object,
  embedding = tail(names(object@reductions), 1),
  coldata = "",
  brush = NULL,
  pt.size   = 1
) {
  # Fetch data from object
  data <- data.frame(
    x = object@reductions[[embedding]]@cell.embeddings[, 1],
    y = object@reductions[[embedding]]@cell.embeddings[, 2],
    row.names = colnames(object)
  )
  # Add metadata/gene expression
  id <- convertFeatures(coldata, object@misc$features, 2, 1)
  if (coldata %in% names(object@meta.data)) {
    data$col <- object@meta.data[[coldata]]
  } else {
    if (id %in% rownames(object)) {
      data$col <- object@assays$RNA@data[id, ]
    } else {
      data$col <- NaN
    }
  }
  
  # Create guide elements
  if (class(data$col) %in% c("numeric", "integer")) {
    color_guide <- ggplot2::guide_colorbar(
      barwidth = 1, barheight = 15, ticks = FALSE, frame.colour = "black"
    )
    ann_cols <- viridis::scale_color_viridis(option = "B", direction = -1)
  } else {
    color_guide <- ggplot2::guide_legend(
      override.aes = list(size = 8)
    )
    ann_cols <- NULL
  }
  
  # Subset data based on cells
  if (is.null(brush)) {
    xlim <- NULL
    ylim <- NULL
  } else {
    xlim <- ggplot2::xlim(brush$xmin, brush$xmax)
    ylim <- ggplot2::ylim(brush$ymin, brush$ymax)
  }
  
  # Order data
  data <- data[order(data$col), ]
  
  # Plot
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, col = col)) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 25) +
    ggplot2::labs(col = NULL, title = coldata) +
    ggplot2::guides(
      color = color_guide
    ) +
    ann_cols +
    xlim +
    ylim
}

#' Convert between gene symbols and Ensembl IDs
#' 
#' @features Character vector of gene symbols or Ensembl IDs
#' @convtab  Table with gene identifiers to be converted
#' @from     Column with current identifiers (matching features)
#' @to       Column with desired identifiers
#' 
convertFeatures <- function(
  features = NULL,
  convtab  = NULL,
  from     = 2,
  to       = 1
) {
  return(convtab[[to]][match(features, convtab[[from]])])
}

# end of document
################################################################################