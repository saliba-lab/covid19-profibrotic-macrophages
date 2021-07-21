################################################################################
# Shiny dashboard - Helper functions

#' Plot low dimensional embedding
#' 
#' @object The Seurat object to pull data from
#' @embedding Name of the low-dimensional embedding
#' @coldata Gene name or metadata slot
#' 
plot.embedding <- function(
  object,
  embedding = tail(names(object@reductions), 1),
  coldata = tail(names(object@meta.data), 1)
) {
  # Fetch data from object
  data <- data.frame(
    x = object@reductions[[embedding]]@cell.embeddings[, 1],
    y = object@reductions[[embedding]]@cell.embeddings[, 2]
  )
  if (coldata %in% names(object@meta.data)) {
    data$col <- object@meta.data[[coldata]]
  } else {
    if (coldata %in% rownames(ds)) {
      data$col <- object@assays$RNA@data[coldata, ]
    } else {
      data$col <- NULL
    }
  }
  
  # Plot
  ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, col = col)) +
    ggplot2::geom_point() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 20)
}

# end of document
################################################################################