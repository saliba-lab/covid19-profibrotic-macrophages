qc.plot <- function(object) {
  data <- data.frame(
    libsize = object@meta.data$nCount_RNA,
    genes = object$nFeature_RNA,
    percent.mito = object$percent.mito,
    donor = object$donor
  )
  data <- tidyr::gather(data, "key", "value", -libsize, -donor)
  ggplot2::ggplot(data, ggplot2::aes(value, libsize, col = donor)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::facet_wrap(~key, scales = "free") +
    ggplot2::theme_classic(20) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    ) +
    ggplot2::labs(x = NULL, y = "Library size")
}
