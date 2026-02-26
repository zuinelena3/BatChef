#' A shifted log transformation
#'
#' @param mat A matrix of data characteristics
#' @param base logarithm base
#'
#' @returns A shifted log transformation matrix.
#'
log_transf <- function(mat, base = exp(1)) {
  apply(mat, 2, function(x) {
    sapply(x, function(val) {
      if (val < 0) {
        -log((-val) + 1, base)
      } else if (val == 0) {
        0
      } else {
        log(val + 1, base)
      }
    })
  })
}

#' Prediction plot
#'
#' @param params A data.frame of data characteristics
#'
#' @returns A ggplot object
#'
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_point
#' @importFrom ggplot2 theme_classic theme labs element_text
#' @importFrom ggvoronoi geom_voronoi
#' @importFrom grDevices chull
#' @importFrom stats predict
#'
predition_plot <- function(params) {
  files <- c("pca.rda", "coords.rda", "palette.rda")
  file_paths <- system.file("extdata", files, package = "BatChef")
  lapply(file_paths, load, envir = .GlobalEnv)

  params <- as.data.frame(t(as.matrix(log_transf(params))))

  new_coords <- predict(pca, params)
  new_coords <- data.frame(PC1 = new_coords[, 1], PC2 = new_coords[, 2])

  hull <- chull(coords[, 1], coords[, 2])
  hull <- c(hull, hull[1])
  polygon_coords <- coords[hull, ]
  poly <- polygon_coords[-11, ]

  p <- ggplot(coords, aes(x = PC1, y = PC2, fill = svm_best)) +
    geom_voronoi(outline = poly, color = 1, linewidth = 0.1) +
    scale_fill_manual(values = pltt) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(fill = "Method")

  p <- p + geom_point(
    data = new_coords, aes(x = PC1, y = PC2),
    fill = "red", color = "red", size = 6
  )

  suppressWarningsByMsg(
    "deprecated",
    print(p)
  )
}
