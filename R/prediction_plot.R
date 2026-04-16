#' A shifted log transformation
#'
#' @param mat A matrix of data characteristics
#' @param base logarithm base
#'
#' @returns A shifted log transformation matrix.
#'
log_transf <- function(mat, base = exp(1)) {
  apply(mat, 2, function(x) {
    vapply(x, function(val) {
      if (val < 0) {
        -log((-val) + 1, base)
      } else if (val == 0) {
        0
      } else {
        log(val + 1, base)
      }
    }, FUN.VALUE = numeric(1))
  })
}

#' Prediction plot
#'
#' @param params A data.frame of data characteristics
#'
#' @returns A ggplot object
#'
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_point geom_sf
#' @importFrom ggplot2 theme_classic theme xlab ylab labs element_text
#' @importFrom grDevices chull
#' @importFrom stats predict
#' @importFrom sf st_as_sf st_sfc st_intersects st_sf st_polygon st_voronoi
#' @importFrom sf st_union st_collection_extract st_intersection st_join st_nearest_feature
#'
prediction_plot <- function(params) {
  file_path <- system.file("extdata", "svm.rda", package = "BatChef")
  load(file_path)

  params <- as.data.frame(t(as.matrix(log_transf(params))))

  new_coords <- predict(res$pca, params)
  new_coords <- data.frame(PC1 = new_coords[, 1], PC2 = new_coords[, 2])

  coords <- res$coords
  hull <- chull(coords[, 1], coords[, 2])
  hull <- c(hull, hull[1])
  polygon_coords <- coords[hull, ]
  poly <- polygon_coords[-11, ]

  points_sf <- st_as_sf(coords, coords = c("PC1", "PC2"))
  points <- st_as_sf(new_coords, coords = c("PC1", "PC2"), crs = NA)
  polygon_sf <- st_sfc(st_polygon(list(as.matrix(polygon_coords[, 1:2]))))
  polygon_sf <- st_sf(polygon_sf)
  inside <- st_intersects(points, polygon_sf, sparse = FALSE)[, 1]

  if (!inside) {
    message("WARNING: data point is outside the boundaries!")
  }

  voronoi <- st_voronoi(st_union(points_sf))
  vor_sf <- st_collection_extract(voronoi)
  vor_clip <- st_intersection(st_as_sf(vor_sf), polygon_sf)

  vor_colored <- st_join(vor_clip, points_sf, join = st_nearest_feature)

  p <- ggplot() +
    geom_sf(data = vor_colored, aes(fill = svm_best), color = "black") +
    geom_sf(data = polygon_sf, fill = NA, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = res$color) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(fill = "Method") + xlab("PC1") + ylab("PC2")

  p <- p + geom_point(
    data = new_coords, aes(x = PC1, y = PC2),
    fill = "red", color = "red", size = 6
  )
  return(p)
}
