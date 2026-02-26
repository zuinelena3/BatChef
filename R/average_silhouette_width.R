#' Average Silhouette Width (ASW)
#'
#' The average silhouette width for each cell is defined as the difference
#' between the average of  within-cluster distances to all cells and the average
#' between-cluster distances of that cell to the closest cluster divided
#' by their maximum.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param label_true A string specifying the ground truth label.
#' @param reduction A string specifying the dimensional reduction.
#' @param metric The metric to use when calculating distance between instances.
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @return A numeric value
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' asw <- average_silhouette_width(
#'   input = sim, label_true = "Group",
#'   reduction = "PCA"
#' )
#'
average_silhouette_width <- function(input, label_true, reduction,
                                     metric = "euclidean") {
  red <- reducedDim(input, reduction)
  dist_matrix <- dist(red, method = metric)

  group <- as.factor(colData(input)[, label_true])
  levels(group) <- seq_along(levels(group))

  group <- as.numeric(group)

  sil <- silhouette(x = group, dist = dist_matrix)
  asw <- mean(sil[, "sil_width"])
}
