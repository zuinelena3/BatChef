#' Adjusted Rand Index (ARI)
#'
#' The Adjusted Rand Index (ARI) is a metric used to measure the similarity
#' between the clustering and ground truth labels.
#'
#' After computing Leiden clustering algorithm, the ARI metric is computed.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param label_true A string specifying the ground truth label.
#' @param reduction A string specifying the dimensional reduction
#' on which the clustering analysis will be performed.
#' @param resolution A numeric value specifying the resolution parameter.
#' @param k An integer scalar specifying the number of nearest neighbors.
#' @param n_iter Number of iterations of the Leiden clustering algorithm.
#' Positive values above 2 define the total number of iterations to perform,
#' -1 has the algorithm run until it reaches its optimal clustering.
#'
#' @export
#' @importFrom mclust adjustedRandIndex
#'
#' @return A numeric value
#' @examples
#' sim <- simulate_data(
#'   n_genes = 500, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' ari <- adjusted_rand_index(
#'   input = sim, label_true = "Group", reduction = "PCA",
#'   resolution = 0.5, n_iter = 2
#' )
#'
adjusted_rand_index <- function(input, label_true, reduction,
                                resolution, k = 15, n_iter = -1) {
  group <- colData(input)[, label_true]

  clust <- leiden_clustering(
    input = input, label_true = label_true,
    reduction = reduction, nmi_compute = FALSE,
    resolution = resolution, k = k, store = FALSE,
    n_iter = n_iter
  )
  ari <- adjustedRandIndex(x = clust, y = group)
}
