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
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed.
#' @param resolution A numeric value specifying the resolution parameter.
#' @param k An integer scalar specifying the number of nearest neighbors.
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
#'   nmi_compute = FALSE, resolution = 0.5
#' )
#'
adjusted_rand_index <- function(input, label_true, reduction, nmi_compute = FALSE,
                                resolution, k = 10) {
  group <- colData(input)[, label_true]
  clust <- leiden_clustering(
    input = input, label_true = label_true,
    reduction = reduction, nmi_compute = nmi_compute,
    resolution = resolution, k = k, store = FALSE
  )

  ari <- adjustedRandIndex(x = clust, y = group)
}
