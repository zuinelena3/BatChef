#' Normalized Mutual Information (NMI)
#'
#' The Normalized Mutual Information (NMI) metric is used to compare
#' the overlap between the true cell-type and clustering labels computed after
#' batch correction.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param label_true A string specifying the ground truth labels.
#' @param reduction A string specifying the dimensional reduction
#' on which the clustering analysis will be performed.
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed (Default: FALSE).
#' @param resolution A numeric value specifying the resolution parameter.
#' @param k An integer scalar specifying the number of nearest neighbors.
#' @param variant How to compute the normalizer in the denominator.
#'
#' @export
#' @importFrom aricode NMI
#'
#' @return A numeric value.
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' nmi <- normalized_mutual_info(
#'   input = sim, label_true = "Group", reduction = "PCA",
#'   nmi_compute = FALSE, resolution = 0.5
#' )
#'
normalized_mutual_info <- function(input, label_true, reduction, nmi_compute = FALSE,
                                   resolution = 1, k = 10, variant = "sum") {
  clust <- leiden_clustering(
    input = input, label_true = label_true,
    reduction = reduction, nmi_compute = nmi_compute,
    resolution = resolution, k = k, store = FALSE
  )

  group <- colData(input)[, label_true]

  nmi <- NMI(c1 = as.vector(group), c2 = as.vector(clust), variant = variant)
}
