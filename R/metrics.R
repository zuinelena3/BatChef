#' Performance evaluation metrics.
#'
#' Performance evaluation metrics: Wasserstein distance, Local Inverse Simpson's
#' Index, Average Silhouette Width, Adjusted Rand Index, and
#' Normalized Mutual Information.
#'
#' This function performs Leiden clustering before computing the evaluation
#' metrics.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param group A string specifying the ground truth labels.
#' @param reduction A string specifying the dimensional reduction.
#' on which the clustering analysis will be performed.
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed (Default: TRUE).
#' @param resolution A numeric value specifying the resolution parameter.
#' @param k An integer scalar specifying the number of nearest neighbors.
#' @param n_iter Number of iterations of the Leiden clustering algorithm to perform.
#' Positive values above 2 define the total number of iterations to perform,
#' -1 has the algorithm run until it reaches its optimal clustering
#' @param rep Number of times the Wasserstein distance is calculated.
#' @param mc_cores The number of cores to use.
#' @param variant How to compute the normalizer in the denominator.
#' @param metric The metric to use when calculating distance between instances.
#' @param meta_data A data frame with one row per cell.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with `RANN:nn2()`.
#'
#' @export
#'
#' @return A data.frame object
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 110),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' metrics <- metrics(
#'   input = sim, batch = "Batch", group = "Group",
#'   reduction = "PCA", rep = 5, n_iter = 2
#' )
#'
metrics <- function(input, batch, group, reduction,
                    rep = 10, mc_cores = 1, nmi_compute = TRUE, resolution = NULL,
                    k = 10, n_iter = -1, metric = "euclidean", variant = "sum",
                    meta_data = colData(input), perplexity = 30, nn_eps = 0) {
  gr <- colData(input)[, group]
  clust <- leiden_clustering(
    input = input, label_true = group,
    reduction = reduction, nmi_compute = nmi_compute,
    resolution = resolution, k = k, store = FALSE,
    n_iter = n_iter
  )

  ari <- adjustedRandIndex(x = clust, y = gr)

  nmi <- NMI(c1 = as.vector(gr), c2 = as.vector(clust), variant = variant)

  casw <- average_silhouette_width(
    input = input, label_true = group,
    reduction = reduction, metric = metric
  )

  iasw <- average_silhouette_width(
    input = input, label_true = batch,
    reduction = reduction, metric = metric
  )

  clisi <- local_inverse_simpson_index(
    input = input, label_true = group,
    reduction = reduction,
    meta_data = meta_data,
    perplexity = perplexity, nn_eps = nn_eps
  )

  ilisi <- local_inverse_simpson_index(
    input = input, label_true = batch,
    reduction = reduction, meta_data = meta_data,
    perplexity = perplexity, nn_eps = nn_eps
  )

  wass <- wasserstein_distance(
    input = input, batch = batch, reduction = reduction,
    rep = rep, mc_cores = mc_cores
  )
  wass <- mean(wass$wasserstein)

  return(data.frame(
    method = reduction, wasserstein = wass, iasw, ilisi,
    ari = ari, nmi = nmi, casw, clisi
  ))
}
