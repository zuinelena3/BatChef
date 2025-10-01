#' LIGER method
#'
#' LIGER is a integrative non-negative matrix factorization method
#'
#' @param input A \linkS4class{liger} object.
#' @param method A string specifyng the batch correction method. Choose from
#' \code{"iNMF"}, \code{"onlineINMF"}, \code{"UINMF"}. Default \code{"iNMF"}.
#' @param quantiles Number of quantiles to use for quantile normalization.
#' @param reference Character, numeric or logical selection of one dataset,
#' out of all available datasets in object, to use as a "reference" for
#' quantile normalization.
#' @param min_cells Minimum number of cells to consider a cluster shared across
#' datasets
#' @param n_neighbors Number of nearest neighbors for within-dataset knn graph.
#' @param use_dims Indices of factors to use for shared nearest factor
#' determination.
#' @param center A logical to center data.
#' @param max_sample Maximum number of cells used for quantile normalization of
#' each cluster and factor
#' @param eps The error bound of the nearest neighbor search.
#' @param refine_knn A logical to increase robustness of cluster assignments
#'  using KNN graph.
#' @param cluster_name Variable name that will store the clustering result
#' in metadata of a liger object or a Seurat object.
#' @param seed Random seed
#' @param verbose Print progress bar/messages
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom rliger runIntegration quantileNorm
#'
#' @returns A \linkS4class{liger} object, which contains the quantile align
#' factor loadings.
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                     compute_pca = TRUE, output_format = "SingleCellExperiment")
#'
#' ll <- lapply(unique(SingleCellExperiment::colData(sim)[, "Batch"]),
#'              function(i) sim[, SingleCellExperiment::colData(sim)[, "Batch"] == i])
#' names(ll) <- unique(SingleCellExperiment::colData(sim)[, "Batch"])
#' # Create a liger object
#' lo <- rliger::createLiger(rawData = ll)
#' lo <- rliger::normalize(object = lo)
#' lo <- rliger::scaleNotCenter(object = lo, features = rownames(sim))
#' liger <- ligerRun(input = lo)
#'
ligerRun <- function(input, method = "iNMF", quantiles = 50, reference = NULL,
                     min_cells = 20, n_neighbors = 20,
                     use_dims = NULL, center = FALSE, max_sample = 1000, eps = 0.9,
                     refine_knn = TRUE, cluster_name = "quantileNorm_cluster",
                     seed = 1, verbose = FALSE, ...) {
  args <- c(list(object = input, method = method),
            capture_params(runIntegration, list(...)))
  out <- do.call(runIntegration, args)

  out <- quantileNorm(object = out, quantiles = quantiles, reference = reference,
                      minCells = min_cells, nNeighbors = n_neighbors,
                      useDims = use_dims, center = center, maxSample = max_sample,
                      eps = eps, refineKNN = refine_knn, clusterName = cluster_name,
                      seed = seed, verbose = verbose)
}
