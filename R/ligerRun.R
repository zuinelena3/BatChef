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
#' @param minCells Minimum number of cells to consider a cluster shared across
#' datasets
#' @param nNeighbors Number of nearest neighbors for within-dataset knn graph.
#' @param useDims Indices of factors to use for shared nearest factor
#' determination.
#' @param center A logical to center data.
#' @param maxSample Maximum number of cells used for quantile normalization of
#' each cluster and factor
#' @param eps The error bound of the nearest neighbor search.
#' @param refineKNN A logical to increase robustness of cluster assignments
#'  using KNN graph.
#' @param clusterName Variable name that will store the clustering result
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
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' ll <- lapply(unique(colData(sim)[, "Batch"]), function(i) sim[, colData(sim)[, "Batch"] == i])
#' names(ll) <- unique(colData(sim)[, "Batch"])
#' # Create a liger object
#' lo <- createLiger(rawData = ll)
#' lo <- normalize(object = lo)
#' lo <- scaleNotCenter(object = lo, features = rownames(sim))
#' liger <- ligerRun(input = lo)
#'
ligerRun <- function(input, method = "iNMF", quantiles = 50, reference = NULL,
                     minCells = 20, nNeighbors = 20,
                     useDims = NULL, center = FALSE, maxSample = 1000, eps = 0.9,
                     refineKNN = TRUE, clusterName = "quantileNorm_cluster",
                     seed = 1, verbose = TRUE, ...) {
  args <- c(list(object = input, method = method),
            capture_params(runIntegration, list(...)))
  out <- do.call(runIntegration, args)

  out <- quantileNorm(object = out, quantiles = quantiles, reference = reference,
                      minCells = minCells, nNeighbors = nNeighbors,
                      useDims = useDims, center = center, maxSample = maxSample,
                      eps = eps, refineKNN = refineKNN, clusterName = clusterName,
                      seed = seed, verbose = verbose)
}
