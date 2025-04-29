#' LIGER method
#'
#' @param input input
#' @param method method
#' @param quantiles quantiles
#' @param reference reference
#' @param minCells minCells
#' @param nNeighbors nNeighbors
#' @param useDims useDims
#' @param center center
#' @param maxSample maxSample
#' @param eps eps
#' @param refineKNN refineKNN
#' @param clusterName clusterName
#' @param seed seed
#' @param verbose verbose
#' @param ... ...
#'
#' @export
#' @importFrom rliger runIntegration quantileNorm
ligerRun <- function(input, method = "iNMF", quantiles = 50, reference = NULL, minCells = 20, nNeighbors = 20,
                     useDims = NULL, center = FALSE, maxSample = 1000, eps = 0.9, refineKNN = TRUE, clusterName = "quantileNorm_cluster",
                     seed = 1, verbose = TRUE, ...) {
  args <- c(list(object = input, method = method), capture_params(runIntegration, list(...)))
  out <- do.call(runIntegration, args)

  out <- quantileNorm(object = out, quantiles = quantiles, reference = reference, minCells = minCells, nNeighbors = nNeighbors,
                      useDims = useDims, center = center, maxSample = maxSample, eps = eps, refineKNN = refineKNN, clusterName = clusterName,
                      seed = seed, verbose = verbose)
}
