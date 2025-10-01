#' Scanorama method
#'
#' Scanorama is a Mutual Nearest Neighbors method. This an algorithm that
#' identifies and merges the shared cell types among all pairs of datasets and
#' accurately integrates heterogeneous collections of scRNA-seq data.
#'
#' @param input List that contains the expression matrices for each batches and
#'  the genes.
#' @param return_dimred A logical to returning integrated low-dimesional
#' embeddings.
#' @param batch_size The batch size used in the alignment vector computation.
#' @param verbose Print progress bars and output.
#' @param ds_names Reports data set names in logging output
#' @param dimred Number of dimensions of corrected embedding.
#' @param approx A logical to use approximate nearest neighbors
#' @param sigma Correction smoothing parameter on Gaussian kernel.
#' @param alpha Alignment score minimum cutoff.
#' @param knn Number of nearest neighbors to use for matching.
#' @param return_dense A logical to return dense matrices.
#' @param hvg Use this number of top highly variable genes based on dispersion.
#' @param union A logical to consider the union of the genes.
#' In alternative, the intersection is considered.
#' @param seed Random seed to use.
#'
#' @export
#' @importFrom basilisk basiliskRun basiliskStart basiliskStop
#' @importFrom reticulate import
#'
#' @return A list that contains the corrected gene expression matrix and
#' the corrected low-dimensional reduction
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' sim <- lapply(unique(SingleCellExperiment::colData(sim)[, "Batch"]),
#'               function(i) sim[, SingleCellExperiment::colData(sim)[, "Batch"] == i])
#' ll <- c(lapply(sim, function(x) t(as.matrix(SingleCellExperiment::logcounts(x)))),
#'         lapply(sim, function(x) rownames(x)))
#' scanorama <- scanoramaRun(input = ll, return_dimred = FALSE)
#'
scanoramaRun <- function(input, return_dimred = FALSE,
                         batch_size = 5000, verbose = TRUE, ds_names = NULL,
                         dimred = 100, approx = TRUE, sigma = 15, alpha = 0.10,
                         knn = 20, return_dense = FALSE, hvg = NULL,
                         union = FALSE, seed = 0) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, return_dimred,
                                                 batch_size, verbose, ds_names,
                                                 dimred, approx, sigma, alpha,
                                                 knn, return_dense, hvg, union,
                                                 seed) {
    scanorama <- import("scanorama")

    n <- length(input)
    n_batch <- n/2

    arg <- c(list(datasets_full = input[1:n_batch],
                  genes_list = input[(n_batch + 1):n],
                  return_dimred = return_dimred,
                  batch_size = as.integer(batch_size), verbose = verbose,
                  ds_names = ds_names, dimred = as.integer(dimred),
                  approx = approx, sigma = sigma, alpha = alpha,
                  knn = as.integer(knn), return_dense = return_dense,
                  hvg = hvg, union = union, seed = as.integer(seed)))
    do.call(scanorama$correct, arg)
  }, input = input, return_dimred = return_dimred, batch_size = batch_size,
  verbose = verbose, ds_names = ds_names, dimred = dimred, approx = approx,
  sigma = sigma, alpha = alpha, knn = knn, return_dense = return_dense,
  hvg = hvg, union = union, seed = seed)
}
