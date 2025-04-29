#' scanorama method
#'
#' @param input input
#' @param return_dimred return_dimred
#' @param batch_size batch_size
#' @param verbose verbose
#' @param ds_names ds_names
#' @param dimred dimred
#' @param approx approx
#' @param sigma sigma
#' @param alpha alpha
#' @param knn knn
#' @param return_dense return_dense
#' @param hvg hvg
#' @param union union
#' @param seed seed
#'
#' @export
#' @importFrom basilisk basiliskRun basiliskStart basiliskStop
#' @importFrom reticulate import
scanoramaRun <- function(input, return_dimred = FALSE,
                         batch_size = 5000, verbose = TRUE, ds_names = NULL,
                         dimred = 100, approx = TRUE, sigma = 15, alpha = 0.10, knn = 20,
                         return_dense = FALSE, hvg = NULL, union = FALSE, seed = 0) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, return_dimred, batch_size, verbose, ds_names,
                                                 dimred, approx, sigma, alpha, knn, return_dense, hvg, union, seed) {
    scanorama <- import("scanorama")

    n <- length(input)
    n_batch <- n/2

    arg <- c(list(datasets_full = input[1:n_batch], genes_list = input[(n_batch + 1):n], return_dimred = return_dimred,
                  batch_size = as.integer(batch_size), verbose = verbose, ds_names = ds_names,
                  dimred = as.integer(dimred), approx = approx, sigma = sigma, alpha = alpha, knn = as.integer(knn), return_dense = return_dense,
                  hvg = hvg, union = union, seed = as.integer(seed)))
    do.call(scanorama$correct, arg)
  }, input = input, return_dimred = return_dimred, batch_size = batch_size, verbose = verbose, ds_names = ds_names,
  dimred = dimred, approx = approx, sigma = sigma, alpha = alpha, knn = knn, return_dense = return_dense,
  hvg = hvg, union = union, seed = seed)
}
