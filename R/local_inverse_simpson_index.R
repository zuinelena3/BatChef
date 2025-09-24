#' Local Inverse Simposon Index (LISI)
#'
#' Local Inverse Simpson’s Index (LISI) is a local level metric based on
#' the kNN algorithm. This metric defines the effective number of datasets in
#' a neighborhood. In other words, LISI determines the number of neighbor cells
#' necessary before one batch is observed twice.
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#' @param meta_data A data frame with one row per cell.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with `RANN:nn2()`.
#'
#' @export
#' @importFrom lisi compute_lisi
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats median
#'
#' @return A numeric value.
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' lisi <- local_inverse_simpson_index(input = sim, label_true = "Group",
#'                                     reduction = "PCA")
#'
local_inverse_simpson_index <- function(input, label_true, reduction,
                                        meta_data = colData(input),
                                        perplexity = 30, nn_eps = 0) {
  red <- reducedDim(input, reduction)
  lisi <- median(compute_lisi(X = red, meta_data = meta_data,
                              label_colnames = label_true,
                              perplexity = perplexity, nn_eps = nn_eps)[, 1])
}

