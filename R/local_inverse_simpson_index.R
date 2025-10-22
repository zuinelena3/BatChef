#' Median of Local Inverse Simposon Index (LISI)
#'
#' Local Inverse Simpson’s Index (LISI) is a local level metric based on
#' the kNN algorithm. This metric defines the effective number of datasets in
#' a neighborhood. In other words, LISI determines the number of neighbor cells
#' necessary before one batch is observed twice.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#' @param meta_data A data frame with one row per cell.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with `RANN:nn2()`.
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats median
#'
#' @return A numeric value.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
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
