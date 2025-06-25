#' Local Inverse Simposon Index
#'
#' @param input A `SingleCellExperiment` object.
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
local_inverse_simpson_index <- function(input, label_true, reduction, meta_data = colData(input), perplexity = 30, nn_eps = 0) {
  red <- reducedDim(input, reduction)
  lisi <- median(compute_lisi(X = red, meta_data = meta_data, label_colnames = label_true, perplexity = perplexity, nn_eps = nn_eps)[, 1])
}

