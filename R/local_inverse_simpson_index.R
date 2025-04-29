#' Local Inverse Simposon Index
#'
#' @param input input
#' @param label_true label_true
#' @param reduction reduction
#' @param meta_data meta_data
#' @param perplexity perplexity
#' @param nn_eps nn_eps
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

