#' fastMNN method
#'
#' @param input A `SingleCellExperiment` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment colData
#'
fastMNNRun <- function(input, batch, ...) {
  args <- c(list(counts = input, batch = colData(input)[, batch]), capture_params(fastMNN, list(...)))
  do.call(fastMNN, args)
}
