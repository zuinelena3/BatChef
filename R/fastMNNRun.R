#' fastMNN method
#'
#' @param input input
#' @param batch batch
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment colData
fastMNNRun <- function(input, batch, ...) {
  args <- c(list(counts = input, batch = colData(input)[, batch]), capture_params(fastMNN, list(...)))
  do.call(fastMNN, args)
}
