#' Harmony method
#'
#' @param input A `SingleCellExperiment` object.
#' @param batch A string specifying the batch variable.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom harmony RunHarmony
#'
harmonyRun <- function(input, batch, ...) {
  args <- c(list(object = input, group.by.vars = batch), capture_params(RunHarmony, list(...)))
  do.call(RunHarmony, args)
}
