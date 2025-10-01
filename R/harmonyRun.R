#' Harmony method
#'
#' Harmony is a mixture-model based method.
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
#' @param batch A string specifying the batch variable.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom harmony RunHarmony
#'
#' @return A \linkS4class{SingleCellExperiment} object which contains a
#' corrected low-dimensional space.
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' harmony <- harmonyRun(input = sim, batch = "Batch")
#'
harmonyRun <- function(input, batch, ...) {
  args <- c(list(object = input, group.by.vars = batch),
            capture_params(RunHarmony, list(...)))
  do.call(RunHarmony, args)
}
