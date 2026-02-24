#' limma method
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param batch A string specifying the batch variable.
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @return A corrected matrix
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = FALSE,
#'   output_format = "SingleCellExperiment"
#' )
#' limma <- limmaRun(input = sim, batch = "Batch")
#'
#' @importFrom limma removeBatchEffect
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData
#'
limmaRun <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(
    list(
      x = as.matrix(assay(input, assay_type)),
      batch = colData(input)[, batch],
      design = matrix(1, ncol(input), 1)
    ),
    capture_params(removeBatchEffect, list(...))
  )
  do.call(removeBatchEffect, args)
}
