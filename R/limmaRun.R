#' limma method
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
#' @param batch A string specifying the batch variable.
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom limma removeBatchEffect
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData
#'
#' @return A corrected matrix
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' limma <- limmaRun(input = sim, batch = "Batch")
#'
limmaRun <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(list(x = assay(input, assay_type), batch = colData(input)[, batch]),
            capture_params(removeBatchEffect, list(...)))
  do.call(removeBatchEffect, args)
}
