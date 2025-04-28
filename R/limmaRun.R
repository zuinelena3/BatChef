#' limma method
#'
#' @param input A SingleCellExperiment object.
#' @param batch A string specifying the batch for each cell.
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom limma removeBatchEffect
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData
#'
limmaRun <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(list(x = assay(input, assay_type), batch = colData(input)[, batch]), capture_params(removeBatchEffect, list(...)))
  do.call(removeBatchEffect, args)
}
