#' scMerge2 method
#'
#' @param input A `SingleCellExperiment` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param assay_type A string specifying the assay to use for correction.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom scMerge scMerge2
#' @importFrom SummarizedExperiment assay colData
#'
scMerge2Run <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(list(exprsMat = assay(input, assay_type), batch = colData(input)[, batch]), capture_params(scMerge2, list(...)))
  do.call(scMerge2, args)
}
