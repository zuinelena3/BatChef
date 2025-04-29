#' scMerge2 method
#'
#' @param input input
#' @param batch batch
#' @param assay_type assay_type
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom scMerge scMerge2
#' @importFrom SummarizedExperiment assay colData
scMerge2Run <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(list(exprsMat = assay(input, assay_type), batch = colData(input)[, batch]), capture_params(scMerge2, list(...)))
  do.call(scMerge2, args)
}
