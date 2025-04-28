#' ComBat method
#'
#' @param input A SingleCellExperiment object.
#' @param batch A string specifying the batch for each cell.
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom sva ComBat_seq
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment counts colData
#'
combatRun <- function(input, batch, assay_type = "counts", ...) {
  args <- c(list(counts = as.matrix(assay(input, assay_type)), batch = colData(input)[, batch]), capture_params(ComBat_seq, list(...)))
  do.call(ComBat_seq, args)
}
