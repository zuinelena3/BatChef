#' batchCorrect
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param params A \linkS4class{BatChefParams} object specifying
#' the batch correction method.
#'
#' @export
#' @import methods
#
setGeneric("batchCorrect", function(input, batch, params) {
  standardGeneric("batchCorrect")
}, signature = "params")
