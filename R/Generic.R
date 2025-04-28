#' batchCorrect
#'
#' Generic function for defining batch effect correction methods.
#'
#' @param input input
#' @param batch batch
#' @param params params
#'
#' @export
#' @import methods
#
setGeneric("batchCorrect", function(input, batch, params)
  standardGeneric("batchCorrect"), signature = "params")
