#' Convert into a SingleCellExperiment object
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param batch A string specifying the batch for each cell.
#'
#' @import methods
#' @rdname sceInput
#'
setGeneric("sceInput", function(input, batch)
  standardGeneric("sceInput"), signature = c("input"))

#' @rdname sceInput
#' @import methods
#' @importFrom Seurat as.SingleCellExperiment
#' @aliases sceInput,Seurat,Seurat-method
#'
setMethod("sceInput", "Seurat",  function(input, batch) {
  stopifnot(batch %in% colnames(input[[]]))
  as.SingleCellExperiment(input)
})

#' @rdname sceInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @aliases sceInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("sceInput", "SingleCellExperiment",  function(input, batch) {
  stopifnot(batch %in% colnames(colData(input)))
  input <- input
})

#' @rdname sceInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @aliases sceInput,AnnDataR6,AnnDataR6-method
#'
setMethod("sceInput", "AnnDataR6",  function(input, batch) {
  stopifnot(batch %in% colnames(input$obs))
  input <- AnnData2SCE(input)
})
