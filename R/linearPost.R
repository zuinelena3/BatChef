#' Convert the output of linear models
#'
#' Convert the output of linear models (limma and ComBat) into a SingleCellExperiment, Seurat or Anndata objects.
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param output Linear models output
#' @param method A string specifying the correction method.
#'
#' @import methods
#' @rdname linearPost
#'
setGeneric("linearPost", function(input, output, method)
  standardGeneric("linearPost"), signature = c("input"))

#' @rdname linearPost
#' @aliases linearPost,Seurat,Seurat-method
#' @import methods
#' @importFrom SeuratObject CreateAssay5Object
#'
setMethod("linearPost", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateAssay5Object(data = output)
  return(input)
})

#' @rdname linearPost
#' @aliases linearPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SummarizedExperiment assay<-
#'
setMethod("linearPost", "SingleCellExperiment",  function(input, output, method) {
  assay(input, method) <- output
  return(input)
})

#' @rdname linearPost
#' @aliases linearPost,AnnDataR6,AnnDataR6-method
#'
setMethod("linearPost", "AnnDataR6",  function(input, output, method) {
  input$layers[method] <- t(output)
  return(input)
})
