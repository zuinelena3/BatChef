#' Title
#'
#' @param input input
#' @param output output
#' @param method method
#'
#' @import methods
#' @rdname linearPost
#'
setGeneric("linearPost", function(input, output, method)
  standardGeneric("linearPost"), signature = c("input"))

#' @rdname linearPost
#' @aliases linearPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject
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
