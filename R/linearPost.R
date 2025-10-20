#' Convert the output of linear models
#'
#' Convert a corrected matrix output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output A corrected matrix.
#' @param method A string specifying the correction method.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
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
  input[[method]] <- CreateAssayObject(counts = output)
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
