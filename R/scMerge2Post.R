#' Convert the scMerge2 output
#'
#' Convert the scMerge2 output into a SingleCellExperiment, Seurat or Anndata objects
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param output scMerge2 output.
#' @param method A string specifying the correction method
#'
#' @import methods
#' @rdname scMerge2Post
#'
setGeneric("scMerge2Post", function(input, output, method)
  standardGeneric("scMerge2Post"), signature = c("input"))

#' @rdname scMerge2Post
#' @aliases scMerge2Post,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject
#'
setMethod("scMerge2Post", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateAssayObject(data = as.matrix(output[["newY"]]),
                                       slot = 'data')
  return(input)
})

#' @rdname scMerge2Post
#' @aliases scMerge2Post,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#'
setMethod("scMerge2Post", "SingleCellExperiment",  function(input, output, method) {
  assay(input, method) <- output[["newY"]]
  return(input)
})

#' @rdname scMerge2Post
#' @aliases scMerge2Post,AnnDataR6,AnnDataR6-method
#'
setMethod("scMerge2Post", "AnnDataR6",  function(input, output, method) {
  input$layers[method] <- t(as.matrix(output[["newY"]]))
  return(input)
})
