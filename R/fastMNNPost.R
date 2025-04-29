#' Title
#'
#' @param input input
#' @param output output
#' @param method method
#'
#' @import methods
#' @rdname fastMNNPost
#'
setGeneric("fastMNNPost", function(input, output, method)
  standardGeneric("fastMNNPost"), signature = c("input"))

#' @rdname fastMNNPost
#' @aliases fastMNNPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay rowData
#'
setMethod("fastMNNPost", "Seurat",  function(input, output, method) {
  input[[paste0(method, "_mat")]] <- CreateAssayObject(data = as.matrix(assay(output, "reconstructed")),
                                                       slot = 'data')
  input[[method]] <- CreateDimReducObject(embeddings = reducedDim(output, "corrected"),
                                          loadings = rowData(output)$rotation,
                                          key = "fastmnn_", assay = DefaultAssay(input))
  return(input)
})

#' @rdname fastMNNPost
#' @aliases fastMNNPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SummarizedExperiment assay<- assay
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("fastMNNPost", "SingleCellExperiment",  function(input, output, method) {
  assay(input, method) <- as.matrix(assay(output, "reconstructed"))
  reducedDim(input, method) <- reducedDim(output, "corrected")
  return(input)
})

#' @rdname fastMNNPost
#' @aliases fastMNNPost,AnnDataR6,AnnDataR6-method
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("fastMNNPost", "AnnDataR6",  function(input, output, method) {
  input$layers[method] <- t(as.matrix(assay(output, "reconstructed")))
  input$obsm[[method]] <- as.matrix(reducedDim(output, "corrected"))
  return(input)
})
