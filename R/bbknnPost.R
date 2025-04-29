#' Title
#'
#' @param input input
#' @param output output
#' @param method method
#'
#' @import methods
#' @rdname bbknnPost
#'
setGeneric("bbknnPost", function(input, output, method)
  standardGeneric("bbknnPost"), signature = c("input"))

#' @rdname bbknnPost
#' @aliases bbknnPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("bbknnPost", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateDimReducObject(embeddings = output,
                                          key = "bbknn_", assay = DefaultAssay(input))
  return(input)
})

#' @rdname bbknnPost
#' @aliases bbknnPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("bbknnPost", "SingleCellExperiment",  function(input, output, method) {
  reducedDim(input, method) <- output
  return(input)
})

#' @rdname bbknnPost
#' @aliases bbknnPost,AnnDataR6,AnnDataR6-method
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("bbknnPost", "AnnDataR6",  function(input, output, method) {
  input$obsm[[method]] <- as.matrix(output)
  return(input)
})
