#' Title
#'
#' @param input input
#' @param output output
#' @param method method
#'
#' @import methods
#' @rdname scVIPost
#'
setGeneric("scVIPost", function(input, output, method)
  standardGeneric("scVIPost"), signature = c("input"))

# #' @rdname scVIPost
# #' @aliases scVIPost,Seurat,Seurat-method
# #' @import methods
# #' @importFrom Seurat CreateDimReducObject DefaultAssay
# #' @importFrom SingleCellExperiment reducedDim
#
# setMethod("scVIPost", "Seurat",  function(input, output, method) {
#   input[[method]] <- CreateDimReducObject(embeddings = as.matrix(output[[2]]),
#                                           key = "scvi_", assay = DefaultAssay(input))
#   return(input)
# })

#' @rdname scVIPost
#' @aliases scVIPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("scVIPost", "SingleCellExperiment",  function(input, output, method) {
  assay(input, method) <- t(output[[1]])
  reducedDim(input, method) <- output[[2]]
  return(input)
})

# #' @rdname scVIPost
# #' @aliases scVIPost,AnnDataR6,AnnDataR6-method
# #' @importFrom SingleCellExperiment reducedDim
#
# # setMethod("scVIPost", "AnnDataR6",  function(input, output, method) {
# #   input$
# #   return(input)
# # })
