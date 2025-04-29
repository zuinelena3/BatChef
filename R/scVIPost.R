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

#' @rdname scVIPost
#' @aliases scVIPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("scVIPost", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateDimReducObject(embeddings = output$obsm["X_scVI"],
                                          key = "scvi_", assay = DefaultAssay(input))
  return(input)
})

#' @rdname scVIPost
#' @aliases scVIPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("scVIPost", "SingleCellExperiment",  function(input, output, method) {
  reducedDim(input, method) <- output$obsm["X_scVI"]
  return(input)
})

#' @rdname scVIPost
#' @aliases scVIPost,AnnDataR6,AnnDataR6-method
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("scVIPost", "AnnDataR6",  function(input, output, method) {
  return(input)
})
