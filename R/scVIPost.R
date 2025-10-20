#' Convert the scVI output
#'
#' Convert the scVI output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output scVI output: a list that contains the corrected
#' gene expression matrix and the corrected low-dimensional space.
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname scVIPost
#'
setGeneric("scVIPost", function(input, output, method)
  standardGeneric("scVIPost"), signature = c("input"))

#' @rdname scVIPost
#' @aliases scVIPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim

setMethod("scVIPost", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateDimReducObject(embeddings = as.matrix(output[[2]]),
                                          key = "scvi_", assay = DefaultAssay(input))
  return(input)
})

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

#' @rdname scVIPost
#' @aliases scVIPost,AnnDataR6,AnnDataR6-method
#' @import methods
#'
setMethod("scVIPost", "AnnDataR6",  function(input, output, method) {
  input$layers["scvi"] <- output[[1]]
  input$obsm[["scvi"]] <- output[[2]]

  return(input)
})
