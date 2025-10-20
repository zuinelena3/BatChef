#' Convert Harmony output
#'
#' Convert Harmony output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output Harmony output: a \link[SingleCellExperiment]{SingleCellExperiment}
#' containing the corrected dimensionality reduction space.
#' @param method A string specifying the correction method.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname harmonyPost
#'
setGeneric("harmonyPost", function(input, output, method)
  standardGeneric("harmonyPost"), signature = c("input"))

#' @rdname harmonyPost
#' @aliases harmonyPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("harmonyPost", "Seurat",  function(input, output, method) {
  input[[method]] <- CreateDimReducObject(embeddings = reducedDim(output, "HARMONY"),
                                          key = "harmony_",
                                          assay = DefaultAssay(input))
  return(input)
})

#' @rdname harmonyPost
#' @aliases harmonyPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("harmonyPost", "SingleCellExperiment",  function(input, output, method) {
  reducedDim(input, method) <- reducedDim(output, "HARMONY")
  return(input)
})

#' @rdname harmonyPost
#' @aliases harmonyPost,AnnDataR6,AnnDataR6-method
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("harmonyPost", "AnnDataR6",  function(input, output, method) {
  input$obsm[[method]] <- as.matrix(reducedDim(output, "HARMONY"))
  return(input)
})
