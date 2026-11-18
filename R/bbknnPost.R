#' Convert BBKNN output
#'
#' Convert BBKNN output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output BBKNN output: a matrix of the coordinates of the points chosen
#' to represent the dissimilarities.
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname bbknnPost
#'
setGeneric("bbknnPost", function(input, output, method) {
  standardGeneric("bbknnPost")
}, signature = c("input"))

#' @rdname bbknnPost
#' @aliases bbknnPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("bbknnPost", "Seurat", function(input, output, method) {
  input[[method]] <- CreateDimReducObject(
    embeddings = output,
    key = "bbknn_",
    assay = DefaultAssay(input)
  )
  return(input)
})

#' @rdname bbknnPost
#' @aliases bbknnPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("bbknnPost", "SingleCellExperiment", function(input, output, method) {
  reducedDim(input, method) <- output
  return(input)
})

#' @rdname bbknnPost
#' @aliases bbknnPost,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("bbknnPost", "ANY", function(input, output, method) {
  if (is(input, "InMemoryAnnData")) {
    input$obsm[[method]] <- as.matrix(output)
    return(input)
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
