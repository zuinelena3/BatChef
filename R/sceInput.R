#' Convert into a SingleCellExperiment object
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch for each cell.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @rdname sceInput
#'
setGeneric("sceInput", function(input, batch) {
  standardGeneric("sceInput")
}, signature = c("input"))

#' @rdname sceInput
#' @import methods
#' @importFrom Seurat VariableFeatures Embeddings GetAssayData
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom S4Vectors DataFrame
#' @aliases sceInput,Seurat,Seurat-method
#'
setMethod("sceInput", "Seurat", function(input, batch) {
  stopifnot(batch %in% colnames(input[[]]))

  hvgs <- rownames(input) %in% VariableFeatures(input)

  pca <- Embeddings(input)
  input <- SingleCellExperiment(
    assays = list(
      counts = GetAssayData(input, assay = "RNA", layer = "counts"),
      data = GetAssayData(input, assay = "RNA", layer = "data")
    ),
    colData = DataFrame(input[[]]),
    rowData = DataFrame(hvgs = hvgs, row.names = rownames(input))
  )
  reducedDim(input, "PCA") <- pca
  return(input)
})

#' @rdname sceInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @aliases sceInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("sceInput", "SingleCellExperiment", function(input, batch) {
  stopifnot(batch %in% colnames(colData(input)))
  input <- input
})

#' @rdname sceInput
#' @import methods
#' @importFrom SingleCellExperiment reducedDims<-
#'
#' @aliases sceInput,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("sceInput", "ANY", function(input, batch) {
  stopifnot(batch %in% colnames(input$obs))
  if (is(input, "InMemoryAnnData")) {
  input <- input$as_SingleCellExperiment(reducedDims_mapping = list(PCA = c(sampleFactors = "X_pca", featureLoadings = "PCs")))
  reducedDim(input, "PCA") <- as.matrix(reducedDim(input, "PCA"))
  return(input)
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
