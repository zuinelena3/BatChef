#' Convert into a SingleCellExperiment object for clustering
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param reduction A string specifying the batch for each cell.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @rdname clustInput
#'
setGeneric("clustInput", function(input, reduction) {
  standardGeneric("clustInput")
}, signature = c("input"))

#' @rdname clustInput
#' @import methods
#' @importFrom Seurat VariableFeatures Embeddings GetAssayData
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom S4Vectors DataFrame
#' @aliases clustInput,Seurat,Seurat-method
#'
setMethod("clustInput", "Seurat", function(input, reduction) {
  hvgs <- rownames(input) %in% VariableFeatures(input)

  red <- Embeddings(object = input, reduction = reduction)
  input <- SingleCellExperiment(
    assays = list(
      counts = GetAssayData(input, assay = "RNA", layer = "counts"),
      logcounts = GetAssayData(input, assay = "RNA", layer = "data")
    ),
    colData = DataFrame(input[[]]),
    rowData = DataFrame(hvgs = hvgs, row.names = rownames(input))
  )
  reducedDim(input, reduction) <- red
  return(input)
})

#' @rdname clustInput
#' @import methods
#' @aliases clustInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("clustInput", "SingleCellExperiment", function(input, reduction) {
  reducedDim(input, "PCA") <- as.matrix(reducedDim(input, "PCA"))
  input
})

#' @rdname clustInput
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
#' @aliases clustInput,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("clustInput", "ANY", function(input, reduction) {
  if (is(input, "InMemoryAnnData")) {
    input <- input$as_SingleCellExperiment(x_mapping = "counts",
                                           assays_mapping = c(logcounts = "normalized"),
                                           reducedDims_mapping = list(
                                             PCA = c(sampleFactors = reduction,
                                                     featureLoadings = names(input$varm))))
    reducedDim(input, "PCA") <- as.matrix(reducedDim(input, "PCA"))
    return(input)
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
