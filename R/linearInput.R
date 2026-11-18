#' Convert into a SingleCellExperiment object for linear model based methods
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch for each cell.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @rdname linearInput
#'
setGeneric("linearInput", function(input, batch) {
  standardGeneric("linearInput")
}, signature = c("input"))

#' @rdname linearInput
#' @import methods
#' @importFrom Seurat VariableFeatures GetAssayData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @aliases linearInput,Seurat,Seurat-method
#'
setMethod("linearInput", "Seurat", function(input, batch) {
  stopifnot(batch %in% colnames(input[[]]))

  hvgs <- rownames(input) %in% VariableFeatures(input)
  input <- SingleCellExperiment(
    assays = list(
      counts = GetAssayData(input, assay = "RNA", layer = "counts"),
      data = GetAssayData(input, assay = "RNA", layer = "data")
    ),
    colData = DataFrame(input[[]]),
    rowData = DataFrame(hvgs = hvgs, row.names = rownames(input))
  )
})

#' @rdname linearInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @aliases linearInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("linearInput", "SingleCellExperiment", function(input, batch) {
  stopifnot(batch %in% colnames(colData(input)))
  input <- input
})

#' @rdname linearInput
#' @import methods
#'
#' @aliases linearInput,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("linearInput", "ANY", function(input, batch) {
  stopifnot(batch %in% colnames(input$obs))

  if (is(input, "InMemoryAnnData")) {
    input <- input$as_SingleCellExperiment()
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
