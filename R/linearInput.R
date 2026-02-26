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
#' @importFrom Seurat VariableFeatures
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @aliases linearInput,Seurat,Seurat-method
#'
setMethod("linearInput", "Seurat", function(input, batch) {
  stopifnot(batch %in% colnames(input[[]]))

  input <- SingleCellExperiment(
    assays = list(
      counts = input@assays$RNA$counts,
      logcounts = input@assays$RNA$data
    ),
    colData = input@meta.data
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
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom reticulate r_to_py
#'
#' @aliases linearInput,AnnDataR6,AnnDataR6-method
#'
setMethod("linearInput", "AnnDataR6", function(input, batch) {
  stopifnot(batch %in% colnames(input$obs))
  input$varm <- NULL
  input$obsm <- NULL
  input <- r_to_py(input, convert = TRUE)
  input <- AnnData2SCE(input, X_name = "counts")
  return(input)
})
