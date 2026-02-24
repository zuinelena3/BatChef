#' Convert into a BBKNN compatible object
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param reduction A string specifying the name of PCA.
#'
#' @import methods
#' @return A BBKNN compatible object
#' @rdname bbknnInput
#'
setGeneric("bbknnInput", function(input, batch, reduction) {
  standardGeneric("bbknnInput")
}, signature = c("input"))

#' @rdname bbknnInput
#' @import methods
#' @importFrom Seurat as.SingleCellExperiment Embeddings Reductions
#' @aliases bbknnInput,Seurat,Seurat-method
#'
setMethod("bbknnInput", "Seurat", function(input, batch, reduction) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot(reduction %in% Reductions(input))

  red <- Embeddings(input, reduction = reduction)
  batch <- input[[batch]][, 1]
  return(list(red = red, batch = batch))
})

#' @rdname bbknnInput
#' @import methods
#' @importFrom SingleCellExperiment colData reducedDim reducedDimNames
#' @aliases bbknnInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("bbknnInput", "SingleCellExperiment", function(input, batch, reduction) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot(reduction %in% reducedDimNames(input))

  red <- reducedDim(input, reduction)
  batch <- colData(input)[, batch]
  return(list(red = red, batch = batch))
})

#' @rdname bbknnInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @aliases bbknnInput,AnnDataR6,AnnDataR6-method
#'
setMethod("bbknnInput", "AnnDataR6", function(input, batch, reduction) {
  stopifnot(batch %in% colnames(input$obs))
  stopifnot(reduction %in% names(input$obsm))

  red <- as.matrix(input$obsm[reduction][[1]])
  batch <- input$obs[batch][, 1]
  return(list(red = red, batch = batch))
})
