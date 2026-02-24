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
#' @importFrom Seurat VariableFeatures Embeddings
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
      counts = input@assays$RNA$counts,
      logcounts = input@assays$RNA$data
    ),
    colData = DataFrame(input@meta.data),
    rowData = data.frame(hvgs = hvgs)
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
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom SingleCellExperiment reducedDims<-
#' @importFrom reticulate r_to_py
#'
#' @aliases sceInput,AnnDataR6,AnnDataR6-method
#'
setMethod("sceInput", "AnnDataR6", function(input, batch) {
  stopifnot(batch %in% colnames(input$obs))

  pca <- input$obsm[["X_pca"]]
  colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
  rownames(pca) <- input$obs_names
  loadings <- input$varm[["PCs"]]
  colnames(loadings) <- paste0("PC_", seq_len(ncol(pca)))
  rownames(loadings) <- input$var_names
  input$varm <- NULL
  input$obsm <- NULL

  input <- r_to_py(input, convert = TRUE)
  input <- AnnData2SCE(input, X_name = "counts")
  reducedDim(input, "PCA") <- pca
  attr(reducedDim(input, "PCA"), "rotation") <- loadings
  return(input)
})
