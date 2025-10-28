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
setGeneric("clustInput", function(input, reduction)
  standardGeneric("clustInput"), signature = c("input"))

#' @rdname clustInput
#' @import methods
#' @importFrom Seurat VariableFeatures Embeddings
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom S4Vectors DataFrame
#' @aliases clustInput,Seurat,Seurat-method
#'
setMethod("clustInput", "Seurat",  function(input, reduction) {
  hvgs <- rownames(input) %in%  VariableFeatures(input)

  red <- Embeddings(object = input, reduction = reduction)
  input <- SingleCellExperiment(assays = list(counts = input@assays$RNA$counts,
                                              logcounts = input@assays$RNA$data),
                                colData = DataFrame(input@meta.data),
                                rowData = data.frame(hvgs = hvgs))
  reducedDim(input, reduction) <- red
  return(input)
})

#' @rdname clustInput
#' @import methods
#' @aliases clustInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("clustInput", "SingleCellExperiment",  function(input, reduction) {
  input <- input
})

#' @rdname clustInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom SingleCellExperiment reducedDims<-
#' @importFrom reticulate r_to_py
#'
#' @aliases clustInput,AnnDataR6,AnnDataR6-method
#'
setMethod("clustInput", "AnnDataR6",  function(input, reduction) {
  red <- input$obsm[[reduction]]
  colnames(red) <- paste0(reduction, 1:ncol(red))
  rownames(red) <- input$obs_names

  input$varm <- NULL
  input$obsm <- NULL

  input <- r_to_py(input, convert = TRUE)
  input <- AnnData2SCE(input, X_name = "counts")
  reducedDim(input, reduction) <- red
  return(input)
})
