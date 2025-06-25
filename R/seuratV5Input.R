#' Convert into a Seurat V3 object
#'
#' @param input input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param batch batch A string specifying the batch variable.
#' @param pca_name A string specifying the PCA.
#'
#' @import methods
#' @rdname seuratv5Input
#'
setGeneric("seuratv5Input", function(input, batch, pca_name = NULL)
  standardGeneric("seuratv5Input"), signature = c("input"))

#' @rdname seuratv5Input
#' @import methods
#' @importFrom Seurat ScaleData DefaultAssay Reductions
#' @aliases seuratv5Input,Seurat,Seurat-method
#'
setMethod("seuratv5Input", "Seurat", function(input, batch, pca_name) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot("Error: 'pca' not found in this Seurat object" = "pca" %in% Reductions(input))
  if (is(input[[DefaultAssay(input)]], "Assay")) {
    input[[DefaultAssay(input)]] <- as(object = input[[DefaultAssay(input)]], Class = "Assay5")
    input[[DefaultAssay(input)]] <- split(x = input[[DefaultAssay(input)]], f = input[[batch]][, 1])
    input <- ScaleData(input, verbose = FALSE)
  }
  else {
    input[[DefaultAssay(input)]] <- split(x = input[[DefaultAssay(input)]], f = input[[batch]][, 1])
    input <- ScaleData(input, verbose = FALSE)
  }
})

#' @rdname seuratv5Input
#' @import methods
#' @importFrom Seurat as.Seurat ScaleData DefaultAssay CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
#' @aliases seuratv5Input,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("seuratv5Input", "SingleCellExperiment", function(input, batch, pca_name) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))
  so <- as.Seurat(input)
  so@reductions[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(input, pca_name),
                                                 loadings = attr(reducedDim(input, pca_name), "rotation"),
                                                 key = "pca_", assay = DefaultAssay(so))
  so[[DefaultAssay(so)]] <- as(object = so[[DefaultAssay(so)]], Class = "Assay5")
  so[[DefaultAssay(so)]] <- split(x = so[[DefaultAssay(so)]], f = so[[batch]][, 1])
  so <- ScaleData(so, verbose = FALSE)
})

#' @rdname seuratv5Input
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom Seurat as.Seurat ScaleData DefaultAssay
#' @aliases seuratv5Input,AnnDataR6,AnnDataR6-method
#'
setMethod("seuratv5Input", "AnnDataR6",  function(input, batch, pca_name) {
  stopifnot(batch %in% colnames(input$obs))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  so <- as.Seurat(AnnData2SCE(input))
  so@reductions[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(input$obsm[pca_name][[1]]),
                                                 loadings = as.matrix(input$varm$PCs),
                                                 key = "pca_", assay = DefaultAssay(so))

  so[[DefaultAssay(so)]] <- as(object = so[[DefaultAssay(so)]], Class = "Assay5")
  so[[DefaultAssay(so)]] <- split(x = so[[DefaultAssay(so)]], f = so[[batch]][, 1])
  so <- ScaleData(so, verbose = FALSE)
})
