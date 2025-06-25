#' Convert into a Seurat V3 object
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param batch A string specifying the batch variable.
#' @param features Vector of features to use.
#' @param pca_name A string specifying the PCA.
#'
#' @import methods
#' @rdname seuratv3Input
#'
setGeneric("seuratv3Input", function(input, batch, features, pca_name = NULL)
  standardGeneric("seuratv3Input"), signature = c("input"))

#' @rdname seuratv3Input
#' @import methods
#' @importFrom Seurat SplitObject VariableFeatures<-
#' @aliases seuratv3Input,Seurat,Seurat-method
#'
setMethod("seuratv3Input", "Seurat", function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot("Error: 'pca' not found in this Seurat object" = "pca" %in% Reductions(input))
  VariableFeatures(input) <- features
  SplitObject(input, split.by = batch)
})

#' @rdname seuratv3Input
#' @import methods
#' @importFrom Seurat as.Seurat SplitObject VariableFeatures<- DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#' @aliases seuratv3Input,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("seuratv3Input", "SingleCellExperiment", function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))
  so <- as.Seurat(input)
  VariableFeatures(so) <- features
  so[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(input, pca_name),
                                      loadings = attr(reducedDim(input, pca_name), "rotation"),
                                      key = "pca_", assay = DefaultAssay(so))
  SplitObject(so, split.by = batch)
})

#' @rdname seuratv3Input
#' @import methods
#' @aliases seuratv3Input,AnnDataR6,AnnDataR6-method
#'
setMethod("seuratv3Input", "AnnDataR6",  function(input, batch, features, pca_name) {
  stopifnot("The input is an AnnData object. Please create a list of AnnData!" = !is(input, "AnnDataR6"))
})

#' @rdname seuratv3Input
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom Seurat CreateDimReducObject VariableFeatures<- ScaleData DefaultAssay
#' @aliases seuratv3Input,list,list-method
#'
setMethod("seuratv3Input", "list",  function(input, batch, features, pca_name) {
  stopifnot("Error: elements inside the list must be an AnnData" = sapply(input, function(x) all(is(x, "AnnDataR6"))))

  stopifnot("Error: 'batch' is not inside the obs" = sapply(input, function(x) batch %in% colnames(x$obs)))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))
  so_ll <- lapply(input, function(x) {
    so <- as.Seurat(AnnData2SCE(x))
    so@reductions[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(x$obsm[pca_name][[1]]),
                                                   loadings = as.matrix(x$varm$PCs),
                                                   key = "pca_", assay = DefaultAssay(so))
    so <- ScaleData(so, verbose = FALSE)
    VariableFeatures(so) <- features
    so
  })
})
