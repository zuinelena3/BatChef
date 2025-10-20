#' Convert to a SeuratV3 compatible object
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param features Vector of features to use.
#' @param pca_name A string specifying the PCA.
#'
#' @import methods
#' @return A list of \link[Seurat]{Seurat} objects.
#' @rdname seuratv3Input
#'
setGeneric("seuratv3Input", function(input, batch, features, pca_name = NULL)
  standardGeneric("seuratv3Input"), signature = c("input"))

#' @rdname seuratv3Input
#' @import methods
#' @importFrom Seurat DefaultAssay SplitObject VariableFeatures<-
#' @importFrom scCustomize Convert_Assay
#' @aliases seuratv3Input,Seurat,Seurat-method
#'
setMethod("seuratv3Input", "Seurat", function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot("Error: 'pca' not found in this Seurat object" = "pca" %in% Reductions(input))

  if (is(input[[DefaultAssay(input)]], "Assay5")) {
    input <-  Convert_Assay(seurat_object = input, assay = DefaultAssay(input),
                            convert_to = "V3")
  }

  else {
    input <- input
  }

  VariableFeatures(input) <- features
  ll <- SplitObject(input, split.by = batch)
})

#' @rdname seuratv3Input
#' @import methods
#' @importFrom Seurat as.Seurat SplitObject VariableFeatures<- DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#' @aliases seuratv3Input,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("seuratv3Input", "SingleCellExperiment", function(input, batch,
                                                            features, pca_name) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  pca <- reducedDim(input, pca_name)
  loadings <- attr(reducedDim(input, pca_name), "rotation")
  reducedDim(input, pca_name) <- NULL
  so <- as.Seurat(input)
  VariableFeatures(so) <- features
  so[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                      loadings = loadings,
                                      key = "pca_", assay = DefaultAssay(so))
  ll <- SplitObject(so, split.by = batch)
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
#' @aliases seuratv3Input,AnnDataR6,AnnDataR6-method
#'
setMethod("seuratv3Input", "AnnDataR6",  function(input, batch, features, pca_name) {
  stopifnot("Error: 'batch' is not inside the obs" =  batch %in% colnames(input$obs))

  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  var <- unique(input$obs[[batch]])
  ll <- lapply(var, function(v) {
    sub <- input[input$obs[[batch]] == v, ]$copy()
  })

  so_ll <- lapply(ll, function(x) {
    pca <- x$obsm[[pca_name]]
    colnames(pca) <- paste0("PC_", 1:ncol(pca))
    rownames(pca) <- x$obs_names
    loadings <- x$varm[["PCs"]]
    colnames(loadings) <- paste0("PC_", 1:ncol(pca))
    rownames(loadings) <- x$var_names
    x$obsm <- NULL
    x$varm <- NULL

    x <- r_to_py(x, convert = TRUE)
    sce <- AnnData2SCE(x, X_name = "counts")
    so <- as.Seurat(sce)
    so@reductions[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                                   loadings = loadings,
                                                   key = "pca_",
                                                   assay = DefaultAssay(so))
    VariableFeatures(so) <- features
    return(so)
  })
})
