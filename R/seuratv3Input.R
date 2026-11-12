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
setGeneric("seuratv3Input", function(input, batch, features, pca_name = NULL) {
  standardGeneric("seuratv3Input")
}, signature = c("input"))

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
    input <- Convert_Assay(
      seurat_object = input, assay = DefaultAssay(input),
      convert_to = "V3"
    )
  } else {
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

  # NOTE: until https://github.com/satijalab/seurat/issues/9850 is fixed.
  suppressWarningsByMsg("deprecated", {
    so <- as.Seurat(input)
  })

  so[["pca"]] <- CreateDimReducObject(
    embeddings = pca,
    loadings = loadings,
    key = "pca_", assay = DefaultAssay(so)
  )
  ll <- SplitObject(so, split.by = batch)
  ll <- lapply(ll, function(x){
    VariableFeatures(x) <- features
    x
  })
})

#' @rdname seuratv3Input
#' @import methods
#' @importFrom Seurat CreateDimReducObject VariableFeatures<- ScaleData DefaultAssay
#' @aliases seuratv3Input,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("seuratv3Input", "ANY", function(input, batch, features, pca_name) {
  stopifnot("Error: 'batch' is not inside the obs" = batch %in% colnames(input$obs))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  if (is(input, "InMemoryAnnData")) {
    input$X <- as(input$X, "dgCMatrix")
    var <- unique(input$obs[[batch]])

    if (is.null(input$layers)) {
      so <- input$as_Seurat(assay_name = "RNA", x_mapping = "counts")
    } else {
      so <- input$as_Seurat(assay_name = "RNA", x_mapping = "counts",
                            layers_mapping = c(data = "normalized"),
                            reduction_mapping = list(pca = c(key = "PC_",
                                                             embeddings = pca_name,
                                                             loadings = names(input$varm))))
    }

    ll <- lapply(var, function(v) {
      x <- so[, so@meta.data[[batch]] == v]
      VariableFeatures(x) <- features
      x
    })
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
