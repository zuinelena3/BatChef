#' Convert to a SeuratV5 compatible object
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param features Vector of features to use.
#' @param pca_name A string specifying the PCA.
#'
#' @import methods
#' @return A \linkS4class{Seurat} object.
#' @rdname seuratv5Input
#'
setGeneric("seuratv5Input", function(input, batch, features = NULL, pca_name = NULL)
  standardGeneric("seuratv5Input"), signature = c("input"))

#' @rdname seuratv5Input
#' @import methods
#' @importFrom Seurat ScaleData DefaultAssay Reductions
#' @importFrom scCustomize Convert_Assay
#' @aliases seuratv5Input,Seurat,Seurat-method
#'
setMethod("seuratv5Input", "Seurat", function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot("Error: 'pca' not found in this Seurat object" = "pca" %in% Reductions(input))

  if (is(input[[DefaultAssay(input)]], "Assay")) {
    input <-  Convert_Assay(seurat_object = input, assay = DefaultAssay(input),
                            convert_to = "V5")
  }
  else {
    input
  }

  input[[DefaultAssay(input)]] <- split(x = input[[DefaultAssay(input)]],
                                        f = input[[batch]][, 1])
  input <- ScaleData(input, verbose = FALSE)
})

#' @rdname seuratv5Input
#' @import methods
#' @importFrom Seurat as.Seurat ScaleData DefaultAssay CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scCustomize Convert_Assay
#' @aliases seuratv5Input,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("seuratv5Input", "SingleCellExperiment", function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  pca <- reducedDim(input, pca_name)
  loadings <- attr(reducedDim(input, pca_name), "rotation")
  reducedDim(input, pca_name) <- NULL

  so <- as.Seurat(input)

  so@reductions[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                                 loadings = loadings,
                                                 key = "pca_",
                                                 assay = DefaultAssay(so))

  so <-  Convert_Assay(seurat_object = so, assay = DefaultAssay(so),
                            convert_to = "V5")
  VariableFeatures(so) <- features
  so[[DefaultAssay(so)]] <- split(x = so[[DefaultAssay(so)]],
                                        f = so[[batch]][, 1])
  so <- ScaleData(so, verbose = FALSE)
})

#' @rdname seuratv5Input
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom Seurat as.Seurat ScaleData DefaultAssay
#' @aliases seuratv5Input,AnnDataR6,AnnDataR6-method
#'
setMethod("seuratv5Input", "AnnDataR6",  function(input, batch, features, pca_name) {
  stopifnot(batch %in% colnames(input$obs))
  stopifnot("Error: specify 'pca_name'" = !is.null(pca_name))

  pca <- input$obsm[[pca_name]]
  colnames(pca) <- paste0("PC_", 1:ncol(pca))
  rownames(pca) <- input$obs_names
  loadings <- input$varm[["PCs"]]
  colnames(loadings) <- paste0("PC_", 1:ncol(pca))
  rownames(loadings) <- input$var_names
  input$obsm <- NULL
  input$varm <- NULL

  input <- r_to_py(input, convert = TRUE)
  sce <- AnnData2SCE(input, X_name = "counts")
  so <- as.Seurat(sce)

  so@reductions[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                                 loadings = loadings,
                                                 key = "pca_",
                                                 assay = DefaultAssay(so))

  so <-  Convert_Assay(seurat_object = so, assay = DefaultAssay(so),
                       convert_to = "V5")
  VariableFeatures(so) <- features
  so[[DefaultAssay(so)]] <- split(x = so[[DefaultAssay(so)]],
                                  f = so[[batch]][, 1])
  so <- ScaleData(so, verbose = FALSE)
})
