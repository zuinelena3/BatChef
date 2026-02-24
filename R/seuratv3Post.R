#' Convert the SeuratV3 output
#'
#' Convert the SeuratV3 output into a  \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A  \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output Seurat V3 output: a Seurat object
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return A  \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname seuratv3Post
#'
setGeneric("seuratv3Post", function(input, output, method) {
  standardGeneric("seuratv3Post")
}, signature = c("input"))

#' @rdname seuratv3Post
#' @aliases seuratv3Post,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject
#' @importFrom SeuratObject LayerData
#'
setMethod("seuratv3Post", "Seurat", function(input, output, method) {
  mat <- LayerData(output)
  mat <- mat[
    order(match(rownames(mat), rownames(input))),
    order(match(colnames(mat), colnames(input)))
  ]
  input[[method]] <- CreateAssayObject(counts = mat)
  return(input)
})

#' @rdname seuratv3Post
#' @aliases seuratv3Post,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<-
#' @importFrom SeuratObject LayerData
#'
setMethod("seuratv3Post", "SingleCellExperiment", function(input, output, method) {
  mat <- LayerData(output)
  alt_se <- SingleCellExperiment(assays = list(counts = mat))

  altExp(input, method) <- alt_se
  return(input)
})

#' @rdname seuratv3Post
#' @aliases seuratv3Post,AnnDataR6,AnnDataR6-method
#' @importFrom SeuratObject LayerData
#' @importFrom anndata AnnData
#'
setMethod("seuratv3Post", "AnnDataR6", function(input, output, method) {
  mat <- t(as.matrix(LayerData(output)))
  list <- list(original = input, corrected = AnnData(X = mat, obs = input$obs))
  return(list)
})
