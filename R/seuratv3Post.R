#' Title
#'
#' @param input input
#' @param output output
#' @param method method
#'
#' @import methods
#' @rdname seuratv3Post
#'
setGeneric("seuratv3Post", function(input, output, method)
  standardGeneric("seuratv3Post"), signature = c("input"))

#' @rdname seuratv3Post
#' @aliases seuratv3Post,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject
#' @importFrom SeuratObject LayerData
#'
setMethod("seuratv3Post", "Seurat",  function(input, output, method) {
  mat <- LayerData(output)
  mat <- mat[, order(match(colnames(mat), colnames(input)))]
  mat <- mat[order(match(rownames(mat), rownames(input))), ]
  input[[method]] <- CreateAssayObject(data = mat, slot = "data")
  return(input)
})

#' @rdname seuratv3Post
#' @aliases seuratv3Post,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SummarizedExperiment assay<-
#' @importFrom SeuratObject LayerData
#'
setMethod("seuratv3Post", "SingleCellExperiment",  function(input, output, method) {
  mat <- LayerData(output)
  mat <- mat[, order(match(colnames(mat), colnames(input)))]
  mat <- mat[order(match(rownames(mat), rownames(input))), ]

  if (length(rownames(input)) != length(rownames(output))) {
    print("Warning: when the length of rownames of output is not identical to input one, a subset SingleCellExperiment object is returned")
    input <- input[rownames(mat), ]
    assay(input, method) <- mat
    return(input)
  }
  else {
    assay(input, method) <- mat
    return(input)
  }
})

#' @rdname seuratv3Post
#' @aliases seuratv3Post,AnnDataR6,AnnDataR6-method
#' @importFrom anndata AnnData
#' @importFrom SeuratObject LayerData
#'
setMethod("seuratv3Post", "list", function(input, output, method) {
  SCE2AnnData(as.SingleCellExperiment(output))
})
