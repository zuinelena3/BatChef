#' Title
#'
#' @param input input
#' @param output output
#'
#' @import methods
#' @rdname seuratv5Post
#'
setGeneric("seuratv5Post", function(input, output)
  standardGeneric("seuratv5Post"), signature = c("input"))

#' @rdname seuratv5Post
#' @aliases seuratv5Post,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject JoinLayers
#'
setMethod("seuratv5Post", "Seurat",  function(input, output) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  return(output)
})

#' @rdname seuratv5Post
#' @aliases seuratv5Post,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment
#' @importFrom SeuratObject JoinLayers
#'
setMethod("seuratv5Post", "SingleCellExperiment",  function(input, output) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  output <- as.SingleCellExperiment(output)
})

#' @rdname seuratv5Post
#' @aliases seuratv5Post,AnnDataR6,AnnDataR6-method
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment
#' @importFrom SeuratObject JoinLayers
#' @importFrom zellkonverter SCE2AnnData
#'
setMethod("seuratv5Post", "AnnDataR6",  function(input, output) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  sce <- as.SingleCellExperiment(output)
  adata <- SCE2AnnData(sce)
  return(adata)
})
