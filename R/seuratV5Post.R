#' Convert the Seurat V5 output
#'
#' Convert the Seurat V5 output into a SingleCellExperiment, Seurat or Anndata objects
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param output Seurat V5 output.
#' @param method A string specifying the correction method.
#' @param name A string specifying the corrected reduce space name.
#'
#' @import methods
#' @rdname seuratv5Post
#'
setGeneric("seuratv5Post", function(input, output, method, name)
  standardGeneric("seuratv5Post"), signature = c("input"))

#' @rdname seuratv5Post
#' @aliases seuratv5Post,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject JoinLayers Embeddings
#'
setMethod("seuratv5Post", "Seurat",  function(input, output, method, name) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  input[[method]] <- CreateDimReducObject(embeddings = Embeddings(output, name),
                                          key = "seuratv5_", assay = DefaultAssay(input))

  return(input)
})

#' @rdname seuratv5Post
#' @aliases seuratv5Post,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment
#' @importFrom SeuratObject JoinLayers Embeddings
#'
setMethod("seuratv5Post", "SingleCellExperiment",  function(input, output, method, name) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  reducedDim(input, method) <- Embeddings(output, name)
  return(input)
})

#' @rdname seuratv5Post
#' @aliases seuratv5Post,AnnDataR6,AnnDataR6-method
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment
#' @importFrom SeuratObject JoinLayers
#' @importFrom zellkonverter SCE2AnnData
#'
setMethod("seuratv5Post", "AnnDataR6",  function(input, output, method, name) {
  output[[DefaultAssay(output)]] <- JoinLayers(output[[DefaultAssay(output)]])
  input$obsm[[method]] <- as.matrix(Embeddings(output, name))
  return(input)
})
