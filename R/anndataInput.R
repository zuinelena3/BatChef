#' Title
#'
#' @param input input
#'
#' @export
#' @import methods
#' @rdname anndataInput
#'
setGeneric("anndataInput", function(input)
  standardGeneric("anndataInput"), signature = c("input"))

#' @rdname anndataInput
#' @import methods
#' @importFrom Seurat as.SingleCellExperiment
#' @aliases anndataInput,Seurat,Seurat-method
#'
setMethod("anndataInput", "Seurat",  function(input) {
  sce <- as.SingleCellExperiment(input)
  adata <- SCE2AnnData(sce)
})

#' @rdname anndataInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @aliases anndataInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("anndataInput", "SingleCellExperiment",  function(input) {
  adata <- SCE2AnnData(input)
})

#' @rdname anndataInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @aliases anndataInput,anndata.AnnData,anndata.AnnData-method
#'
setMethod("anndataInput", "AnnDataR6",  function(input) {
  input <- input
})
