#' scIntegration
#'
#' Generic function for defining batch effect correction methods.
#'
#' @param obj a \linkS4class{SingleCellExperiment} object or list containing single-cell gene expression matrices
#' @param batch a string specifying the batch
#' @param assay a string specifying the assay to use for correction
#' @param hvgs a vector specifying which features to use for correction
#' @param dims number of dimensions to use for dimension reduction
#' @param reduction a string specifying the dimension reduction to use for correction
#' @param anchor a string specifying the anchor integration type for \linkS4class{Seurat} (V3 and V5) method
#' @param k_anchor number of anchors
#' @param genelist list containing the genes in each batch
#' @param cell_type string specifying the cell-type labels
#' @param METHOD a \code{\link{batchMethod}} object specifying the batch correction method
#' @param alt_out alternative output: a \code{\link{saveClass}} class
#'
#' @export
#' @import  methods
#'
setGeneric("scIntegration", function(obj, batch = NULL, assay = NULL,
                                     hvgs = NULL, dims = NULL, reduction = NULL,
                                     anchor = NULL, k_anchor = NULL, genelist = NULL,
                                     cell_type = NULL, METHOD, alt_out = FALSE)
             standardGeneric("scIntegration"), signature = c("METHOD"))


