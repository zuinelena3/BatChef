#' scMerge2 method
#'
#' scMerge2 is based on pseudo-replicates to remove unwanted batch effects.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param batch A string specifying the batch variable.
#' @param assay_type A string specifying the assay to use for correction.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom scMerge scMerge2
#' @importFrom SummarizedExperiment assay colData
#'
#' @return A list that contains the corrected gene expression matrix
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(250, 200),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' scmerge2 <- scMerge2Run(input = sim, batch = "Batch", assay_type = "logcounts")
#'
scMerge2Run <- function(input, batch, assay_type = "logcounts", ...) {
  args <- c(
    list(exprsMat = assay(input, assay_type), batch = colData(input)[, batch]),
    capture_params(scMerge2, list(...))
  )
  do.call(scMerge2, args)
}
