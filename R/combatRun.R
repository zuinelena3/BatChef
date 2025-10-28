#' ComBat method
#'
#' ComBat allows users to adjust for batch effects in datasets using an
#' empirical Bayes framework.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param batch A string specifying the batch for each cell.
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @importFrom sva ComBat_seq
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment counts colData
#'
#' @return A corrected matrix
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = FALSE,
#'                      output_format = "SingleCellExperiment")
#' combat <- combatRun(input = sim, batch = "Batch")
#'
combatRun <- function(input, batch, assay_type = "counts", ...) {
  args <- c(list(counts = as.matrix(assay(input, assay_type)),
                 batch = colData(input)[, batch]),
            capture_params(ComBat_seq, list(...)))
  do.call(ComBat_seq, args)
}
