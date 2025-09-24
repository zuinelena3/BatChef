#' ComBat method
#'
#' ComBat allows users to adjust for batch effects in datasets using an
#' empirical Bayes framework.
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
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
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' combat <- combatRun(input = sim, batch = "Batch")
#'
combatRun <- function(input, batch, assay_type = "counts", ...) {
  args <- c(list(counts = as.matrix(assay(input, assay_type)),
                 batch = colData(input)[, batch]),
            capture_params(ComBat_seq, list(...)))
  do.call(ComBat_seq, args)
}
