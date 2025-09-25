#' fastMNN method
#'
#' Correct for batch effects in single-cell expression data using a fast
#' version of the mutual nearest neighbors (MNN) method.
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
#' @param batch A string specifying the batch variable.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment colData
#'
#' @returns A \linkS4class{SingleCellExperiment} object is returned where
#' each row is a gene and each column is a cell. This contains a corrected
#' low-dimensional coordinates and a reconstructed matrix in the assays slot.
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 50),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' fastmnn <- fastMNNRun(input = sim, batch = "Batch")
#'
fastMNNRun <- function(input, batch, ...) {
  args <- c(list(counts = input, batch = colData(input)[, batch]),
            capture_params(fastMNN, list(...)))
  do.call(fastMNN, args)
}
