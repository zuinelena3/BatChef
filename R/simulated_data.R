#' Simulated data (with splatter package)
#'
#' @param nGenes number of genes
#' @param batchCells number of cells for each batch
#' @param group.prob percent of cell types
#' @param ncomp number of principal components
#'
#' @importFrom splatter newSplatParams splatSimulateGroups
#' @importFrom SummarizedExperiment assays<- assays
#' @importFrom scuttle logNormCounts
#' @importFrom scater runPCA
#'
#' @export
#'
simulated_data <- function(nGenes, batchCells, group.prob, ncomp) {
  params <- newSplatParams(nGenes = nGenes, batchCells = batchCells,
                           group.prob = group.prob, seed = 33)
  sim <- splatSimulateGroups(params, verbose = FALSE)
  assays(sim) <- assays(sim)["counts"]
  sim <- logNormCounts(sim)
  sim <- runPCA(sim,  ncomponent = ncomp)
  return(sim)
}
