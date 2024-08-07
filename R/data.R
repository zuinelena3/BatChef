#' Simulated data
#'
#' @param nGenes number of genes
#' @param batchCells number of cells for each batch
#' @param group.prob percent of cell types
#'
#' @importFrom splatter newSplatParams splatSimulateGroups
#' @importFrom scuttle logNormCounts
#' @importFrom scater runPCA
#'
simulated_data <- function(nGenes, batchCells, group.prob) {
  params <- newSplatParams(nGenes = nGenes, batchCells = batchCells,
                           group.prob = group.prob, seed = 33)
  sim <- splatSimulateGroups(params, verbose = FALSE)
  sim <- logNormCounts(sim)
  sim <- runPCA(sim,  ncomponent = 10)
  return(sim)
}
