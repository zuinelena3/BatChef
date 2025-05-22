#' Function to simulate data
#'
#' @param nGenes nGenes
#' @param batchCells batchCells
#' @param batch.facLoc batch.facLoc
#' @param batch.facScale batch.facScale
#' @param batch.rmEffect batch.rmEffect
#' @param mean.rate mean.rate
#' @param mean.shape mean.shape
#' @param lib.loc lib.loc
#' @param lib.scale lib.scale
#' @param lib.norm lib.norm
#' @param out.prob out.prob
#' @param out.facLoc out.facLoc
#' @param out.facScale out.facScale
#' @param group.prob group.prob
#' @param de.prob de.prob
#' @param de.downProb de.downProb
#' @param de.facLoc de.facLoc
#' @param de.facScale de.facScale
#' @param bcv.common bcv.common
#' @param bcv.df bcv.df
#' @param dropout.type dropout.type
#' @param dropout.mid dropout.mid
#' @param dropout.shape dropout.shape
#' @param path.from path.from
#' @param path.nSteps path.nSteps
#' @param path.skew path.skew
#' @param path.nonlinearProb path.nonlinearProb
#' @param path.sigmaFac path.sigmaFac
#' @param seed seed
#' @param ncomp ncomp
#' @param n_hvgs n_hvgs
#'
#' @export
#'
#' @importFrom splatter newSplatParams splatSimulateGroups
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom scrapper centerSizeFactors normalizeCounts modelGeneVariances chooseHighlyVariableGenes
#' @importFrom SingleCellExperiment counts
#' @importFrom DelayedArray DelayedArray
#' @importFrom scater runPCA
#' @importFrom SparseArray colSums
#'
simulated_data <- function(nGenes = 10000, batchCells = 100, batch.facLoc = 0.1, batch.facScale = 0.1, batch.rmEffect = FALSE,
                           mean.rate = 0.3, mean.shape = 0.6, lib.loc = 11, lib.scale = 0.2, lib.norm = FALSE,
                           out.prob = 0.05, out.facLoc = 4, out.facScale = 0.5,
                           group.prob = 1, de.prob = 0.1, de.downProb = 0.5, de.facLoc = 0.1, de.facScale = 0.4,
                           bcv.common = 0.1, bcv.df = 60, dropout.type = "none", dropout.mid = 0, dropout.shape = -1,
                           path.from = 0, path.nSteps = 100, path.skew = 0.5, path.nonlinearProb = 0.1, path.sigmaFac = 0.8, seed = 333, ncomp = 10, n_hvgs = 500) {
  params <- newSplatParams(nGenes = nGenes, batchCells = batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale, batch.rmEffect = batch.rmEffect,
                           mean.rate = mean.rate, mean.shape = mean.shape, lib.loc = lib.loc, lib.scale = lib.scale, lib.norm = lib.norm,
                           out.prob = out.prob, out.facLoc = out.facLoc, out.facScale = out.facScale,
                           group.prob = group.prob, de.prob = de.prob, de.downProb = de.downProb, de.facLoc = de.facLoc, de.facScale = de.facScale,
                           bcv.common = bcv.common, bcv.df = bcv.df, dropout.type = dropout.type, dropout.mid = dropout.mid, dropout.shape = dropout.shape,
                           path.from = path.from, path.nSteps = path.nSteps, path.skew = path.skew, path.nonlinearProb = path.nonlinearProb, path.sigmaFac = path.sigmaFac, seed = seed)
  sim <- splatSimulateGroups(params, verbose = FALSE)
  assays(sim) <- assays(sim)["counts"]
  size_factors <- centerSizeFactors(colSums(counts(sim)))
  normalized <- normalizeCounts(DelayedArray(assay(sim)), size_factors)
  assay(sim, "logcounts") <- as(normalized, "dgCMatrix")

  gene_var <- modelGeneVariances(normalized, num.threads = 2)
  top_hvgs <- chooseHighlyVariableGenes(gene_var$statistics$residuals, top = n_hvgs)

  sim <- sim[top_hvgs, ]
  sim <- runPCA(sim, ntop = n_hvgs, ncomponent = ncomp)
  return(sim)
}

