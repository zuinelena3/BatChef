#' Function to simulate data
#'
#' The function allows to simulated single-cell RNA-seq data using Splatter
#' package, normalize data, select highly variable genes and
#' compute Principal Component Analysis.
#'
#' @param nGenes Number of genes.
#' @param batchCells Number of cells per batch.
#' @param batch.facLoc Batch factor location parameter.
#' @param batch.facScale Batch factor scale parameter.
#' @param batch.rmEffect Remove batch effect.
#' @param mean.rate Mean rate.
#' @param mean.shape Mean shape.
#' @param lib.loc Library size location parameter.
#' @param lib.scale Library size scale parameter.
#' @param lib.norm Library size distribution.
#' @param out.prob Expression outlier probability
#' @param out.facLoc Expression outlier factor location
#' @param out.facScale Expression outlier factor scale
#' @param group.prob Group probabilities
#' @param de.prob Differential expression probability
#' @param de.downProb Down-regulation probability
#' @param de.facLoc DE factor location
#' @param de.facScale DE factor scale
#' @param bcv.common Common Biological Coefficient of Variation
#' @param bcv.df BCV Degrees of Freedom
#' @param dropout.type Dropout type
#' @param dropout.mid Dropout mid point
#' @param dropout.shape Dropout shape
#' @param path.from Path origin
#' @param path.nSteps Number of steps
#' @param path.skew Path skew
#' @param path.nonlinearProb Non-linear probability
#' @param path.sigmaFac Path skew
#' @param seed Random seed
#' @param ncomp Number of principal component
#' @param n_hvgs Number of higly variable genes.
#'
#' @export
#' @importFrom splatter newSplatParams splatSimulateGroups
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom scrapper centerSizeFactors normalizeCounts modelGeneVariances chooseHighlyVariableGenes
#' @importFrom SingleCellExperiment counts
#' @importFrom DelayedArray DelayedArray
#' @importFrom scater runPCA
#' @importFrom SparseArray colSums
#'
#' @return A \linkS4class{SingleCellExperiment} object.
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#'
simulated_data <- function(nGenes = 10000, batchCells = 100, batch.facLoc = 0.1,
                           batch.facScale = 0.1, batch.rmEffect = FALSE,
                           mean.rate = 0.3, mean.shape = 0.6, lib.loc = 11,
                           lib.scale = 0.2, lib.norm = FALSE,
                           out.prob = 0.05, out.facLoc = 4, out.facScale = 0.5,
                           group.prob = 1, de.prob = 0.1, de.downProb = 0.5,
                           de.facLoc = 0.1, de.facScale = 0.4,
                           bcv.common = 0.1, bcv.df = 60, dropout.type = "none",
                           dropout.mid = 0, dropout.shape = -1,
                           path.from = 0, path.nSteps = 100, path.skew = 0.5,
                           path.nonlinearProb = 0.1, path.sigmaFac = 0.8,
                           seed = 333, ncomp = 10, n_hvgs = 500) {
  params <- newSplatParams(nGenes = nGenes, batchCells = batchCells,
                           batch.facLoc = batch.facLoc,
                           batch.facScale = batch.facScale,
                           batch.rmEffect = batch.rmEffect,
                           mean.rate = mean.rate, mean.shape = mean.shape,
                           lib.loc = lib.loc, lib.scale = lib.scale,
                           lib.norm = lib.norm,
                           out.prob = out.prob, out.facLoc = out.facLoc,
                           out.facScale = out.facScale,
                           group.prob = group.prob, de.prob = de.prob,
                           de.downProb = de.downProb, de.facLoc = de.facLoc,
                           de.facScale = de.facScale,
                           bcv.common = bcv.common, bcv.df = bcv.df,
                           dropout.type = dropout.type,
                           dropout.mid = dropout.mid,
                           dropout.shape = dropout.shape,
                           path.from = path.from, path.nSteps = path.nSteps,
                           path.skew = path.skew,
                           path.nonlinearProb = path.nonlinearProb,
                           path.sigmaFac = path.sigmaFac, seed = seed)
  sim <- splatSimulateGroups(params, verbose = FALSE)
  assays(sim) <- assays(sim)["counts"]
  size_factors <- centerSizeFactors(colSums(counts(sim)))
  normalized <- normalizeCounts(DelayedArray(assay(sim)), size_factors)
  assay(sim, "logcounts") <- as(normalized, "dgCMatrix")

  gene_var <- modelGeneVariances(normalized, num.threads = 2)
  top_hvgs <- chooseHighlyVariableGenes(gene_var$statistics$residuals,
                                        top = n_hvgs)

  sim <- sim[top_hvgs, ]
  sim <- runPCA(sim, ntop = n_hvgs, ncomponent = ncomp)
  return(sim)
}
