#' Function to simulate data
#'
#' The function allows to simulated single-cell RNA-seq data using Splatter
#' package, normalize data, select highly variable genes and
#' compute Principal Component Analysis.
#'
#' @param n_genes Number of genes.
#' @param batch_cells Number of cells per batch.
#' @param batch_fac_loc Batch factor location parameter.
#' @param batch_fac_scale Batch factor scale parameter.
#' @param batch_rm_effect Remove batch effect.
#' @param mean_rate Mean rate.
#' @param mean_shape Mean shape.
#' @param lib_loc Library size location parameter.
#' @param lib_scale Library size scale parameter.
#' @param lib_norm Library size distribution.
#' @param out_prob Expression outlier probability.
#' @param out_fac_loc Expression outlier factor location.
#' @param out_fac_scale Expression outlier factor scale.
#' @param group_prob Group probabilities.
#' @param de_prob Differential expression probability.
#' @param de_down_prob Down-regulation probability.
#' @param de_fac_loc DE factor location.
#' @param de_fac_scale DE factor scale.
#' @param bcv_common Common biological coefficient of variation.
#' @param bcv_df BCV degrees of freedom.
#' @param dropout_type Dropout type.
#' @param dropout_mid Dropout mid point.
#' @param dropout_shape Dropout shape.
#' @param path_from Path origin.
#' @param path_n_steps Number of steps.
#' @param path_skew Path skew.
#' @param path_nonlinear_prob Non-linear probability.
#' @param path_sigma_fac Path skew.
#' @param compute_hvgs Boolean value. If TRUE, highly variable genes
#' will be selected Default is FALSE.
#' @param n_hvgs Number of highly variable genes.
#' @param num_threads Integer scalar specifying the number of threads to use.
#' @param compute_pca Boolean value. If TRUE, Principal Component Analysis (PCA)
#' will be computed. Default is FALSE.
#' @param pca_ncomp Number of principal component.
#' @param output_format A \linkS4class{SingleCellExperiment},
#' \linkS4class{Seurat}, or `AnnData` object.
#' @param seed Random seed.
#'
#' @export
#' @importFrom splatter newSplatParams splatSimulateGroups
#' @importFrom SummarizedExperiment assays
#' @importFrom scrapper centerSizeFactors normalizeCounts modelGeneVariances chooseHighlyVariableGenes runPca
#' @importFrom Seurat CreateSeuratObject VariableFeatures<- CreateDimReducObject DefaultAssay
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom anndata AnnData
#' @importFrom SparseArray colSums
#' @importFrom S4Vectors SimpleList
#'
#' @return A \linkS4class{SingleCellExperiment},
#' \linkS4class{Seurat}, or `AnnData` object.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 1000,
#'                      compute_pca = FALSE)
#'
simulate_data <- function(n_genes = 10000, batch_cells = 100,
                          batch_fac_loc = 0.1, batch_fac_scale = 0.1,
                          batch_rm_effect = FALSE,
                          mean_rate = 0.3, mean_shape = 0.6,
                          lib_loc = 11, lib_scale = 0.2, lib_norm = FALSE,
                          out_prob = 0.05, out_fac_loc = 4, out_fac_scale = 0.5,
                          group_prob = 1, de_prob = 0.1, de_down_prob = 0.5,
                          de_fac_loc = 0.1, de_fac_scale = 0.4,
                          bcv_common = 0.1, bcv_df = 60, dropout_type = "none",
                          dropout_mid = 0, dropout_shape = -1,
                          path_from = 0, path_n_steps = 100, path_skew = 0.5,
                          path_nonlinear_prob = 0.1, path_sigma_fac = 0.8,
                          compute_hvgs = FALSE, n_hvgs = 1000, num_threads = 1,
                          compute_pca = FALSE, pca_ncomp = 10,
                          output_format = "SingleCellExperiment", seed = 333) {
  params <- newSplatParams(nGenes = n_genes, batchCells = batch_cells,
                           batch.facLoc = batch_fac_loc,
                           batch.facScale = batch_fac_scale,
                           batch.rmEffect = batch_rm_effect,
                           mean.rate = mean_rate, mean.shape = mean_shape,
                           lib.loc = lib_loc, lib.scale = lib_scale,
                           lib.norm = lib_norm, out.prob = out_prob,
                           out.facLoc = out_fac_loc, out.facScale = out_fac_scale,
                           group.prob = group_prob, de.prob = de_prob,
                           de.downProb = de_down_prob, de.facLoc = de_fac_loc,
                           de.facScale = de_fac_scale, bcv.common = bcv_common,
                           bcv.df = bcv_df, dropout.type = dropout_type,
                           dropout.mid = dropout_mid, dropout.shape = dropout_shape,
                           path.from = path_from, path.nSteps = path_n_steps,
                           path.skew = path_skew,
                           path.nonlinearProb = path_nonlinear_prob,
                           path.sigmaFac = path_sigma_fac, seed = seed)
  sim <- splatSimulateGroups(params, verbose = FALSE)

  counts <- as(assay(sim, "counts"), "dgCMatrix")
  coldata <- colData(sim)

  size_factors <- centerSizeFactors(colSums(counts))
  normalized <- as(normalizeCounts(counts, size_factors), "dgCMatrix")

  gene_var <- modelGeneVariances(normalized, num.threads = num_threads)
  top_hvgs <- chooseHighlyVariableGenes(gene_var$statistics$residuals,
                                        top = n_hvgs)

  rowdata <- data.frame(genes = rownames(normalized),
                        hvgs = seq_along(rownames(normalized)) %in% top_hvgs)

  pca <- runPca(normalized[top_hvgs, ], num.threads = num_threads, number = pca_ncomp)
  loadings <- pca[[2]]
  colnames(loadings) <- paste0("PC_", 1:pca_ncomp)
  rownames(loadings) <- rownames(sim)[top_hvgs]

  pca <- t(pca[[1]])
  colnames(pca) <- paste0("PC_", 1:pca_ncomp)
  rownames(pca) <- colnames(sim)

  if (output_format == "Seurat") {
    assay <- CreateAssay5Object(counts = counts, data = normalized)
    so <- CreateSeuratObject(assay, meta.data = as.data.frame(coldata))
    VariableFeatures(so) <- rowdata$genes[rowdata$hvgs]

    if (compute_pca == TRUE) {
      so[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                          loadings = loadings,
                                          key = "pca_", assay = DefaultAssay(so))
      return(so)
    }
    else return(so)
  }

  else if (output_format == "SingleCellExperiment") {
    sce <- SingleCellExperiment(assays = list(counts = counts, logcounts = normalized),
                                colData = coldata, rowData = rowdata)

    if (compute_pca == TRUE) {
      reducedDim(sce, "PCA") <- pca
      attr(reducedDim(sce, "PCA"), "rotation") <- loadings
      return(sce)
    }
    else return(sce)
  }

  else {
    counts <- counts[top_hvgs, ]
    normalized <- normalized[top_hvgs, ]
    obs <- as.data.frame(coldata)
    rownames(obs) <- colnames(counts)
    var <- rowdata[top_hvgs, ]
    rownames(var) <- rownames(counts)
    adata <- AnnData(X = t(as.matrix(counts)),
                     obs = obs,
                     var = var)
    adata$layers[["logcounts"]] <- t(as.matrix(normalized))

    if (compute_pca == TRUE) {
      adata$obsm[["X_pca"]] <- pca
      adata$varm[["PCs"]] <- loadings
      return(adata)
    }
    else return(adata)
  }
}
