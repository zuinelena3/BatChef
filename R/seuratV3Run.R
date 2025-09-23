#' Seurat V3 method
#'
#' SeuratV3 is an anchor-based method.
#'
#' @param input A list of \linkS4class{Seurat} objects.
#' @param assay A vector of assay names specifying which assay to use when
#' constructing anchors.
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration.
#' @param anchor.features Number of features to be used in anchor finding.
#' @param scale A logical to scale the features provided.
#' @param normalization.method Name of normalization method used:
#' LogNormalize (default) or SCT.
#' @param sct.clip.range Numeric of length two specifying the min and max values
#'  the Pearson residual will be clipped to.
#' @param reduction Dimensional reduction to perform when finding anchors.
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings
#' after dimensional reduction.
#' @param dims Number of dimensions.
#' @param k.anchor Number of neighbors (k) to use when picking anchors.
#' @param k.filter Number of neighbors (k) to use when filtering anchors.
#' @param k.score Number of neighbors (k) to use when scoring anchors.
#' @param max.features The maximum number of features to use when specifying
#' the neighborhood search space in the anchor filtering.
#' @param nn.method Method for nearest neighbor finding.
#' @param n.trees More trees gives higher precision when using annoy
#' approximate nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#' @param verbose Print progress bars and output.
#' @param new.assay.name Name for the new assay containing
#' the integrated data.
#' @param features Vector of features to use.
#' @param features.to.integrate Vector of features to integrate.
#' @param k.weight Number of neighbors to consider when weighting anchors.
#' @param weight.reduction Dimension reduction to use when calculating
#' anchor weights.
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample.tree Specify the order of integration.
#' @param preserve.order Do not reorder objects based on size for each
#' pairwise integration.
#'
#' @export
#' @importFrom Seurat FindIntegrationAnchors IntegrateData
#' @return A \linkS4class{Seurat} object that contains the corrected matrix.#'
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' so <- as.Seurat(sim)
#' VariableFeatures(so) <- rownames(so)
#' so[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(sim, "PCA"),
#'                                     loadings = attr(reducedDim(sim, "PCA"), "rotation"),
#'                                     key = "pca_", assay = DefaultAssay(so))
#' ll <- SplitObject(so, split.by = "Batch")
#' seuv3 <- seuratV3Run(input = ll, reduction = "cca")
#'
seuratV3Run <- function(input, assay = NULL, reference = NULL,
                        anchor.features = 2000, scale = TRUE,
                        normalization.method = "LogNormalize",
                        sct.clip.range = NULL, reduction = "cca", l2.norm = TRUE,
                        dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30,
                        max.features = 200, nn.method = "annoy", n.trees = 50,
                        eps = 0, verbose = TRUE, new.assay.name = "integrated",
                        features = NULL, features.to.integrate = NULL,
                        k.weight = 100, weight.reduction = NULL, sd.weight = 1,
                        sample.tree = NULL, preserve.order = FALSE) {

  anchorset <- FindIntegrationAnchors(object.list = input, assay = assay,
                                      reference = reference,
                                      anchor.features = anchor.features,
                                      scale = scale,
                                      normalization.method = normalization.method,
                                      sct.clip.range = sct.clip.range,
                                      reduction = reduction, l2.norm = l2.norm,
                                      dims = dims, k.anchor = k.anchor,
                                      k.filter = k.filter,
                                      k.score = k.score,
                                      max.features = max.features,
                                      nn.method = nn.method, n.trees = n.trees,
                                      eps = eps, verbose = verbose)
  out <- IntegrateData(anchorset = anchorset, new.assay.name = new.assay.name,
                       normalization.method = normalization.method,
                       features = features,
                       features.to.integrate = features.to.integrate,
                       dims = dims, k.weight = k.weight,
                       weight.reduction = weight.reduction,
                       sd.weight = sd.weight, sample.tree = sample.tree,
                       preserve.order = preserve.order, eps = eps,
                       verbose = verbose)
}
