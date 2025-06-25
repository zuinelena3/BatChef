#' Seurat V5 method
#'
#' @param input A `Seurat` object.
#' @param method Integration method function.
#' @param orig.reduction Name of dimensional reduction for correction.
#' @param assay Name of assay for integration.
#' @param features A vector of features to use for integration.
#' @param layers Names of normalized layers in assay.
#' @param scale.layer Name(s) of scaled layer(s) in assay.
#' @param new.reduction Name of new integrated dimensional reduction.
#' @param reference A reference Seurat object.
#' @param normalization.method Name of normalization method used: LogNormalize or SCT.
#' @param dims Number of dimensions of dimensional reduction.
#' @param k.filter Number of anchors to filter.
#' @param dims.to.integrate Number of dimensions to return integrated values for.
#' @param k.weight Number of neighbors to consider when weighting anchors.
#' @param weight.reduction Dimension reduction to use when calculating anchor weights.
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample.tree Specify the order of integration.
#' @param preserve.order Do not reorder objects based on size for each pairwise integration.
#' @param verbose Print progress bars and output.
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after dimensional reduction.
#' @param k.anchor Number of neighbors (k) to use when picking anchors.
#' @param k.score Number of neighbors (k) to use when scoring anchors.
#' @param max.features The maximum number of features to use when specifying the neighborhood search space in the anchor filtering.
#' @param nn.method Method for nearest neighbor finding.
#' @param n.trees More trees gives higher precision when using annoy approximate nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#' @export
#' @importFrom Seurat IntegrateLayers
#'
seuratV5Run <- function(input, method = "CCAIntegration", orig.reduction = "pca", assay = NULL, features = NULL, layers = NULL, scale.layer = "scale.data",
                        new.reduction = "integrated.dr", reference = NULL, normalization.method = "LogNormalize",
                        dims = 1:30, k.filter = NA, dims.to.integrate = NULL, k.weight = 100, weight.reduction = NULL,
                        sd.weight = 1, sample.tree = NULL, preserve.order = FALSE, verbose = TRUE, l2.norm = TRUE, k.anchor = 5,
                        k.score = 30, max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0) {

  IntegrateLayers(object = input, method = method, orig.reduction = orig.reduction, assay = assay, features = features, layers = layers, scale.layer = scale.layer,
                  new.reduction = new.reduction, reference = reference, normalization.method = normalization.method,
                  dims = dims, k.filter = k.filter, dims.to.integrate = dims.to.integrate, k.weight = k.weight, weight.reduction = weight.reduction,
                  sd.weight = sd.weight, sample.tree = sample.tree, preserve.order = preserve.order, verbose = verbose, l2.norm = l2.norm, k.anchor = k.anchor,
                  k.score = k.score, max.features = max.features, nn.method = nn.method, n.trees = n.trees, eps = eps)
}
