#' Seurat V5 method
#'
#' @param input input
#' @param method method
#' @param orig.reduction orig.reduction
#' @param assay assay
#' @param features features
#' @param layers layers
#' @param scale.layer scale.layer
#' @param new.reduction new.reduction
#' @param reference reference
#' @param normalization.method normalization.method
#' @param dims dims
#' @param k.filter k.filter
#' @param dims.to.integrate dims.to.integrate
#' @param k.weight k.weight
#' @param weight.reduction weight.reduction
#' @param sd.weight sd.weight
#' @param sample.tree sample.tree
#' @param preserve.order preserve.order
#' @param verbose verbose
#' @param l2.norm l2.norm
#' @param k.anchor k.anchor
#' @param k.score k.score
#' @param max.features max.features
#' @param nn.method nn.method
#' @param n.trees n.trees
#' @param eps eps
#'
#' @export
#' @importFrom Seurat IntegrateLayers
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
