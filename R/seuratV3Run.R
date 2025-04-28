#' Seurat V3 method
#'
#' @param input A list of Seurat objects.
#' @param assay A vector of assay names specifying which assay to use when constructing anchors. If NULL, the current default assay for each object is used.
#' @param reference A vector specifying the object/s to be used as a reference during integration.
#' @param anchor.features Number of features to be used in anchor finding.
#' @param scale A logical
#' @param normalization.method Name of normalization method used: LogNormalize (default) or SCT.
#' @param sct.clip.range sct.clip.range
#' @param reduction reduction
#' @param l2.norm l2.norm
#' @param dims dims
#' @param k.anchor k.anchor
#' @param k.filter k.filter
#' @param k.score k.score
#' @param max.features max.features
#' @param nn.method nn.method
#' @param n.trees n.trees
#' @param eps eps
#' @param verbose verbose
#' @param new.assay.name new.assay.name
#' @param features features
#' @param features.to.integrate features.to.integrate
#' @param k.weight k.weight
#' @param weight.reduction weight.reduction
#' @param sd.weight sd.weight
#' @param sample.tree sample.tree
#' @param preserve.order preserve.order
#'
#' @export
#' @importFrom Seurat FindIntegrationAnchors IntegrateData
seuratV3Run <- function(input, assay = NULL, reference = NULL, anchor.features = 2000, scale = TRUE, normalization.method = "LogNormalize",
                        sct.clip.range = NULL, reduction = "cca", l2.norm = TRUE, dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30,
                        max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0, verbose = TRUE, new.assay.name = "integrated",
                        features = NULL, features.to.integrate = NULL, k.weight = 100, weight.reduction = NULL, sd.weight = 1, sample.tree = NULL,
                        preserve.order = FALSE) {

  anchorset <- FindIntegrationAnchors(object.list = input, assay = assay, reference = reference, anchor.features = anchor.features, scale = scale,
                                      normalization.method = normalization.method,
                                      sct.clip.range = sct.clip.range, reduction = reduction, l2.norm = l2.norm, dims = dims, k.anchor = k.anchor, k.filter = k.filter,
                                      k.score = k.score,
                                      max.features = max.features, nn.method = nn.method, n.trees = n.trees, eps = eps, verbose = verbose)
  out <- IntegrateData(anchorset = anchorset, new.assay.name = new.assay.name,
                       normalization.method = normalization.method, features = features, features.to.integrate = features.to.integrate, dims = dims,
                       k.weight = k.weight, weight.reduction = weight.reduction, sd.weight = sd.weight, sample.tree = sample.tree,
                       preserve.order = preserve.order, eps = eps, verbose = verbose)
}
