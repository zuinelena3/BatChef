#' BatChefParams methods
#'
#' Constructors and methods for the params parameter classes.
#' BatChefParams objects contain method specific parameters to pass to the batchCorrect generic.
#'
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
LimmaParams <- function(assay_type = "logcounts", ...) {
  new("LimmaParams", assay_type = assay_type, extra = list(...))
}

#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
CombatParams <- function(assay_type = "counts", ...) {
  new("CombatParams", assay_type = assay_type, extra = list(...))
}

#' @param pca_name pca_name
#' @param assay assay
#' @param reference reference
#' @param anchor.features anchor.features
#' @param scale scale
#' @param normalization.method normalization.method
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
#' @rdname BatChefParams
#' @importFrom methods new
SeuratV3Params <- function(features, pca_name = NULL, assay = NULL, reference = NULL, anchor.features = 2000, scale = TRUE, normalization.method = "LogNormalize",
                           sct.clip.range = NULL, reduction = "cca", l2.norm = TRUE, dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30,
                           max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0, verbose = TRUE, new.assay.name = "integrated",
                           features.to.integrate = NULL, k.weight = 100, weight.reduction = NULL, sd.weight = 1, sample.tree = NULL,
                           preserve.order = FALSE) {
  new("SeuratV3Params", features = features, pca_name = pca_name, assay = assay, reference = reference, anchor.features = anchor.features, scale = scale, normalization.method = normalization.method,
      sct.clip.range = sct.clip.range, reduction = reduction, l2.norm = l2.norm, dims = dims, k.anchor = k.anchor, k.filter = k.filter, k.score = k.score,
      max.features = max.features, nn.method = nn.method, n.trees = n.trees, eps = eps, verbose = verbose, new.assay.name = new.assay.name,
      features.to.integrate = features.to.integrate, k.weight = k.weight, weight.reduction = weight.reduction, sd.weight = sd.weight, sample.tree = sample.tree,
      preserve.order = preserve.order)
}


#' @param pca_name pca_name
#' @param method method
#' @param orig.reduction orig.reduction
#' @param layers layers
#' @param scale.layer scale.layer
#' @param normalization.method normalization.method
#' @param new.reduction new.reduction
#' @param reference reference
#' @param dims dims
#' @param dims.to.integrate dims.to.integrate
#' @param k.anchor k.anchor
#' @param k.filter k.filter
#' @param k.score k.score
#' @param max.features max.features
#' @param nn.method nn.method
#' @param n.trees n.trees
#' @param eps eps
#' @param verbose verbose
#' @param l2.norm l2.norm
#' @param features features
#' @param k.weight k.weight
#' @param weight.reduction weight.reduction
#' @param sd.weight sd.weight
#' @param sample.tree sample.tree
#' @param preserve.order preserve.order
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
SeuratV5Params <- function(pca_name = NULL, method = "CCAIntegration", orig.reduction = "pca", assay = NULL, features = NULL, layers = NULL, scale.layer = "scale.data",
                           new.reduction = "integrated.dr", reference = NULL, normalization.method = "LogNormalize",
                           dims = 1:30, k.filter = NA, dims.to.integrate = NULL, k.weight = 100, weight.reduction = NULL,
                           sd.weight = 1, sample.tree = NULL, preserve.order = FALSE, verbose = TRUE, l2.norm = TRUE, k.anchor = 5,
                           k.score = 30, max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0) {
  new("SeuratV5Params", pca_name = pca_name, method = method, orig.reduction = orig.reduction, assay = assay, features = features, layers = layers, scale.layer = scale.layer,
      new.reduction = new.reduction, reference = reference, normalization.method = normalization.method,
      dims = dims, k.filter = k.filter, dims.to.integrate = dims.to.integrate, k.weight = k.weight, weight.reduction = weight.reduction,
      sd.weight = sd.weight, sample.tree = sample.tree, preserve.order = preserve.order, verbose = verbose, l2.norm = l2.norm, k.anchor = k.anchor,
      k.score = k.score, max.features = max.features, nn.method = nn.method, n.trees = n.trees, eps = eps)
}
