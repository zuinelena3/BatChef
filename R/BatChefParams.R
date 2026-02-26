#' BatChefParams methods
#'
#' Constructors and methods for the params parameter classes.
#' BatChefParams objects contain method specific parameters
#' to pass to the batchCorrect generic.
#'
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @return
#' A BatChefParams object of the specified subclass, containing parameter
#' settings for the corresponding batch correction method.
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

#' @param features Vector of features to use.
#' @param pca_name A string specifying the Principal Component Analysis name.
#' @param assay A vector of assay names specifying which assay to use
#' when constructing anchors.
#' @param reference A vector specifying the object/s to be used
#' as a reference during integration.
#' @param anchor_features Number of features to be used in anchor finding.
#' @param scale A logical to scale the features provided.
#' @param normalization_method Name of normalization method used:
#' LogNormalize (default) or SCT.
#' @param sct_clip_range Numeric of length two specifying the min and max values
#' the Pearson residual will be clipped to.
#' @param reduction Dimensional reduction to perform when finding anchors.
#' @param l2_norm Perform L2 normalization on the CCA cell embeddings
#' after dimensional reduction.
#' @param dims Number of dimensions.
#' @param k_anchor Number of neighbors (k) to use when picking anchors.
#' @param k_filter Number of neighbors (k) to use when filtering anchors.
#' @param k_score Number of neighbors (k) to use when scoring anchors.
#' @param max_features The maximum number of features to use when specifying
#' the neighborhood search space in the anchor filtering.
#' @param nn_method Method for nearest neighbor finding.
#' @param n_trees More trees gives higher precision when using annoy
#' approximate nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#' @param verbose Print progress bars and output.
#' @param new_assay_name Name for the new assay containing the integrated data.
#' @param features_to_integrate Vector of features to integrate.
#' @param k_weight Number of neighbors to consider when weighting anchors.
#' @param weight_reduction Dimension reduction to use when calculating anchor
#'  weights.
#' @param sd_weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample_tree Specify the order of integration.
#' @param preserve_order Do not reorder objects based on size for each pairwise
#' integration.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
#'
SeuratV3Params <- function(features, pca_name = NULL, assay = NULL,
                           reference = NULL, anchor_features = 2000,
                           scale = TRUE, normalization_method = "LogNormalize",
                           sct_clip_range = NULL, reduction = "cca",
                           l2_norm = TRUE, dims = 1:30, k_anchor = 5,
                           k_filter = 200, k_score = 30,
                           max_features = 200, nn_method = "annoy", n_trees = 50,
                           eps = 0, verbose = FALSE, new_assay_name = "integrated",
                           features_to_integrate = NULL, k_weight = 100,
                           weight_reduction = NULL, sd_weight = 1,
                           sample_tree = NULL, preserve_order = FALSE) {
  new("SeuratV3Params",
    features = features, pca_name = pca_name, assay = assay,
    reference = reference, anchor_features = anchor_features, scale = scale,
    normalization_method = normalization_method,
    sct_clip_range = sct_clip_range, reduction = reduction, l2_norm = l2_norm,
    dims = dims, k_anchor = k_anchor, k_filter = k_filter, k_score = k_score,
    max_features = max_features, nn_method = nn_method, n_trees = n_trees,
    eps = eps, verbose = verbose, new_assay_name = new_assay_name,
    features_to_integrate = features_to_integrate, k_weight = k_weight,
    weight_reduction = weight_reduction, sd_weight = sd_weight,
    sample_tree = sample_tree, preserve_order = preserve_order
  )
}

#' @param pca_name  A string specifying the PCA name.
#' @param method Integration method function.
#' @param orig_reduction Name of dimensional reduction for correction.
#' @param assay Name of assay for integration.
#' @param features A vector of features to use for integration.
#' @param layers Names of normalized layers in assay.
#' @param scale_layer Name(s) of scaled layer(s) in assay.
#' @param new_reduction Name of new integrated dimensional reduction.
#' @param reference A reference Seurat object.
#' @param normalization_method Name of normalization method used:
#' LogNormalize or SCT.
#' @param dims Number of dimensions of dimensional reduction.
#' @param k_filter Number of anchors to filter.
#' @param dims_to_integrate Number of dimensions to return integrated values for.
#' @param k_weight Number of neighbors to consider when weighting anchors.
#' @param weight_reduction Dimension reduction to use when calculating
#' anchor weights.
#' @param sd_weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample_tree Specify the order of integration.
#' @param preserve_order Do not reorder objects based on size for each pairwise
#' integration.
#' @param verbose Print progress bars and output.
#' @param l2_norm Perform L2 normalization on the CCA cell embeddings
#' after dimensional reduction.
#' @param k_anchor Number of neighbors (k) to use when picking anchors.
#' @param k_score Number of neighbors (k) to use when scoring anchors.
#' @param max_features The maximum number of features to use when specifying
#' the neighborhood search space in the anchor filtering.
#' @param nn_method Method for nearest neighbor finding.
#' @param n_trees More trees gives higher precision when using annoy
#' approximate nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
#'
SeuratV5Params <- function(pca_name = NULL, method = "CCAIntegration",
                           orig_reduction = "pca", assay = NULL,
                           features = NULL, layers = NULL,
                           scale_layer = "scale.data",
                           new_reduction = "integrated.dr", reference = NULL,
                           normalization_method = "LogNormalize",
                           dims = 1:30, k_filter = NA, dims_to_integrate = NULL,
                           k_weight = 100, weight_reduction = NULL,
                           sd_weight = 1, sample_tree = NULL,
                           preserve_order = FALSE, verbose = FALSE,
                           l2_norm = TRUE, k_anchor = 5,
                           k_score = 30, max_features = 200, nn_method = "annoy",
                           n_trees = 50, eps = 0) {
  new("SeuratV5Params",
    pca_name = pca_name, method = method,
    orig_reduction = orig_reduction, assay = assay, features = features,
    layers = layers, scale_layer = scale_layer,
    new_reduction = new_reduction, reference = reference,
    normalization_method = normalization_method,
    dims = dims, k_filter = k_filter, dims_to_integrate = dims_to_integrate,
    k_weight = k_weight, weight_reduction = weight_reduction,
    sd_weight = sd_weight, sample_tree = sample_tree,
    preserve_order = preserve_order,
    verbose = verbose, l2_norm = l2_norm, k_anchor = k_anchor,
    k_score = k_score, max_features = max_features, nn_method = nn_method,
    n_trees = n_trees, eps = eps
  )
}

#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
FastMNNParams <- function(...) {
  new("FastMNNParams", extra = list(...))
}

#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
HarmonyParams <- function(...) {
  new("HarmonyParams", extra = list(...))
}

#' @param assay_type A string specifying the assay to use for correction.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
ScMerge2Params <- function(assay_type = "logcounts", ...) {
  new("ScMerge2Params", assay_type = assay_type, extra = list(...))
}

#' @param features Vector of features to use.
#' @param method iNMF variant algorithm to use for integration.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
LigerParams <- function(features, method = "iNMF", ...) {
  new("LigerParams", features = features, method = method, extra = list(...))
}
