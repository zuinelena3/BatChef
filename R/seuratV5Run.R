#' Seurat V5 method
#'
#' SeuratV5 is an anchor-based method.
#'
#' @param input A \linkS4class{Seurat} object.
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
#' @param preserve_order Do not reorder objects based on size for each
#' pairwise integration.
#' @param verbose Print progress bars and output.
#' @param l2_norm Perform L2 normalization on the CCA cell embeddings after
#'  dimensional reduction.
#' @param k_anchor Number of neighbors (k) to use when picking anchors.
#' @param k_score Number of neighbors (k) to use when scoring anchors.
#' @param max_features The maximum number of features to use when specifying
#' the neighborhood search space in the anchor filtering.
#' @param nn_method Method for nearest neighbor finding.
#' @param n_trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#'
#' @export
#' @importFrom Seurat IntegrateLayers
#'
#' @return A \linkS4class{Seurat} object.
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(250, 200),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE,
#'                      output_format = "Seurat")
#' sim[[SeuratObject::DefaultAssay(sim)]] <- split(x = sim[[SeuratObject::DefaultAssay(sim)]],
#'                                                 f = sim[["Batch"]][, 1])
#' sim <- Seurat::ScaleData(sim, verbose = FALSE)
#' seuv5 <- seuratV5Run(input = sim, method = "CCAIntegration",
#'                      features = rownames(sim))
#'
seuratV5Run <- function(input, method = "CCAIntegration", orig_reduction = "pca",
                        assay = NULL, features = NULL, layers = NULL,
                        scale_layer = "scale.data",
                        new_reduction = "integrated.dr", reference = NULL,
                        normalization_method = "LogNormalize",
                        dims = 1:30, k_filter = NA, dims_to_integrate = NULL,
                        k_weight = 100, weight_reduction = NULL,
                        sd_weight = 1, sample_tree = NULL, preserve_order = FALSE,
                        verbose = FALSE, l2_norm = TRUE, k_anchor = 5,
                        k_score = 30, max_features = 200, nn_method = "annoy",
                        n_trees = 50, eps = 0) {

  suppressWarningsByMsg(
    c("deprecate", "without the associated assay"),
    IntegrateLayers(object = input, method = method, orig.reduction = orig_reduction,
                    assay = assay, features = features, layers = layers,
                    scale.layer = scale_layer,
                    new.reduction = new_reduction, reference = reference,
                    normalization.method = normalization_method,
                    dims = dims, k.filter = k_filter,
                    dims.to.integrate = dims_to_integrate,
                    k.weight = k_weight, weight.reduction = weight_reduction,
                    sd.weight = sd_weight, sample.tree = sample_tree,
                    preserve.order = preserve_order, verbose = verbose,
                    l2.norm = l2_norm, k.anchor = k_anchor,
                    k.score = k_score, max.features = max_features,
                    nn.method = nn_method, n.trees = n_trees, eps = eps)
  )
}
