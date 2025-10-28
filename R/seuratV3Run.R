#' Seurat V3 method
#'
#' SeuratV3 is an anchor-based method.
#'
#' @param input A list of \link[Seurat]{Seurat} objects.
#' @param assay A vector of assay names specifying which assay to use when
#' constructing anchors.
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration.
#' @param anchor_features Number of features to be used in anchor finding.
#' @param scale A logical to scale the features provided.
#' @param normalization_method Name of normalization method used:
#' LogNormalize (default) or SCT.
#' @param sct_clip_range Numeric of length two specifying the min and max values
#'  the Pearson residual will be clipped to.
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
#' @param new_assay_name Name for the new assay containing
#' the integrated data.
#' @param features Vector of features to use.
#' @param features_to_integrate Vector of features to integrate.
#' @param k_weight Number of neighbors to consider when weighting anchors.
#' @param weight_reduction Dimension reduction to use when calculating
#' anchor weights.
#' @param sd_weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample_tree Specify the order of integration.
#' @param preserve_order Do not reorder objects based on size for each
#' pairwise integration.
#'
#' @export
#' @importFrom Seurat FindIntegrationAnchors IntegrateData
#'
#' @return A \link[Seurat]{Seurat} object that contains the corrected matrix.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(200, 200),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = FALSE,
#'                      output_format = "Seurat")
#' feat <- Seurat::VariableFeatures(sim)
#' sim <- Seurat::SplitObject(sim, split.by = "Batch")
#' seuv3 <- seuratV3Run(input = sim, reduction = "cca",
#'                      features = feat)
#'
seuratV3Run <- function(input, assay = NULL, reference = NULL,
                        anchor_features = 2000, scale = TRUE,
                        normalization_method = "LogNormalize",
                        sct_clip_range = NULL, reduction = "cca", l2_norm = TRUE,
                        dims = 1:30, k_anchor = 5, k_filter = 200, k_score = 30,
                        max_features = 200, nn_method = "annoy", n_trees = 50,
                        eps = 0, verbose = FALSE, new_assay_name = "integrated",
                        features = NULL, features_to_integrate = NULL,
                        k_weight = 100, weight_reduction = NULL, sd_weight = 1,
                        sample_tree = NULL, preserve_order = FALSE) {

  # NOTE: until https://github.com/satijalab/seurat/issues/9850 is fixed.
  suppressWarningsByMsg("deprecated",
                        anchorset <- FindIntegrationAnchors(
                          object.list = input, assay = assay,
                          reference = reference, anchor.features = anchor_features,
                          scale = scale, normalization.method = normalization_method,
                          sct.clip.range = sct_clip_range, reduction = reduction,
                          l2.norm = l2_norm, dims = dims, k.anchor = k_anchor,
                          k.filter = k_filter, k.score = k_score,
                          max.features = max_features, nn.method = nn_method,
                          n.trees = n_trees, eps = eps, verbose = verbose)
  )

  # NOTE: until https://github.com/satijalab/seurat/issues/8938 is fixed.
  suppressWarningsByMsg(
    "counts isn't present",
    IntegrateData(anchorset = anchorset, new.assay.name = new_assay_name,
                  normalization.method = normalization_method,
                  features = features,
                  features.to.integrate = features_to_integrate,
                  dims = dims, k.weight = k_weight,
                  weight.reduction = weight_reduction,
                  sd.weight = sd_weight, sample.tree = sample_tree,
                  preserve.order = preserve_order, eps = eps,
                  verbose = verbose)
  )
}
