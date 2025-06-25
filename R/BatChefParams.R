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

#' @param features Vector of features to use.
#' @param pca_name A string specifying the PCA name.
#' @param assay A vector of assay names specifying which assay to use when constructing anchors.
#' @param reference A vector specifying the object/s to be used as a reference during integration.
#' @param anchor.features Number of features to be used in anchor finding.
#' @param scale A logical to scale the features provided.
#' @param normalization.method Name of normalization method used: LogNormalize (default) or SCT.
#' @param sct.clip.range Numeric of length two specifying the min and max values the Pearson residual will be clipped to.
#' @param reduction Dimensional reduction to perform when finding anchors.
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after dimensional reduction.
#' @param dims Number of dimensions.
#' @param k.anchor Number of neighbors (k) to use when picking anchors.
#' @param k.filter Number of neighbors (k) to use when filtering anchors.
#' @param k.score Number of neighbors (k) to use when scoring anchors.
#' @param max.features The maximum number of features to use when specifying the neighborhood search space in the anchor filtering.
#' @param nn.method Method for nearest neighbor finding.
#' @param n.trees More trees gives higher precision when using annoy approximate nearest neighbor search.
#' @param eps Error bound on the neighbor finding algorithm.
#' @param verbose Print progress bars and output.
#' @param new.assay.name Name for the new assay containing the integrated data.
#' @param features.to.integrate Vector of features to integrate.
#' @param k.weight Number of neighbors to consider when weighting anchors.
#' @param weight.reduction Dimension reduction to use when calculating anchor weights.
#' @param sd.weight Controls the bandwidth of the Gaussian kernel for weighting.
#' @param sample.tree Specify the order of integration.
#' @param preserve.order Do not reorder objects based on size for each pairwise integration.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
#'
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

#' @param pca_name  A string specifying the PCA name.
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
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
#'
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

#' @param assay_type assay_type
#' @param return_dimred return_dimred
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
ScanoramaParams <- function(assay_type = NULL, return_dimred = FALSE, ...) {
  new("ScanoramaParams", assay_type = assay_type, return_dimred = return_dimred, extra = list(...))
}

#' @param layer layer
#' @param labels_key labels_key
#' @param size_factor_key size_factor_key
#' @param categorical_covariate_keys categorical_covariate_keys
#' @param continuous_covariate_keys continuous_covariate_keys
#' @param n_hidden n_hidden
#' @param n_latent n_latent
#' @param n_layers n_layers
#' @param dropout_rate dropout_rate
#' @param dispersion dispersion
#' @param gene_likelihood gene_likelihood
#' @param latent_distribution latent_distribution
#' @param max_epochs max_epochs
#' @param accelerator accelerator
#' @param devices devices
#' @param train_size train_size
#' @param validation_size validation_size
#' @param shuffle_set_split shuffle_set_split
#' @param load_sparse_tensor load_sparse_tensor
#' @param batch_size batch_size
#' @param early_stopping early_stopping
#' @param datasplitter_kwargs datasplitter_kwargs
#' @param plan_kwargs plan_kwargs
#' @param datamodule datamodule
#' @param indices indices
#' @param give_mean give_mean
#' @param mc_samples mc_samples
#' @param return_dist return_dist
#' @param dataloader dataloader
#' @param transform_batch transform_batch
#' @param gene_list gene_list
#' @param library_size library_size
#' @param n_samples n_samples
#' @param n_samples_overall n_samples_overall
#' @param weights weights
#' @param return_mean return_mean
#' @param return_numpy return_numpy
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
ScVIParams <- function(layer = NULL, labels_key = NULL, size_factor_key = NULL, categorical_covariate_keys = NULL, continuous_covariate_keys = NULL,
                       n_hidden = 128, n_latent = 10, n_layers = 1, dropout_rate = 0.1,
                       dispersion = "gene", gene_likelihood = "zinb", latent_distribution = "normal",
                       max_epochs = 400, accelerator = "auto", devices = 1, train_size = 0.25, validation_size = NULL, shuffle_set_split = TRUE,
                       load_sparse_tensor = FALSE, batch_size = 128, early_stopping = FALSE, datasplitter_kwargs = NULL, plan_kwargs = NULL, datamodule = NULL,
                       indices = NULL, give_mean = TRUE, mc_samples = 5000, return_dist = FALSE, dataloader = NULL,
                       transform_batch = NULL, gene_list = NULL, library_size = 1, n_samples = 1, n_samples_overall = NULL, weights = NULL, return_mean = TRUE, return_numpy = NULL) {
  new("SCVIParams", layer = layer, labels_key = labels_key, size_factor_key = size_factor_key, categorical_covariate_keys = categorical_covariate_keys, continuous_covariate_keys = continuous_covariate_keys,
      n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers, dropout_rate = dropout_rate,
      dispersion = dispersion, gene_likelihood = gene_likelihood, latent_distribution = latent_distribution,
      max_epochs = max_epochs, accelerator = accelerator, devices = devices, train_size = train_size, validation_size = validation_size, shuffle_set_split = shuffle_set_split,
      load_sparse_tensor = load_sparse_tensor, batch_size = batch_size, early_stopping = early_stopping, datasplitter_kwargs = datasplitter_kwargs, plan_kwargs = plan_kwargs, datamodule = datamodule,
      indices = indices, give_mean = give_mean, mc_samples = mc_samples, return_dist = return_dist, dataloader = dataloader,
      transform_batch = transform_batch, gene_list = gene_list, library_size = library_size, n_samples = n_samples, n_samples_overall = n_samples_overall, weights = weights,
      return_mean = return_mean, return_numpy = return_numpy)
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

#' @param reduction reduction
#' @param ... Named arguments to pass to individual methods upon dispatch

#' @export
#' @rdname BatChefParams
#' @importFrom methods new
BBKNNParams <- function(reduction, ...) {
  new("BBKNNParams", reduction = reduction, extra = list(...))
}

#' @param dataloader dataloader
#' @param useDatasets useDatasets
#' @param verbose verbose
#' @param format.type format.type
#' @param remove.missing remove.missing
#' @param method method
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
LigerParams <- function(features, useDatasets = NULL,
                        verbose = TRUE, format.type = NULL, remove.missing = NULL, method = "iNMF", ...) {
  new("LigerParams", features = features, useDatasets = useDatasets,
      verbose = verbose, format.type = format.type, remove.missing = remove.missing, extra_input = list(...), method = method, extra = list(...))
}
