#' BatChefParams class
#'
#' @export
#' @importClassesFrom S4Vectors SimpleList
#' @import methods
#'
#' @rdname BatChefParams
setClass("BatChefParams", contains = "SimpleList")

#' @export
#' @rdname BatChefParams
setClass("LimmaParams", contains = "BatChefParams", slots = c(assay_type = "character", extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("CombatParams", contains = "BatChefParams", slots = c(assay_type = "character", extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("SeuratV3Params", contains = "BatChefParams", slots = c(features = "character", pca_name = "ANY", assay = "ANY", reference = "ANY",
                                                                 anchor.features = "numeric", scale = "logical",
                                                                 normalization.method = "character",
                                                                 sct.clip.range = "ANY", reduction = "character",
                                                                 l2.norm = "logical", dims = "numeric", k.anchor = "numeric",
                                                                 k.filter = "numeric", k.score = "numeric",
                                                                 max.features = "numeric", nn.method = "character",
                                                                 n.trees = "numeric", eps = "numeric", verbose = "logical",
                                                                 new.assay.name = "character",
                                                                 features.to.integrate = "ANY", k.weight = "numeric",
                                                                 weight.reduction = "ANY", sd.weight = "numeric", sample.tree = "ANY",
                                                                 preserve.order = "logical"))

#' @export
#' @rdname BatChefParams
setClass("SeuratV5Params", contains = "BatChefParams", slots = c(pca_name = "ANY", method = "character", orig.reduction = "character",
                                                                 assay = "ANY", features = "ANY", layers = "ANY", scale.layer = "character",
                                                                 new.reduction = "character", reference = "ANY",
                                                                 normalization.method = "character", dims = "numeric",
                                                                 k.filter = "ANY", dims.to.integrate = "ANY", k.weight = "numeric",
                                                                 weight.reduction = "ANY", sd.weight = "numeric", sample.tree = "ANY",
                                                                 preserve.order = "logical", verbose = "logical",
                                                                 l2.norm = "logical", k.anchor = "numeric", k.score = "numeric",
                                                                 max.features = "numeric", nn.method = "character",
                                                                 n.trees = "numeric", eps = "numeric"))

#' @export
#' @rdname BatChefParams
setClass("FastMNNParams", contains = "BatChefParams", slots = c(extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("HarmonyParams", contains = "BatChefParams", slots = c(extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("ScanoramaParams", contains = "BatChefParams", slots = c(assay_type = "ANY",
                                                                  return_dimred = "logical", extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("SCVIParams", contains = "BatChefParams", slots = c(layer = "ANY", labels_key = "ANY", size_factor_key = "ANY", categorical_covariate_keys = "ANY",
                                                             continuous_covariate_keys = "ANY",
                                                             n_hidden = "numeric", n_latent = "numeric", n_layers = "numeric", dropout_rate = "numeric",
                                                             dispersion = "character", gene_likelihood = "character", latent_distribution = "character",
                                                             max_epochs = "numeric", accelerator = "character", devices = "numeric", train_size = "numeric", validation_size = "ANY", shuffle_set_split = "logical",
                                                             load_sparse_tensor = "logical", batch_size = "numeric", early_stopping = "logical", datasplitter_kwargs = "ANY", plan_kwargs = "ANY", datamodule = "ANY",
                                                             indices = "ANY", give_mean = "logical", mc_samples = "numeric", return_dist = "logical", dataloader = "ANY"))

