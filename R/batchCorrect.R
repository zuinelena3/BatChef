#' Batch correction methods
#'
#' A common interface for single-cell batch correction methods.
#'
#' @param input A SingleCellExperiment, Seurat or AnnData objects can be supplied.
#' @param batch A string specifying the batch for each cell.
#' @param params A \linkS4class{BatChefParams} object specifying the batch correction method.
#'
#' @import methods
#' @rdname batchCorrect
setMethod("batchCorrect", "LimmaParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch, assay_type = params@assay_type), params@extra, "LimmaParams")
  out <- do.call(limmaRun, args)

  res <- linearPost(input = input, output = out, method = "limma")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "CombatParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch, assay_type = params@assay_type), params@extra, "CombatParams")
  out <- do.call(combatRun, args)

  res <- linearPost(input = input, output = out, method = "combat")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "SeuratV3Params", function(input, batch, params) {
  options(Seurat.object.assay.version = "v3")
  ll <- seuratv3Input(input = input, batch = batch, features = params@features, pca_name = params@pca_name)

  out <- seuratV3Run(input = ll, assay = params@assay, reference = params@reference, anchor.features = params@anchor.features,
                     scale = params@scale, normalization.method = params@normalization.method,
                     sct.clip.range = params@sct.clip.range, reduction = params@reduction, l2.norm = params@l2.norm, dims = params@dims,
                     k.anchor = params@k.anchor, k.filter = params@k.filter, k.score = params@k.score,
                     max.features = params@max.features, nn.method = params@nn.method, n.trees = params@n.trees, eps = params@eps,
                     verbose = params@verbose, new.assay.name = params@new.assay.name,
                     features = params@features, features.to.integrate = params@features.to.integrate, k.weight = params@k.weight,
                     weight.reduction = params@weight.reduction, sd.weight = params@sd.weight, sample.tree = params@sample.tree,
                     preserve.order = params@preserve.order)

  res <- seuratv3Post(input = input, output = out, method = paste0("seuratv3_", params@reduction))
})

#' @rdname batchCorrect
setMethod("batchCorrect", "SeuratV5Params", function(input, batch, params) {
  so <- seuratv5Input(input = input, batch = batch, pca_name = params@pca_name)

  out <- seuratV5Run(input = so, method = params@method, orig.reduction = params@orig.reduction, assay = params@assay,
                     features = params@features, layers = params@layers, scale.layer = params@scale.layer,
                     new.reduction = params@new.reduction, reference = params@reference, normalization.method = params@normalization.method,
                     dims = params@dims, k.filter = params@k.filter, dims.to.integrate = params@dims.to.integrate,
                     k.weight = params@k.weight, weight.reduction = params@weight.reduction,
                     sd.weight = params@sd.weight, sample.tree = params@sample.tree, preserve.order = params@preserve.order,
                     verbose = params@verbose, l2.norm = params@l2.norm, k.anchor = params@k.anchor,
                     k.score = params@k.score, max.features = params@max.features, nn.method = params@nn.method, n.trees = params@n.trees,
                     eps = params@eps)

  res <- seuratv5Post(input = input, output = out, method = paste0("seuratv5_", tolower(sub("\\I.*", "", params@method))), name = params@new.reduction)
})

#' @rdname batchCorrect
setMethod("batchCorrect", "FastMNNParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch), params@extra, "FastMNNParams")
  out <- do.call(fastMNNRun, args)

  res <- fastMNNPost(input = input, output = out, method = "fastmnn")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "HarmonyParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch), params@extra, "HarmonyParams")
  out <- do.call(harmonyRun, args)

  res <- harmonyPost(input = input, output = out, method = "harmony")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "ScanoramaParams", function(input, batch, params) {
  ll <- scanoramaInput(input = input, batch = batch, assay.type = params@assay_type)

  args <- merge_params(list(input = ll, return_dimred = params@return_dimred), params@extra, "ScanoramaParams")
  out <- do.call(scanoramaRun, args)

  res <- scanoramaPost(input = input, list = ll, output = out, return_dimred = params@return_dimred, method = "scanorama")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "SCVIParams", function(input, batch, params) {
  out <- scVIRun(input = input, batch = batch, layer = params@layer, labels_key = params@labels_key, size_factor_key = params@size_factor_key,
                 categorical_covariate_keys = params@categorical_covariate_keys, continuous_covariate_keys = params@continuous_covariate_keys,
                 n_hidden = params@n_hidden, n_latent = params@n_latent, n_layers = params@n_layers, dropout_rate = params@dropout_rate,
                 dispersion = params@dispersion, gene_likelihood = params@gene_likelihood, latent_distribution = params@latent_distribution,
                 max_epochs = params@max_epochs, accelerator = params@accelerator, devices = params@devices, train_size = params@train_size,
                 validation_size = params@validation_size, shuffle_set_split = params@shuffle_set_split,
                 load_sparse_tensor = params@load_sparse_tensor, batch_size = params@batch_size, early_stopping = params@early_stopping,
                 datasplitter_kwargs = params@datasplitter_kwargs, plan_kwargs = params@plan_kwargs, datamodule = params@datamodule,
                 indices = params@indices, give_mean = params@give_mean, mc_samples = params@mc_samples, return_dist = params@return_dist,
                 dataloader = params@dataloader)

  res <- scVIPost(input = input, output = out, method = "scvi")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "ScMerge2Params", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch, assay_type = params@assay_type), params@extra, "ScMerge2Params")
  out <- do.call(scMerge2Run, args)

  res <- scMerge2Post(input = input, output = out, method = "scmerge2")
})


#' @rdname batchCorrect
setMethod("batchCorrect", "BBKNNParams", function(input, batch, params) {
  ll <- bbknnInput(input = input, batch = batch, reduction = params@reduction)

  args <- merge_params(list(input = ll$red, batch = ll$batch), params@extra, "BBKNNParams")
  out <- do.call(bbknnRun, args)

  res <- bbknnPost(input = input, output = out, method = "bbknn")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "LigerParams", function(input, batch, params) {
  lo <- ligerInput(input = input, batch = batch, features = params@features, useDatasets = params@useDatasets,
                   format.type = params@format.type, verbose = params@verbose, remove.missing = params@remove.missing, params@extra_input)

  args <- merge_params(list(input = lo, method = params@method), params@extra, "LigerParams")
  out <- do.call(ligerRun, args)

  res <- ligerPost(input = input, output = out, method = "liger")
})
