#' Batch correction methods
#'
#' A common interface for single-cell batch correction methods.
#'
#' Users can pass parameters to each method via the constructors for params.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param params A \linkS4class{BatChefParams} object specifying
#' the batch correction method to use and the parameters for its execution.
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object, where the output of the method
#' (such as the corrected gene expression matrix and/or the corrected
#' dimensional reduction space) is stored within the original input object.
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' out <- batchCorrect(input = sim, batch = "Batch", params = HarmonyParams())
#' @rdname batchCorrect
#'
setMethod("batchCorrect", "LimmaParams", function(input, batch, params) {
  sce <- linearInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch,
                            assay_type = params@assay_type),
                       params@extra, "LimmaParams")
  out <- do.call(limmaRun, args)

  res <- linearPost(input = input, output = out, method = "limma")
})

#' @rdname batchCorrect
#'
setMethod("batchCorrect", "CombatParams", function(input, batch, params) {
  sce <- linearInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch,
                            assay_type = params@assay_type),
                       params@extra, "CombatParams")
  out <- do.call(combatRun, args)

  res <- linearPost(input = input, output = out, method = "combat")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "SeuratV3Params", function(input, batch, params) {
  options(Seurat.object.assay.version = "v3")
  ll <- seuratv3Input(input = input, batch = batch, features = params@features,
                      pca_name = params@pca_name)

  out <- seuratV3Run(input = ll, assay = params@assay,
                     reference = params@reference,
                     anchor_features = params@anchor_features,
                     scale = params@scale,
                     normalization_method = params@normalization_method,
                     sct_clip_range = params@sct_clip_range,
                     reduction = params@reduction,
                     l2_norm = params@l2_norm, dims = params@dims,
                     k_anchor = params@k_anchor, k_filter = params@k_filter,
                     k_score = params@k_score,
                     max_features = params@max_features,
                     nn_method = params@nn_method, n_trees = params@n_trees,
                     eps = params@eps, verbose = params@verbose,
                     new_assay_name = params@new_assay_name,
                     features = params@features,
                     features_to_integrate = params@features_to_integrate,
                     k_weight = params@k_weight,
                     weight_reduction = params@weight_reduction,
                     sd_weight = params@sd_weight,
                     sample_tree = params@sample_tree,
                     preserve_order = params@preserve_order)

  res <- seuratv3Post(input = input, output = out,
                      method = paste0("seuratv3_", params@reduction))
})

#' @rdname batchCorrect
setMethod("batchCorrect", "SeuratV5Params", function(input, batch, params) {
  so <- seuratv5Input(input = input, batch = batch, features = params@features,
                      pca_name = params@pca_name)

  out <- seuratV5Run(input = so, method = params@method,
                     orig_reduction = params@orig_reduction,
                     assay = params@assay,
                     features = params@features, layers = params@layers,
                     scale_layer = params@scale_layer,
                     new_reduction = params@new_reduction,
                     reference = params@reference,
                     normalization_method = params@normalization_method,
                     dims = params@dims, k_filter = params@k_filter,
                     dims_to_integrate = params@dims_to_integrate,
                     k_weight = params@k_weight,
                     weight_reduction = params@weight_reduction,
                     sd_weight = params@sd_weight,
                     sample_tree = params@sample_tree,
                     preserve_order = params@preserve_order,
                     verbose = params@verbose, l2_norm = params@l2_norm,
                     k_anchor = params@k_anchor, k_score = params@k_score,
                     max_features = params@max_features,
                     nn_method = params@nn_method,
                     n_trees = params@n_trees, eps = params@eps)

  res <- seuratv5Post(input = input, output = out,
                      method = paste0("seuratv5_", tolower(sub("\\I.*", "", params@method))),
                      name = params@new_reduction)
})

#' @rdname batchCorrect
setMethod("batchCorrect", "FastMNNParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch),
                       params@extra, "FastMNNParams")
  out <- do.call(fastMNNRun, args)

  res <- fastMNNPost(input = input, output = out, method = "fastmnn")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "HarmonyParams", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch),
                       params@extra, "HarmonyParams")
  out <- do.call(harmonyRun, args)

  res <- harmonyPost(input = input, output = out, method = "harmony")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "ScMerge2Params", function(input, batch, params) {
  sce <- sceInput(input = input, batch = batch)

  args <- merge_params(list(input = sce, batch = batch,
                            assay_type = params@assay_type),
                       params@extra, "ScMerge2Params")
  out <- do.call(scMerge2Run, args)

  res <- scMerge2Post(input = input, output = out, method = "scmerge2")
})

#' @rdname batchCorrect
setMethod("batchCorrect", "LigerParams", function(input, batch, params) {
  lo <- ligerInput(input = input, batch = batch, features = params@features)

  args <- merge_params(list(input = lo, method = params@method),
                       params@extra, "LigerParams")
  out <- do.call(ligerRun, args)

  res <- ligerPost(input = input, output = out, method = "liger")
})
