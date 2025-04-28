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

  res <- seuratv3Post(input = input, output = out, method = "seuratv3")
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

  res <- seuratv5Post(input = input, output = out)
})
