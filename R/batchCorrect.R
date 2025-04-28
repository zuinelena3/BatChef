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
