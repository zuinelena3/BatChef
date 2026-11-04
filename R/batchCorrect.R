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
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
# out <- batchCorrect(input = sim, batch = "Batch",
#                     params = ScanoramaParams(assay_type = "logcounts"))
#' @rdname batchCorrect
setMethod("batchCorrect", "ScanoramaParams", function(input, batch, params) {
  ll <- scanoramaInput(input = input, batch = batch, assay_type = params@assay_type)

  args <- merge_params(
    list(input = ll, return_dimred = params@return_dimred),
    params@extra, "ScanoramaParams"
  )
  out <- do.call(scanoramaRun, args)

  res <- scanoramaPost(
    input = input, list = ll, output = out,
    return_dimred = params@return_dimred, method = "scanorama"
  )
})

#' @rdname batchCorrect
setMethod("batchCorrect", "BBKNNParams", function(input, batch, params) {
  ll <- bbknnInput(input = input, batch = batch, reduction = params@reduction)

  args <- merge_params(
    list(input = ll$red, batch = ll$batch),
    params@extra, "BBKNNParams"
  )
  out <- do.call(bbknnRun, args)

  res <- bbknnPost(input = input, output = out, method = "bbknn")
})

