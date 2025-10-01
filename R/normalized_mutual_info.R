#' Normalized Mutual Information (NMI)
#'
#' The Normalized Mutual Information (NMI) metric is used to compare
#' the overlap between the true cell-type and clustering labels computed after
#' batch correction.
#'
#' After computing clustering using Leiden algorithm, the ARI metric is computed
#' using the `normalized_mutual_info_score` function from the Python sklearn package.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param label_true A string specifying cell types.
#' @param average_method How to compute the normalizer in the denominator.
#' @param reduction A string specifying the dimensional reduction.
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @return A numeric value.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' nmi <- normalized_mutual_info(input = sim, label_true = "Group",
#'                               reduction = "PCA")
#'
normalized_mutual_info <- function(input, label_true,
                                   average_method = "arithmetic", reduction) {
  adata <- clustering(input = input, label_true = label_true, reduction = reduction)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true,
                                                 average_method, reduction) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    nmi <- sklearn$metrics$normalized_mutual_info_score(
      labels_true = input$obs[, label_true],
      labels_pred = input$obs$cluster,
      average_method = average_method)
  }, input = adata, label_true = label_true, average_method = average_method,
  reduction = reduction)
}
