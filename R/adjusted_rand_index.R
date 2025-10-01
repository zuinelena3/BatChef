#' Adjusted Rand Index (ARI)
#'
#' The Adjusted Rand Index (ARI) is a metric used to measure the similarity
#' between the clustering and ground truth labels.
#'
#' After computing clustering using Leiden algorithm, the ARI metric is computed
#' using the `adjusted_rand_score` function from the Python sklearn package.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the name of dimensional reduction.
#'
#' @export
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#'
#' @return A numeric value
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' ari <- adjusted_rand_index(input = sim, label_true = "Group",
#'                            reduction = "PCA")
#'
adjusted_rand_index <- function(input, label_true, reduction) {
  adata <- clustering(input = input, label_true = label_true,
                      reduction = reduction)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    ari <- sklearn$metrics$adjusted_rand_score(labels_true = input$obs[, label_true],
                                               labels_pred = input$obs$cluster)
  }, input = adata, label_true = label_true, reduction = reduction)
}
