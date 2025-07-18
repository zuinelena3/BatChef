#' Adjusted Rand Index (ARI)
#'
#' Before ARI, the Leiden clustering is computed.
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#'
#' @export
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#'
adjusted_rand_index <- function(input, label_true, reduction) {
  adata <- clustering(input = input, label_true = label_true, reduction = reduction)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    ari <- sklearn$metrics$adjusted_rand_score(labels_true = input$obs[, label_true],
                                               labels_pred = input$obs$cluster)
  }, input = adata, label_true = label_true, reduction = reduction)
}
