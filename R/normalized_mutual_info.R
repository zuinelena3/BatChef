#' Normalized Mutual Information (NMI)
#'
#' Before NMI, the Leiden clustering is computed.
#'
#' @param input A `SingleCellExperiment`, `Seurat` or `AnnData` objects can be supplied.
#' @param label_true A string specifying cell types.
#' @param average_method How to compute the normalizer in the denominator.
#' @param reduction A string specifying the dimensional reduction.
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @export
#'
normalized_mutual_info <- function(input, label_true, average_method = "arithmetic", reduction) {
  adata <- clustering(input = input, label_true = label_true, reduction = reduction)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, average_method, reduction) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = input$obs[, label_true], labels_pred = input$obs$cluster,
                                                        average_method = average_method)
  }, input = adata, label_true = label_true, average_method = average_method, reduction = reduction)
}
