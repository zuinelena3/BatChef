#' Normalized Mutual Information
#'
#' @param input input
#' @param label_true label_true
#' @param average_method average_method
#' @param reduction reduction
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
