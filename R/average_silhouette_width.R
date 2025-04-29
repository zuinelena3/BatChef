#' Average Silhouette Width
#'
#' @param input input
#' @param label_true label_true
#' @param reduction reduction
#' @param metric metric
#' @param sample_size sample_size
#' @param random_state random_state
#'
#' @export
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom SingleCellExperiment reducedDim colData
average_silhouette_width <- function(input, label_true, reduction, metric = "euclidean", sample_size = NULL, random_state = NULL) {

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction, metric, sample_size, random_state) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    red <- reducedDim(input, reduction)
    label <- colData(input)[, label_true]

    asw <- sklearn$metrics$silhouette_score(X = red, labels = label,
                                            metric = metric, sample_size = sample_size, random_state = random_state)
  }, input = input, label_true = label_true, reduction = reduction,  metric = metric, sample_size = sample_size, random_state = random_state)
}
