#' Average Silhouette Width (ASW)
#'
#' The average silhouette width for each cell is defined as the difference
#' between the average of  within-cluster distances to all cells and the average
#' between-cluster distances of that cell to the closest cluster divided
#' by their maximum.
#'
#' The ASW metric is computed using the `silhouette_score` function from the
#' Python sklearn package.
#'
#' @param input A \linkS4class{SingleCellExperiment} object.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#' @param metric The metric to use when calculating distance between instances.
#' @param sample_size The size of the sample to use when computing
#' the Silhouette Coefficient on a random subset of the data.
#' @param random_state Determines random number generation for selecting
#' a subset of samples.
#'
#' @export
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom SingleCellExperiment reducedDim colData
#'
#' @return A numeric value
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 50),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' asw <- average_silhouette_width(input = sim, label_true = "Group",
#'                                 reduction = "PCA")
average_silhouette_width <- function(input, label_true, reduction,
                                     metric = "euclidean", sample_size = NULL,
                                     random_state = NULL) {

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction,
                                                 metric, sample_size,
                                                 random_state) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    red <- reducedDim(input, reduction)
    label <- colData(input)[, label_true]

    asw <- sklearn$metrics$silhouette_score(X = red, labels = label,
                                            metric = metric,
                                            sample_size = sample_size,
                                            random_state = random_state)
  }, input = input, label_true = label_true, reduction = reduction,
  metric = metric, sample_size = sample_size, random_state = random_state)
}
