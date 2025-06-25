#' Performance evaluation metrics.
#'
#' Wasserstein distance, Local Inverse Simpson's Index, Average Silhouette Width,
#' Adjusted Rand Index, and Normalized Mutual Information were computed
#'
#' @param input A `SingleCellExperiment` object can be supplied.
#' @param batch A string specifying batch variable.
#' @param group A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#' @param rep Number of times the Wasserstein distance is calculated.
#' @param mc.cores The number of cores to use.
#' @param average_method How to compute the normalizer in the denominator.
#' @param metric The metric to use when calculating distance between instances.
#' @param sample_size The size of the sample to use when computing the Silhouette Coefficient on a random subset of the data.
#' @param random_state Determines random number generation for selecting a subset of samples.
#' @param meta_data A data frame with one row per cell.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with `RANN:nn2()`.
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
metrics <- function(input, batch, group, reduction, rep = 10, mc.cores = 2,
                    average_method = "arithmetic", metric = "euclidean",
                    sample_size = NULL, random_state = NULL, meta_data = colData(input),
                    perplexity = 30, nn_eps = 0) {
  adata <- clustering(input = input, label_true = group, reduction = reduction)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction, average_method) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = input$obs[, label_true], labels_pred = input$obs$cluster,
                                                        average_method = average_method)
    ari <- sklearn$metrics$adjusted_rand_score(labels_true = input$obs[, label_true], labels_pred = input$obs$cluster)
    return(data.frame(nmi, ari))
  }, input = adata, label_true = group, reduction = reduction, average_method = average_method)

  casw <- average_silhouette_width(input = input, label_true = group, reduction = reduction, metric = metric, sample_size = sample_size,
                                   random_state = random_state)
  iasw <- average_silhouette_width(input = input, label_true = batch, reduction = reduction, metric = metric, sample_size = sample_size,
                                   random_state = random_state)
  clisi <- local_inverse_simpson_index(input = input, label_true = group, reduction = reduction, meta_data = meta_data, perplexity = perplexity,
                                       nn_eps = nn_eps)
  ilisi <- local_inverse_simpson_index(input = input, label_true = batch, reduction = reduction, meta_data = meta_data, perplexity = perplexity,
                                       nn_eps = nn_eps)
  wass <- wasserstein_distance(input = input, batch = batch, reduction = reduction, rep = rep, mc.cores = mc.cores)
  wass <- mean(wass$wasserstein)

  return(data.frame(method = reduction, wasserstein = wass, iasw, ilisi, ari = out$ari, nmi = out$nmi, casw, clisi))
}
