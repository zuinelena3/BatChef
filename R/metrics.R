#' Evaluation metrics function
#'
#' @param input input
#' @param batch batch
#' @param group group
#' @param reduction reduction
#' @param rep rep
#' @param mc.cores mc.cores
#' @param average_method average_method
#' @param metric metric
#' @param sample_size sample_size
#' @param random_state random_state
#' @param meta_data meta_data
#' @param perplexity perplexity
#' @param nn_eps nn_eps
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
metrics <- function(input, batch, group, reduction, rep = 10, mc.cores = 2, average_method = "arithmetic",
                    metric = "euclidean", sample_size = NULL, random_state = NULL, meta_data = colData(input), perplexity = 30, nn_eps = 0) {
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
