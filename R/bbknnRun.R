#' BBKNN method
#'
#' @param input input
#' @param batch batch
#' @param neighbors_within_batch neighbors_within_batch
#' @param n_pcs n_pcs
#' @param trim trim
#' @param computation computation
#' @param annoy_n_trees annoy_n_trees
#' @param pynndescent_n_neighbors pynndescent_n_neighbors
#' @param pynndescent_random_state pynndescent_random_state
#' @param metric metric
#' @param set_op_mix_ratio set_op_mix_ratio
#' @param local_connectivity local_connectivity
#' @param approx approx
#' @param use_annoy use_annoy
#' @param use_faiss use_faiss
#' @param scanpy_logging scanpy_logging
#'
#' @export
#' @importFrom stats cmdscale
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
bbknnRun <- function(input, batch, neighbors_within_batch = 3, n_pcs = 50, trim = NULL,
                     computation = "annoy", annoy_n_trees = 10, pynndescent_n_neighbors = 30,
                     pynndescent_random_state = 0, metric = "euclidean", set_op_mix_ratio = 1,
                     local_connectivity = 1, approx = NULL, use_annoy = NULL, use_faiss = NULL,
                     scanpy_logging = FALSE) {

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, batch, neighbors_within_batch = 3, n_pcs = 50, trim = NULL,
                                                 computation = "annoy", annoy_n_trees = 10, pynndescent_n_neighbors = 30,
                                                 pynndescent_random_state = 0, metric = "euclidean", set_op_mix_ratio = 1,
                                                 local_connectivity = 1, approx = NULL, use_annoy = NULL, use_faiss = NULL,
                                                 scanpy_logging = FALSE) {
    bbknn <- import("bbknn")

    args <- c(list(pca = input, batch_list = batch, neighbors_within_batch = as.integer(neighbors_within_batch), n_pcs = as.integer(n_pcs),
                   trim = trim, computation = computation, annoy_n_trees = as.integer(annoy_n_trees), pynndescent_n_neighbors = as.integer(pynndescent_n_neighbors),
                   pynndescent_random_state = as.integer(pynndescent_random_state), metric = metric, set_op_mix_ratio = as.integer(set_op_mix_ratio),
                   local_connectivity = as.integer(local_connectivity), approx = approx, use_annoy = use_annoy, use_faiss = use_faiss,
                   scanpy_logging = scanpy_logging))
    out <- do.call(bbknn$matrix$bbknn, args)
    mds <- cmdscale(out[[1]], k = n_pcs)
    rownames(mds) <- rownames(input)
    colnames(mds) <- paste0("bbknn_", seq(n_pcs))
    return(mds)
  }, input = input, batch = batch, neighbors_within_batch = neighbors_within_batch, n_pcs = n_pcs,
  trim = trim, computation = computation, annoy_n_trees = annoy_n_trees, pynndescent_n_neighbors = pynndescent_n_neighbors,
  pynndescent_random_state = pynndescent_random_state, metric = metric, set_op_mix_ratio = set_op_mix_ratio,
  local_connectivity = local_connectivity, approx = approx, use_annoy = use_annoy, use_faiss = use_faiss,
  scanpy_logging = scanpy_logging)
}
