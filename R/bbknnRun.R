#' BBKNN method
#'
#' BBKNN is a graph-based batch correction method.
#'
#' @param input Dimensional reduction coordinates (such as Principal Component
#' Analysis) for each cell, with cells as rows.
#' @param batch Vector of batch assignments for each cell.
#' @param neighbors_within_batch Number of many top neighbours to
#' report for each batch.
#' @param n_pcs Number of dimensions.
#' @param trim Trim the neighbours of each cell to these many top connectivities.
#' @param computation A string specifying the KNN algorithm to use.
#' @param annoy_n_trees The number of trees to construct in the annoy forest.
#' @param pynndescent_n_neighbors The number of neighbours to include
#' in the approximate neighbour graph.
#' @param pynndescent_random_state The RNG seed to use when creating the graph.
#' @param metric A string specifying the distance metric to use.
#' @param set_op_mix_ratio UMAP connectivity computation parameter
#' @param local_connectivity UMAP connectivity computation parameter,
#' @param approx A logical to approximate nearest neighbor method selection.
#' @param use_annoy Flag to explicitly enable or disable Annoy for neighbor search.
#' @param use_faiss Flag to explicitly enable or disable
#' FAISS (Facebook AI Similarity Search) for neighbor search.
#' @param scanpy_logging If True, enables logging compatible with Scanpy’s
#' logging system for better integration.
#'
#' @export
#' @importFrom stats cmdscale
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @return A matrix of the coordinates of the points chosen to represent
#' the dissimilarities.
#' @examples
#' sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
#'                       group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
#' bbknn <- bbknnRun(input = SingleCellExperiment::reducedDim(sim, "PCA"),
#'                   batch = sim$Batch, n_pcs = 10)
#'
bbknnRun <- function(input, batch, neighbors_within_batch = 3, n_pcs = 50,
                     trim = NULL, computation = "annoy", annoy_n_trees = 10,
                     pynndescent_n_neighbors = 30, pynndescent_random_state = 0,
                     metric = "euclidean", set_op_mix_ratio = 1,
                     local_connectivity = 1, approx = NULL, use_annoy = NULL,
                     use_faiss = NULL, scanpy_logging = FALSE) {

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, batch,
                                                 neighbors_within_batch = 3,
                                                 n_pcs = 50, trim = NULL,
                                                 computation = "annoy",
                                                 annoy_n_trees = 10,
                                                 pynndescent_n_neighbors = 30,
                                                 pynndescent_random_state = 0,
                                                 metric = "euclidean",
                                                 set_op_mix_ratio = 1,
                                                 local_connectivity = 1,
                                                 approx = NULL, use_annoy = NULL,
                                                 use_faiss = NULL,
                                                 scanpy_logging = FALSE) {
    bbknn <- import("bbknn")

    args <- c(list(pca = input, batch_list = batch,
                   neighbors_within_batch = as.integer(neighbors_within_batch),
                   n_pcs = as.integer(n_pcs), trim = trim, computation = computation,
                   annoy_n_trees = as.integer(annoy_n_trees),
                   pynndescent_n_neighbors = as.integer(pynndescent_n_neighbors),
                   pynndescent_random_state = as.integer(pynndescent_random_state),
                   metric = metric, set_op_mix_ratio = as.integer(set_op_mix_ratio),
                   local_connectivity = as.integer(local_connectivity),
                   approx = approx, use_annoy = use_annoy, use_faiss = use_faiss,
                   scanpy_logging = scanpy_logging))
    out <- do.call(bbknn$matrix$bbknn, args)
    mds <- cmdscale(out[[1]], k = n_pcs)
    rownames(mds) <- rownames(input)
    colnames(mds) <- paste0("bbknn_", seq(n_pcs))
    return(mds)
  }, input = input, batch = batch, neighbors_within_batch = neighbors_within_batch,
  n_pcs = n_pcs, trim = trim, computation = computation,
  annoy_n_trees = annoy_n_trees, pynndescent_n_neighbors = pynndescent_n_neighbors,
  pynndescent_random_state = pynndescent_random_state, metric = metric,
  set_op_mix_ratio = set_op_mix_ratio, local_connectivity = local_connectivity,
  approx = approx, use_annoy = use_annoy, use_faiss = use_faiss,
  scanpy_logging = scanpy_logging)
}
