#' Normalized Mutual Information (NMI)
#'
#' The Normalized Mutual Information (NMI) metric is used to compare
#' the overlap between the true cell-type and clustering labels computed after
#' batch correction.
#'
#' After computing clustering using Leiden algorithm, the NMI metric is computed
#' using the `normalized_mutual_info_score` function from the Python sklearn package.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param label_true A string specifying the ground truth labels.
#' @param reduction A string specifying the dimensional reduction
#' on which the clustering analysis will be performed.
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed (Default: FALSE).
#' @param resolution A numeric value specifying the resolution parameter.
#' @param average_method How to compute the normalizer in the denominator.
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @return A numeric value.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' nmi <- normalized_mutual_info(input = sim, label_true = "Group", reduction = "PCA",
#'                            nmi_compute = FALSE, resolution = 0.5)
#'
normalized_mutual_info <- function(input, label_true, reduction, nmi_compute = FALSE,
                                   resolution, average_method = "arithmetic") {

  adata <- leiden_clustering(input = input, label_true = label_true,
                             reduction = reduction, nmi_compute = nmi_compute,
                             resolution = resolution)

  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true,
                                                 average_method, reduction) {
    scanpy <- reticulate::import("scanpy")
    sklearn <- reticulate::import("sklearn")
    np <- reticulate::import("numpy")

    nmi <- sklearn$metrics$normalized_mutual_info_score(
      labels_true = np$array(input$obs[, label_true]),
      labels_pred = np$array(input$obs$cluster),
      average_method = average_method)

  }, input = adata, label_true = label_true, average_method = average_method,
  reduction = reduction)
}
