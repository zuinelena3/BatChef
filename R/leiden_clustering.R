#' Leiden clustering
#'
#' Leiden clustering algorithm is implemented and executed using
#' the Scanpy Python package.
#'
#' The clustering algorithm can be executed by specifying either a single
#' resolution parameter or range of resolution parameters (from 0.1 to 2).
#' In the case of multiple resolutions, the clustering outcome that corresponds
#' to the highest Normalized Mutual Information (NMI) score is selected. Finding
#' the optimal clustering can be useful to compute the performance evaluation of
#' batch correction methods.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param label_true A string specifying the ground truth label.
#' @param reduction A string specifying the dimensional reduction
#' on which the clustering analysis will be performed.
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed (Default: TRUE).
#' @param resolution A numeric value specifying the resolution parameter.
#'
#' @export
#' @import methods
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @return A `AnnData` object that contains the clustering labels.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' clust <- leiden_clustering(input = sim, reduction = "PCA",
#'                            nmi_compute = FALSE, resolution = 1)
#'
leiden_clustering <- function(input, label_true = NULL, reduction,
                              nmi_compute = TRUE, resolution = NULL) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true,
                                                 reduction, nmi_compute,
                                                 resolution) {
    scanpy <- reticulate::import("scanpy")
    sklearn <- reticulate::import("sklearn")
    np <- reticulate::import("numpy")

    if (inherits(input, "Seurat")) {
      input <- Seurat::as.SingleCellExperiment(input)
    }
    if (inherits(input, "SingleCellExperiment")) {
      adata <- zellkonverter::SCE2AnnData(input, X_name = "counts")
      adata <- reticulate::r_to_py(adata, convert = TRUE)
    } else {
      adata <- input
    }

    scanpy$pp$neighbors(adata, use_rep = reduction)

    if (nmi_compute == FALSE) {
      stopifnot("Specify the resolution parameter" = !is.null(resolution))

      scanpy$tl$leiden(adata, resolution = resolution, key_added = "cluster",
                       flavor = "leidenalg")
      return(adata)
    }

    else {
      stopifnot("Specify the true label" = !is.null(label_true))

      res <- seq(0.1, 2, by = 0.1)
      max <- 0

      for (i in seq_along(res)) {
        scanpy$tl$leiden(adata, resolution = res[i], key_added = "cluster",
                         flavor = "leidenalg")

        nmi <- sklearn$metrics$normalized_mutual_info_score(
          labels_true = np$array(adata$obs[, label_true]),
          labels_pred = np$array(adata$obs$cluster))

        if (nmi > max) {
          max <- nmi
          names(max) <- res[i]
        }
      }
      scanpy$tl$leiden(adata, resolution = as.numeric(names(max)),
                       key_added = "cluster", flavor = "leidenalg")
      return(adata)
    }
  }, input = input, label_true = label_true, reduction = reduction,
  nmi_compute = nmi_compute, resolution = resolution)
}
