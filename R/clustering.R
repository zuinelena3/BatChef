#' Leiden clustering
#'
#' Leiden clustering was computed with different resolution parameters
#' (from 0.1 to 2). The clustering output with the highest Normalized Mutual
#' Information (NMI) was used to compute the other metrics.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param label_true A string specifying cell types.
#' @param reduction A string specifying the dimensional reduction.
#'
#' @export
#' @import methods
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#'
#' @return A `AnnData` object.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' clust <- clustering(input = sim, label_true = "Group", reduction = "PCA")
#'
clustering <- function(input, label_true, reduction) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, label_true, reduction) {
    scanpy <- reticulate::import("scanpy")
    sklearn <- reticulate::import("sklearn")

    adata <- anndataInput(input = input)
    scanpy$pp$neighbors(adata, use_rep = reduction)
    res <- seq(0.1, 2, by = 0.1)
    MAX <- 0

    for (i in seq_along(res)) {
      scanpy$tl$leiden(adata, resolution = res[i], key_added = "cluster",
                       flavor = "leidenalg")
      nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = adata$obs[, label_true],
                                                          labels_pred = adata$obs$cluster)

      if (nmi > MAX) {
        MAX <- nmi
        names(MAX) <- res[i]
      }
    }

    scanpy$tl$leiden(adata, resolution = as.numeric(names(MAX)),
                     key_added = "cluster", flavor = "leidenalg")

    return(adata)
  }, input = input, label_true = label_true, reduction = reduction)
}
