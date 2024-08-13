#' Clustering after integration reduction
#'
#' Compute Leiden clustering. The resolution was chosen using the NMI metric.
#'
#' @param obj a SingleCellExperiment object
#' @param cell_type a string specifying cell-type labels
#' @param reduction a string specifying the dimension reduction name
#'
#' @import methods
#'
#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom SingleCellExperiment reducedDimNames reducedDimNames<-
#' @export
#'
clustering <- function(obj, cell_type, reduction) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, cell_type, reduction, features) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")
    anndata <- import("anndata")

    reducedDimNames(obj) <- paste0("X_", tolower(reducedDimNames(obj)))
    adata <- SCE2AnnData(obj)
    adata$layers["counts"] <- adata$X

    scanpy$pp$neighbors(adata, use_rep = paste0("X_", tolower(reduction)))
    res <- seq(0.1, 2, by = 0.1)
    MAX <- 0

    for (i in seq_along(res)) {
      scanpy$tl$leiden(adata, resolution = res[i], key_added = "cluster", flavor = "leidenalg")
      nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = adata$obs[, cell_type], labels_pred = adata$obs$cluster)

      if (nmi > MAX) {
        MAX <- nmi
        names(MAX) <- res[i]
      }
    }

    scanpy$tl$leiden(adata, resolution = as.numeric(names(MAX)), key_added = "cluster", flavor = "leidenalg")

    return(adata)
  }, obj = obj, cell_type = cell_type, reduction = reduction)
}
