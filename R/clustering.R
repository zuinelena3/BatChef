#' Clustering after integration reduction
#'
#' Compute Leiden clustering. The resolution was chosen using the NMI metric.
#'
#' @param input input
#' @param label_true label_true
#' @param reduction reduction
#'
#' @import methods
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @export
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
      scanpy$tl$leiden(adata, resolution = res[i], key_added = "cluster", flavor = "leidenalg")
      nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = adata$obs[, label_true], labels_pred = adata$obs$cluster)

      if (nmi > MAX) {
        MAX <- nmi
        names(MAX) <- res[i]
      }
    }

    scanpy$tl$leiden(adata, resolution = as.numeric(names(MAX)), key_added = "cluster", flavor = "leidenalg")

    return(adata)
  }, input = input, label_true = label_true, reduction = reduction)
}
