#' Common quantitative metrics
#'
#' @param anndata anndata object, computed by clustering function
#' @param cell_type string specifying cell-type labels
#' @param batch string specifying batch
#' @param reduction string specifying the integration reduction
#' @param only_ARI compute only ARI metrics
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom lisi compute_lisi
#' @importFrom stats median
#' @export
#'
metrics <- function(anndata, cell_type, batch, reduction, only_ARI = FALSE) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(anndata, cell_type, batch, reduction) {
    scanpy <- import("scanpy")
    sklearn <- import("sklearn")

    if(only_ARI == FALSE){
      nmi <- sklearn$metrics$normalized_mutual_info_score(labels_true = anndata$obs[, cell_type], labels_pred = anndata$obs$cluster)
      ari <- sklearn$metrics$adjusted_rand_score(labels_true = anndata$obs[, cell_type], labels_pred = anndata$obs$cluster)
      casw <- sklearn$metrics$silhouette_score(anndata$obsm[[paste0("X_", reduction)]], anndata$obs[, cell_type])
      clisi <- median(compute_lisi(anndata$obsm[[paste0("X_", reduction)]], data.frame(ct = anndata$obs[, cell_type]), "ct")$ct)
      homogeneity <- sklearn$metrics$homogeneity_score(anndata$obs[, cell_type], anndata$obs$cluster)

      iasw <- sklearn$metrics$silhouette_score(anndata$obsm[[paste0("X_", reduction)]], anndata$obs[, batch])
      ilisi <- median(compute_lisi(anndata$obsm[[paste0("X_", reduction)]], data.frame(batch = anndata$obs[, batch]), "batch")$batch)

      df <- data.frame(nmi, ari, casw, clisi, homogeneity, iasw, ilisi)
      return(df)
    }
    else ari <- sklearn$metrics$adjusted_rand_score(labels_true = anndata$obs[, cell_type], labels_pred = anndata$obs$cluster)

  }, anndata = anndata, cell_type = cell_type, batch = batch, reduction = reduction)
}
