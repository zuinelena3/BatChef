#' Title
#'
#' @param obj a SingleCellExperiment object
#' @param batch a string specifying the batch
#' @param reduction a string specifying the dimension reduction name
#' @param n number of random selection cells
#' @param rep number of random selection replicates
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom transport pp wasserstein
#' @export
#'
wasserstein_metric <- function(obj, batch, reduction, n, rep) {
  ll <- lapply(unique(colData(obj)[, batch]), function(i) obj[, colData(obj)[, batch] == i])
  emb <- lapply(ll, function(x) reducedDim(x, reduction))

  wass <- lapply(1:rep, function(i) {
    random <- lapply(emb, function(x) x[sample(nrow(x), n), ])
    obj_pp <- lapply(random, function(x) pp(as.matrix(x)))
    wass_metric <- wasserstein(obj_pp[[1]], obj_pp[[2]], p = 2, prob = TRUE)
  })

  wass <- as.data.frame(do.call(rbind, wass))
  colnames(wass)[1] <- "wasserstein"
  wass$method <- reduction
  return(wass)
}
