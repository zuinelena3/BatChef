#' Wasserstein distance
#'
#' @param obj a SingleCellExperiment object
#' @param batch a string specifying the batch
#' @param reduction a string specifying the dimension reduction name
#' @param rep number of random selection replicates
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom transport pp wasserstein
#' @export
#'
wasserstein_metric <- function(obj, batch, reduction, rep) {
  ll <- lapply(unique(colData(obj)[, batch]), function(i) obj[, colData(obj)[, batch] == i])
  emb <- lapply(ll, function(x) reducedDim(x, reduction))

  n <- min(table(colData(obj)[, batch]))

  wass <- lapply(1:rep, function(i) {
    random <- lapply(emb, function(x) x[sample(nrow(x), n), ])
    obj_pp <- lapply(random, function(x) pp(as.matrix(x)))

    ll <- list()
    pair <- numeric(0)
    ind <- 0
    for (i in 1:length(obj_pp)) {
      for(j in 1:length(obj_pp)) {
        ind <- ind +1
        pair[ind] <- paste(i, j, sep = "_")
        ll[[ind]] <- wasserstein(obj_pp[[i]], obj_pp[[j]], p = 2, prob = TRUE)
      }
    }
    df <- data.frame(wass = do.call(rbind, ll), pair = pair)
  })

  wass <- as.data.frame(do.call(rbind, wass))
  colnames(wass)[1] <- "wasserstein"
  wass$method <- reduction
  wass <- wass[wass$wasserstein > 0, ]
  return(wass)
}
