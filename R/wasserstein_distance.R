#' Wasserstein distance
#'
#' @param input A `SingleCellExperiment` object.
#' @param batch A string specifying batch variable.
#' @param reduction A string specifying the dimensional reduction.
#' @param rep Number of times the Wasserstein distance is calculated.
#' @param mc.cores The number of cores to use.
#'
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom parallel mclapply
#' @importFrom transport pp wasserstein
#'
wasserstein_distance <- function(input, batch, reduction, rep, mc.cores) {
  stopifnot(batch %in% colnames(colData(input)))

  stopifnot("Error: N. of batches has to be greater than 1" = length(unique(colData(input)[, batch])) > 1)

  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  embed <- lapply(ll, function(x) reducedDim(x, reduction))

  n <- min(sapply(embed, function(x) length(rownames(x))))
  stopifnot("Error: N. of cells for a batch too small" = n >= 100)

  wass <- mclapply(1:rep, function(i) {
    random <- lapply(embed, function(x) x[sample(nrow(x), n), ])
    pp <- lapply(random, function(x) pp(as.matrix(x)))

    ll <- list()
    pair <- numeric(0)
    ind <- 0
    for (i in 1:length(pp)) {
      for(j in 1:length(pp)) {
        ind <- ind +1
        pair[ind] <- paste(i, j, sep = "_")
        ll[[ind]] <- wasserstein(pp[[i]], pp[[j]], p = 2, prob = TRUE)
      }
    }
    names(ll) <- pair
    return(ll)
  }, mc.cores = mc.cores)
  wass <- lapply(wass, function(x) {
    df <- data.frame(wasserstein = do.call(rbind, x), pair = names(x), method = reduction)
    rownames(df) <- NULL
    df
  })
  wass <- do.call(rbind, wass)
  wass <- wass[wass$wasserstein > 0, ]
}
