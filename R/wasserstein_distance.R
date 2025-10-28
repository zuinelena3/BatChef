#' Wasserstein distance
#'
#' The Wasserstein distance measures the minimal transport needed to shift one
#' distribution to match another.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param batch A string specifying batch variable.
#' @param reduction A string specifying the dimensional reduction.
#' @param rep Number of times the Wasserstein distance is calculated.
#' @param mc_cores The number of cores to use.
#'
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom parallel mclapply
#' @importFrom transport pp wasserstein
#'
#' @return A numeric value.
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(130, 110),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                     compute_pca = TRUE, output_format = "SingleCellExperiment")
#' wass <- wasserstein_distance(input = sim, batch = "Batch",
#'                              reduction = "PCA", rep = 1, mc_cores = 1)
#'
wasserstein_distance <- function(input, batch, reduction, rep, mc_cores) {
  stopifnot(batch %in% colnames(colData(input)))

  stopifnot("Error: N. of batches has to be greater than 1" =
              length(unique(colData(input)[, batch])) > 1)

  ll <- lapply(unique(colData(input)[, batch]),
               function(i) input[, colData(input)[, batch] == i])
  embed <- lapply(ll, function(x) reducedDim(x, reduction))

  n <- min(sapply(embed, function(x) length(rownames(x))))

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
  }, mc.cores = mc_cores)
  wass <- lapply(wass, function(x) {
    df <- data.frame(wasserstein = do.call(rbind, x), pair = names(x),
                     method = reduction)
    rownames(df) <- NULL
    df
  })
  wass <- do.call(rbind, wass)
  wass <- wass[wass$wasserstein > 0, ]
}
