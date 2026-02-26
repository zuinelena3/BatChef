#' Normalization
#'
#' @param counts Raw counts matrix.
#'
#' @returns A normalized matrix.
#' @import Matrix
#'
normalized <- function(counts) {
  lib_sizes <- colSums(counts)
  lib_med <- median(lib_sizes)
  norm <- t(t(counts) / lib_sizes * lib_med)
}

#' Winsorization
#'
#' It's a transformation by limiting extreme values in data
#' to reduce the effect of outliers.
#'
#' @param norm_counts A normalized matrix
#' @returns Winsorized numeric vector
#'
#' @importFrom stats quantile
#'
winsorization <- function(norm_counts) {
  means <- rowMeans(norm_counts)
  means <- means[means != 0]

  q <- 0.1

  lohi <- quantile(means, c(q, 1 - q), na.rm = TRUE)

  if (diff(lohi) < 0) {
    lohi <- rev(lohi)
  }

  means[!is.na(means) & means < lohi[1]] <- lohi[1]
  means[!is.na(means) & means > lohi[2]] <- lohi[2]

  return(means)
}

#' Gamma parameters estimation of batch factors
#'
#' @param factors Batch factors
#'
#' @returns a data.frame with estimated shape and rate parameters
#'
#' @importFrom fitdistrplus fitdist
#'
params_btc_factors <- function(factors) {
  params <- lapply(factors, function(x) {
    fit <- fitdistrplus::fitdist(x, "gamma", method = "mme")
    params <- data.frame(
      shape = unname(fit$estimate["shape"]),
      rate = unname(fit$estimate["rate"])
    )
  })
  params <- do.call(rbind, params)
}

#' Kullback-Leibler (KL) divergence) between two Gamma distributions
#'
#' @param p Gamma parameters (shape and rate) of probability distribution P
#' @param q Gamma parameters (shape and rate) of probability distribution Q
#'
#' @returns Kullback-Leibler (KL) divergence) between two Gamma distributions
#'
kl_gamma <- function(p, q) {
  a1 <- p$shape
  b1 <- p$rate

  a2 <- q$shape
  b2 <- q$rate

  a2 * log(b1 / b2) - (lgamma(a1) - lgamma(a2)) +
    (a1 - a2) * digamma(a1) - (b1 - b2) * (a1 / b1)
}

#' Batch effects strength
#'
#' @param norm_counts Normalized counts matrix.
#' @param batch A string specifying the batch variable.
#'
#' @returns Median of symmetric Kullback-Leibler (KL) divergence
#' @importFrom utils combn
#'
batch_params <- function(norm_counts, batch) {
  over_means <- winsorization(norm_counts = norm_counts)

  batches <- unique(batch)

  batch_means <- lapply(batches, function(b) {
    winsorization(norm_counts = norm_counts[, batch == b])
  })

  over_means <- lapply(batch_means, function(x) over_means[names(over_means) %in% names(x)])

  batch_factors <- lapply(1:length(batches), function(x) batch_means[[x]] / over_means[[x]])

  params <- params_btc_factors(factors = batch_factors)

  combs <- combn(1:nrow(params), 2)

  sym_kl <- numeric(ncol(combs))

  for (i in 1:ncol(combs)) {
    kl_pq <- kl_gamma(p = params[combs[1, i], ], q = params[combs[2, i], ])
    kl_qp <- kl_gamma(p = params[combs[2, i], ], q = params[combs[1, i], ])
    sym_kl[i] <- kl_pq + kl_qp
  }

  params <- data.frame(batch_strength = median(sym_kl))
}

#' Library size parameters
#'
#' @param counts Raw counts matrix.
#'
#' @returns Median of sequencing depth
#' @importFrom fitdistrplus fitdist
#'
lib_size_params <- function(counts) {
  lib_sizes <- colSums(counts)

  params <- data.frame(seq_depth = median(lib_sizes))

  return(params)
}

#' Highly espressed genes probability
#'
#' @param norm_counts A normalized matrix
#'
#' @returns Highly espressed genes probability
#'
outlier_params <- function(norm_counts) {
  means <- rowMeans(norm_counts)
  lmeans <- log(means)

  med <- median(lmeans)
  mad <- mad(lmeans)

  bound <- med + 2 * mad

  outs <- which(lmeans > bound)

  out_prob <- length(outs) / nrow(norm_counts)

  params <- data.frame(out_prob)

  return(params)
}
