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

#' Batch params
#'
#' @param normalized Normalized counts matrix.
#' @param batch A string specifying the batch variable.
#'
#' @returns Batch factor location.
#' @importFrom fitdistrplus fitdist
#'
batch_params <- function(normalized, batch) {
  gene_means_overall <- rowMeans(normalized)
  gene_means_overall <- gene_means_overall[gene_means_overall > 0]

  levels <- unique(batch)
  batch_facLoc <- numeric(length(levels))
  batch_facScale <- numeric(length(levels))

  batch_facloc <- lapply(seq_along(levels), function(i) {
    batch_cells <- which(batch == levels[i])
    gene_means_batch <- rowMeans(normalized[, batch_cells])

    batch_factors <- gene_means_batch / gene_means_overall
    # If batch_factors is zero
    min_val <- 1e-6
    batch_factors <- pmax(batch_factors, min_val)

    fit <- fitdist(batch_factors, "lnorm")
    facloc <- unname(fit$estimate["meanlog"])
  })
  batch_facloc <- do.call(rbind, batch_facloc)
  batch_facloc <- median(batch_facloc[, 1])
  return(data.frame(batch_facloc = batch_facloc))
}

#' Mean parameters
#'
#' @param normalized Normalized counts matrix.

#'
#' @returns Mean rate parameter.
#' @importFrom fitdistrplus fitdist
#' @importFrom datawizard winsorize
#'
mean_params <- function(normalized) {
  means <- rowMeans(normalized)
  means <- means[means != 0]
  means <- winsorize(means, q = 0.1)

  fit <- fitdist(data = means, distr = "gamma",
                 method = "mge", gof = "CvM")
  if (fit$convergence > 0) {
    warning(
      "Fitting means using the Goodness of Fit method failed, ",
      "using the Method of Moments instead"
    )
    fit <- fitdist(means, "gamma", method = "mme")
    params <- data.frame(mean_rate = unname(fit$estimate["rate"]))
    return(params)
  }

  else {
    params <- data.frame(mean_rate = unname(fit$estimate["rate"]))
    return(params)
  }
}

#' Library size parameters
#'
#' @param counts Raw counts matrix.
#'
#' @returns Library size parameters.
#' @importFrom fitdistrplus fitdist
#' @importFrom stats shapiro.test
#'
lib_size_params <- function(counts) {
  lib_sizes <- colSums(counts)

  if (length(lib_sizes) > 5000) {
    lib_sizes_sampled <- sample(lib_sizes, 5000, replace = FALSE)
  } else {
    lib_sizes_sampled <- lib_sizes
  }

  norm_test <- shapiro.test(lib_sizes_sampled)
  p_val <- norm_test$p.value

  if (p_val >= 0.05) {
    fit <- fitdist(lib_sizes, "norm")
    lib_loc <- unname(fit$estimate["mean"])
    lib_scale <- unname(fit$estimate["sd"])
  } else {
    fit <- fitdist(lib_sizes, "lnorm")
    lib_loc <- unname(fit$estimate["meanlog"])
    lib_scale <- unname(fit$estimate["sdlog"])
  }

  params <- data.frame(seq_depth = median(lib_sizes), lib_loc, lib_scale)
  return(params)
}

#' Expression outlier parameters
#'
#' @param normalized Normalized counts matrix.
#'
#' @returns Probability that genes will be selected to be expression outliers.
#'
outlier_params <- function(normalized) {
  means <- rowMeans(normalized)
  lmeans <- log(means)

  med <- median(lmeans)
  mad <- mad(lmeans)

  bound <- med + 2 * mad
  outs <- which(lmeans > bound)

  prob <- length(outs) / nrow(normalized)

  params <- data.frame(prob)
  return(params)
}

#' Groups parameters
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param normalized A normalized counts matrix.
#'
#' @returns The number of clusters.
#'
#' @importFrom scrapper modelGeneVariances chooseHighlyVariableGenes runPca
#'
groups_params <- function(input, normalized) {
  genes <- rownames(normalized)
  gene_var <- modelGeneVariances(normalized, num.threads = 1)
  top_hvgs <- chooseHighlyVariableGenes(gene_var$statistics$residuals,
                                        top = 1000)
  top_hvgs <- rownames(gene_var$statistics[top_hvgs, ])
  normalized <- normalized[top_hvgs, ]

  pca <- runPca(normalized, num.threads = 1, number = 10)
  pca <- t(pca[[1]])
  colnames(pca) <- paste0("PC_", 1:10)
  rownames(pca) <- colnames(input)
  reducedDim(input, "PCA") <- pca

  clust <- leiden_clustering(input = input, reduction = "PCA",
                             nmi_compute = FALSE, resolution = 1, store = FALSE)
  nclust <- data.frame(ngroup = length(unique(clust)))
}
