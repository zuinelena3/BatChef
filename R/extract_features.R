#' Extract data characteristics
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment} object can be supplied.
#' @param batch A string specifying the batch variable.
#'
#' @returns A data.frame that contains the features of input data.
#' @importFrom SingleCellExperiment colData counts
#'
extract_features <- function(input, batch) {
  stopifnot("Specify batch label!" = !missing(batch))

  ncells <- ncol(input)
  ngenes <- nrow(input)
  nbatch <- length(unique(colData(input)[, batch]))

  counts <- counts(input)
  counts <- counts[rowSums(counts) > 0, ]
  norm_counts <- normalized(counts)

  batch_params <- batch_params(
    normalized = norm_counts,
    batch = colData(input)[, batch]
  )

  mean <- mean_params(norm_counts)

  lib_size <- lib_size_params(counts)

  out <- outlier_params(norm_counts)

  groups <- groups_params(input = input, normalized = norm_counts)

  params <- data.frame(
    ncells = ncells, ngenes = ngenes, nbatch = nbatch,
    batch_params, mean, lib_size, out, groups
  )
  colnames(params) <- paste0("feature_", colnames(params))
  return(params)
}
