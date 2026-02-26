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

  input <- input[rowSums(counts(input)) > 0, colSums(counts(input)) > 0]
  n_cells <- as.numeric(ncol(input))

  batch <- colData(input)[, batch]
  n_batch <- length(unique(batch))

  counts <- counts(input)
  norm_counts <- normalized(counts)

  batch_params <- batch_params(norm_counts = norm_counts, batch = batch)
  lib_size <- lib_size_params(counts)
  out_prob <- outlier_params(norm_counts = norm_counts)

  params <- data.frame(
    n_cells = n_cells, n_batch = n_batch, lib_size,
    batch_params, out_prob
  )
  return(params)
}
