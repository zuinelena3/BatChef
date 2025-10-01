# NOTE: The code is copied from [https://github.com/immunogenomics/LISI/blob/master/R/utils.R].
# Proper credit is given to the original author.

#' Use this function to compute LISI scores of one or more labels.
#'
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell.
#' @param label_colnames Which variables to compute LISI for.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with \code{RANN:nn2()}.
#' Default of 0.0 implies exact nearest neighbor search.
#'
#' @export
#' @importFrom RANN nn2
#'
#' @return A data frame of LISI values. Each row is a cell and each
#' column is a different label variable.
#'
#' @examples
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' red <- SingleCellExperiment::reducedDim(sim, "PCA")
#' lisi <- compute_lisi(X = red, meta_data = SingleCellExperiment::colData(sim),
#'                     label_colnames = "Group")
#'
compute_lisi <- function(X, meta_data, label_colnames, perplexity = 30,
                         nn_eps = 0) {
  N <- nrow(meta_data)
  dknn <- nn2(X, k = perplexity * 3, eps = nn_eps)
  lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
  lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message('Cannot compute LISI on missing values')
      return(rep(NA, N))
    } else {
      ## don't count yourself in your neighborhood
      dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
      dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
      labels <- as.integer(factor(labels)) - 1
      n_batches <- length(unique(labels))
      simpson <- compute_simpson_index(
        t(dknn$nn.dists),
        t(dknn$nn.idx) - 1,
        labels,
        n_batches,
        perplexity
      )
      return(1 / simpson)
    }
  }))
  lisi_df <- as.data.frame(lisi_df)
  colnames(lisi_df) <- label_colnames
  row.names(lisi_df) <- row.names(meta_data)
  return(lisi_df)
}
