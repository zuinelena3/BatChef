#' Suggested method prediction
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#'
#' @returns A list containing two elements:a string specifying the recommended method;
#' a ggplot object that visualizes the data points in a two-dimensional space
#' derived from the characteristics of 130 datasets.
#'
#' @export
#'
#' @import e1071
#'
#' @examples
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 1000,
#'   compute_pca = FALSE, output_format = "SingleCellExperiment"
#' )
#' pred <- suggested_method(input = sim, batch = "Batch")
#'
suggested_method <- function(input, batch) {
  sce <- linearInput(input = input, batch = batch)
  params <- extract_features(input = sce, batch = batch)

  file_path <- system.file("extdata", "svm.rda", package = "BatChef")
  load(file_path)

  svm_best <- res$svm
  pred <- as.character(predict(svm_best, params))

  msg <- paste0("Optimal method: ", pred)
  message(msg)
  prediction_plot(params = params)
}
