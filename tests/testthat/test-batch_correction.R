test_that("Scanorama method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(
      n_genes = 500, batch_cells = c(50, 50),
      group_prob = c(0.5, 0.5), n_hvgs = 500,
      compute_pca = FALSE,
      output_format = "SingleCellExperiment"
    )
    sim <- batchCorrect(
      input = sim, batch = "Batch",
      params = ScanoramaParams(assay_type = "logcounts")
    )
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("BBKNN method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(
      n_genes = 500, batch_cells = c(100, 50),
      group_prob = c(0.5, 0.5), n_hvgs = 500,
      compute_pca = TRUE,
      output_format = "SingleCellExperiment"
    )
    sim <- batchCorrect(
      input = sim, batch = "Batch",
      params = BBKNNParams(reduction = "PCA")
    )
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})
