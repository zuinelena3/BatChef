test_that("Performance evaluation", {
  suppressMessages(capture.output({
  sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
                       group_prob = c(0.5, 0.5), n_hvgs = 500,
                       compute_pca = TRUE, output_format = "SingleCellExperiment")
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = HarmonyParams())
  metrics <- metrics(input = sim, batch = "Batch", group = "Group",
                     reduction = "harmony", rep = 5, n_iter = 2)
  expect_s3_class(metrics, "data.frame")
  expect_match(colnames(metrics), "wasserstein", all = FALSE)
  expect_match(colnames(metrics), "ari", all = FALSE)
  expect_match(colnames(metrics), "iasw", all = FALSE)
  expect_match(colnames(metrics), "ilisi", all = FALSE)
  expect_match(colnames(metrics), "nmi", all = FALSE)
  expect_match(colnames(metrics), "casw", all = FALSE)
  expect_match(colnames(metrics), "clisi", all = FALSE)
  }))
})
