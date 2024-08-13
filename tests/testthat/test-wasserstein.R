test_that("wasserstein function works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)
  wass <- wasserstein_metric(data, batch = "Batch", reduction = "PCA", n = 10, rep = 2)
  expect_type(wass, "list")
})
