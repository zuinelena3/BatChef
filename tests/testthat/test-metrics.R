test_that("metrics function works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  out_pre <- clustering(data, cell_type = "Group", reduction = "pca")
  metrics_pre <- metrics(anndata = out_pre, cell_type = "Group", batch = "Batch", reduction = "pca")
  expect_type(metrics_pre, "list")
})
