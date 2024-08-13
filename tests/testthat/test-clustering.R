test_that("clustering function works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  clust_before <- clustering(data, cell_type = "Group", reduction = "pca")
  expect_type(clust_before, "environment")
})
