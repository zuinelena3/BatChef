test_that("fake_data function works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2))
  expect_s4_class(data, "SingleCellExperiment")
})
