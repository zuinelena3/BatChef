test_that("multiplication works", {
  data <- fake_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2))
  expect_s4_class(data, "SingleCellExperiment")
})
