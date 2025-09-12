test_that("Performance evaluation", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = HarmonyParams())
  metrics <- metrics(input = sim, batch = "Batch", group = "Group",
                     reduction = "harmony", rep = 2)
  expect_s3_class(metrics, "data.frame")
  expect_match(colnames(metrics), "wasserstein", all = FALSE)
  expect_match(colnames(metrics), "ari", all = FALSE)
  expect_match(colnames(metrics), "iasw", all = FALSE)
  expect_match(colnames(metrics), "ilisi", all = FALSE)
  expect_match(colnames(metrics), "nmi", all = FALSE)
  expect_match(colnames(metrics), "casw", all = FALSE)
  expect_match(colnames(metrics), "clisi", all = FALSE)
})
