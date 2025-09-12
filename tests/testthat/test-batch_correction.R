test_that("Limma method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- as.Seurat(sim)
  sim <- batchCorrect(input = sim, batch = "Batch", params = LimmaParams())
  expect_s4_class(sim, "Seurat")
})

test_that("ComBat method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- as.Seurat(sim)
  sim <- batchCorrect(input = sim, batch = "Batch", params = CombatParams())
  expect_s4_class(sim, "Seurat")
})

test_that("SeuratV3 method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = SeuratV3Params(features = rownames(sim),
                                              pca_name = "PCA"))
  expect_s4_class(sim, "SingleCellExperiment")
})

test_that("SeuratV5 method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = SeuratV5Params(pca_name = "PCA",
                                              features = rownames(sim)))
  expect_s4_class(sim, "SingleCellExperiment")
})

test_that("fastMNN method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- as.Seurat(sim)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = FastMNNParams())
  expect_s4_class(sim, "Seurat")
})

test_that("Harmony method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- as.Seurat(sim)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = HarmonyParams())
  expect_s4_class(sim, "Seurat")
})

test_that("Scanorama method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = ScanoramaParams(return_dimred = TRUE,
                                               assay_type = "logcounts"))
  expect_s4_class(sim, "SingleCellExperiment")
})

test_that("scMerge2 method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- as.Seurat(sim)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = ScMerge2Params())
  expect_s4_class(sim, "Seurat")
})

test_that("scVI method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = ScVIParams(max_epochs = 3))
  expect_s4_class(sim, "SingleCellExperiment")
})

test_that("BBKNN method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = BBKNNParams(reduction = "PCA"))
  expect_s4_class(sim, "SingleCellExperiment")
})


test_that("LIGER method works", {
  sim <- simulated_data(nGenes = 1000, batchCells = c(150, 200),
                        group.prob = c(0.5, 0.5), n_hvgs = 1000, ncomp = 10)
  sim <- batchCorrect(input = sim, batch = "Batch",
                      params = LigerParams(features = rownames(sim)))
  expect_s4_class(sim, "SingleCellExperiment")
})
