test_that("Limma method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = FALSE, output_format = "Seurat")
    sim <- batchCorrect(input = sim, batch = "Batch", params = LimmaParams())
    expect_s4_class(sim, "Seurat")
  }))
})

test_that("ComBat method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = FALSE, output_format = "AnnData")
    sim <- batchCorrect(input = sim, batch = "Batch", params = CombatParams())
    expect_true(inherits(sim, "AnnDataR6"))
  }))
})

test_that("SeuratV3 method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 1000, batch_cells = c(250, 200),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE, output_format = "AnnData")
    sim <- batchCorrect(input = sim, batch = "Batch",
                        params = SeuratV3Params(features = sim$var_names,
                                                pca_name = "X_pca"))
    expect_true(inherits(sim, "list"))
  }))
})

test_that("SeuratV5 method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(250, 200),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE, pca_ncomp = 10,
                         output_format = "SingleCellExperiment")
    sim <- batchCorrect(input = sim, batch = "Batch",
                        params = SeuratV5Params(pca_name = "PCA",
                                                features = rownames(sim)[rowData(sim)$hvgs]))
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("fastMNN method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(150, 100),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE, pca_ncomp = 10,
                         output_format = "Seurat")
    sim <- batchCorrect(input = sim, batch = "Batch", params = FastMNNParams())
    expect_s4_class(sim, "Seurat")
  }))
})

test_that("Harmony method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(150, 100),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE, pca_ncomp = 10,
                         output_format = "AnnData")

    sim <- batchCorrect(input = sim, batch = "Batch", params = HarmonyParams())
    expect_true(inherits(sim, "AnnDataR6"))
  }))
})

test_that("Scanorama method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(50, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = FALSE,
                         output_format = "SingleCellExperiment")
    sim <- batchCorrect(input = sim, batch = "Batch",
                          params = ScanoramaParams(assay_type = "logcounts"))
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("scMerge2 method works", {
  R.devices::suppressGraphics(suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 1000, batch_cells = c(200, 200),
                         group_prob = c(0.5, 0.5), n_hvgs = 1000,
                         compute_pca = TRUE, pca_ncomp = 10,
                         output_format = "Seurat")

    sim <- batchCorrect(input = sim, batch = "Batch", params = ScMerge2Params())
    expect_s4_class(sim, "Seurat")
  })))
})

test_that("scVI method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(100, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         output_format = "SingleCellExperiment")
    sim <- batchCorrect(input = sim, batch = "Batch",
                          params = ScVIParams(n_latent = 10, max_epochs = 1))
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("BBKNN method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(100, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE,
                         output_format = "SingleCellExperiment")
    sim <- batchCorrect(input = sim, batch = "Batch",
                          params = BBKNNParams(reduction = "PCA"))
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("LIGER method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(n_genes = 500, batch_cells = c(100, 50),
                         group_prob = c(0.5, 0.5), n_hvgs = 500,
                         compute_pca = TRUE, output_format = "AnnData")
    sim <- batchCorrect(input = sim, batch = "Batch",
                          params = LigerParams(features = sim$var_names,
                                               verbose = FALSE))
    expect_true(inherits(sim, "AnnDataR6"))
  }))
})
