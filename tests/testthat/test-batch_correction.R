# test_that("Limma method works", {
#   suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 1000, batch_cells = c(150, 50),
#       group_prob = c(0.5, 0.5), n_hvgs = 500,
#       compute_pca = FALSE, output_format = "Seurat"
#     )
#     sim <- batchCorrect(input = sim, batch = "Batch",
#                         params = LimmaParams(assay_type = "data"))
#     expect_s4_class(sim, "Seurat")
#   }))
# })
#
# test_that("ComBat method works", {
#   suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 1000, batch_cells = c(150, 50),
#       group_prob = c(0.5, 0.5), n_hvgs = 500,
#       compute_pca = FALSE, output_format = "Seurat"
#     )
#     sim <- batchCorrect(input = sim, batch = "Batch",
#                         params = CombatParams(assay_type = "data"))
#     expect_s4_class(sim, "Seurat")
#   }))
# })

test_that("SeuratV3 method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(
      n_genes = 1000, batch_cells = c(250, 200),
      group_prob = c(0.5, 0.5), n_hvgs = 500,
      compute_pca = TRUE, output_format = "SingleCellExperiment"
    )
    sim <- batchCorrect(
      input = sim, batch = "Batch",
      params = SeuratV3Params(
        features = rownames(sim)[rowData(sim)$hvg],
        pca_name = "PCA"
      )
    )
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

test_that("SeuratV5 method works", {
  suppressMessages(capture.output({
    sim <- simulate_data(
      n_genes = 500, batch_cells = c(250, 200),
      group_prob = c(0.5, 0.5), n_hvgs = 500,
      compute_pca = TRUE, pca_ncomp = 10,
      output_format = "SingleCellExperiment"
    )
    sim <- batchCorrect(
      input = sim, batch = "Batch",
      params = SeuratV5Params(
        pca_name = "PCA",
        features = rownames(sim)[rowData(sim)$hvg]
      )
    )
    expect_s4_class(sim, "SingleCellExperiment")
  }))
})

# test_that("fastMNN method works", {
#   suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 500, batch_cells = c(150, 100),
#       group_prob = c(0.5, 0.5), n_hvgs = 500,
#       compute_pca = TRUE, pca_ncomp = 10,
#       output_format = "SingleCellExperiment"
#     )
#     sim <- batchCorrect(input = sim, batch = "Batch", params = FastMNNParams())
#     expect_s4_class(sim, "SingleCellExperiment")
#   }))
# })
#
# test_that("Harmony method works", {
#   suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 500, batch_cells = c(150, 100),
#       group_prob = c(0.5, 0.5), n_hvgs = 500,
#       compute_pca = TRUE, pca_ncomp = 10,
#       output_format = "Seurat"
#     )
#
#     sim <- batchCorrect(input = sim, batch = "Batch", params = HarmonyParams())
#     expect_s4_class(sim, "Seurat")
#   }))
# })
#
# test_that("scMerge2 method works", {
#   R.devices::suppressGraphics(suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 1000, batch_cells = c(200, 200),
#       group_prob = c(0.5, 0.5), n_hvgs = 1000,
#       compute_pca = TRUE, pca_ncomp = 10,
#       output_format = "Seurat"
#     )
#
#     sim <- batchCorrect(input = sim, batch = "Batch",
#                         params = ScMerge2Params(assay_type = "data", cosineNorm = FALSE))
#     expect_s4_class(sim, "Seurat")
#   })))
# })
#
# test_that("LIGER method works", {
#   suppressMessages(capture.output({
#     sim <- simulate_data(
#       n_genes = 500, batch_cells = c(100, 50),
#       group_prob = c(0.5, 0.5), n_hvgs = 500,
#       compute_pca = TRUE, output_format = "SingleCellExperiment"
#     )
#     sim <- batchCorrect(
#       input = sim, batch = "Batch",
#       params = LigerParams(
#         features = rownames(sim),
#         verbose = FALSE
#       )
#     )
#     expect_s4_class(sim, "SingleCellExperiment")
#   }))
# })
