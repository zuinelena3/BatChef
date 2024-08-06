test_that("limma integration method works", {
  data("data")

  limma_out <- scIntegration(data, batch = "Batch", METHOD = limmaMethod())
  expect_type(limma_out, "double")

  limma_out <- scIntegration(obj = data, batch = "Batch", METHOD = limmaMethod(), alt_out = TRUE, cell_type = "Group")
  expect_s4_class(limma_out, "AltOutput")
})

test_that("ComBat integration method works", {
  data("data")

  combat_out <- scIntegration(obj = data, assay = "logcounts", batch = "Batch", cell_type = "Group", METHOD = combatMethod())
  expect_s4_class(combat_out, "SingleCellExperiment")

  combat_out <- scIntegration(obj = data, assay = "logcounts", batch = "Batch", cell_type = "Group", METHOD = combatMethod(), alt_out = TRUE)
  expect_s4_class(combat_out, "AltOutput")
})

test_that("seuratv3 integration method works", {
  data("data")

  seuratv3_out <- scIntegration(data, batch = "Batch", cell_type = "Group", anchor = "cca", k_anchor = 5, dims = 10,
                                hvgs = rownames(data), reduction = "PCA", METHOD = seuratv3Method())
  expect_s4_class(seuratv3_out, "Seurat")

  seuratv3_out <- scIntegration(data, batch = "Batch", cell_type = "Group", anchor = "cca", k_anchor = 5, dims = 10,
                                hvgs = rownames(data), reduction = "PCA", METHOD = seuratv3Method(), alt_out = TRUE)
  expect_s4_class(seuratv3_out, "AltOutput")
})

test_that("SeuratV5 integration method works", {
  data("data")

  seuratv5_out <- scIntegration(data, anchor = "CCAIntegration", batch = "Batch", cell_type = "Group",
                                k_anchor = 5, dims = 10, hvgs = rownames(data), reduction = "PCA", METHOD = seuratv5Method())
  expect_s4_class(seuratv5_out, "Seurat")

  seuratv5_out <- scIntegration(data, anchor = "CCAIntegration", batch = "Batch", cell_type = "Group",
                                k_anchor = 5, dims = 10, hvgs = rownames(data), reduction = "PCA", METHOD = seuratv5Method(), alt_out = TRUE)
  expect_s4_class(seuratv5_out, "AltOutput")
})
