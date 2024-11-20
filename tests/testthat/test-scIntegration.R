test_that("limma integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  limma_out <- scIntegration(data, batch = "Batch", METHOD = limmaMethod())
  expect_type(limma_out, "double")

  limma_out <- scIntegration(obj = data, batch = "Batch", METHOD = limmaMethod(), alt_out = TRUE, cell_type = "Group")
  expect_s4_class(limma_out, "AltOutput")
})

test_that("ComBat integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  combat_out <- scIntegration(obj = data, assay = "logcounts", batch = "Batch", cell_type = "Group", METHOD = combatMethod())
  expect_type(combat_out, "double")

  combat_out <- scIntegration(obj = data, assay = "logcounts", batch = "Batch", cell_type = "Group", METHOD = combatMethod(), alt_out = TRUE)
  expect_s4_class(combat_out, "AltOutput")
})

test_that("seuratv3 integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  seuratv3_out <- scIntegration(data, batch = "Batch", cell_type = "Group", anchor = "cca", k_anchor = 5, dims = 10,
                                hvgs = rownames(data), reduction = "PCA", METHOD = seuratv3Method())
  expect_s4_class(seuratv3_out, "Seurat")

  seuratv3_out <- scIntegration(data, batch = "Batch", cell_type = "Group", anchor = "cca", k_anchor = 5, dims = 10,
                                hvgs = rownames(data), reduction = "PCA", METHOD = seuratv3Method(), alt_out = TRUE)
  expect_s4_class(seuratv3_out, "AltOutput")
})

test_that("SeuratV5 integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  seuratv5_out <- scIntegration(data, anchor = "CCAIntegration", batch = "Batch", cell_type = "Group",
                                k_anchor = 5, dims = 10, hvgs = rownames(data), reduction = "PCA", METHOD = seuratv5Method())
  expect_s4_class(seuratv5_out, "Seurat")

  seuratv5_out <- scIntegration(data, anchor = "CCAIntegration", batch = "Batch", cell_type = "Group",
                                k_anchor = 5, dims = 10, hvgs = rownames(data), reduction = "PCA", METHOD = seuratv5Method(), alt_out = TRUE)
  expect_s4_class(seuratv5_out, "AltOutput")
})

test_that("fastMNN integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  mnn_out <- scIntegration(data, batch = "Batch", hvgs = rownames(data), dims = 10, METHOD = fastMNNMethod())
  expect_s4_class(mnn_out, "SingleCellExperiment")

 mnn_out <- scIntegration(data, batch = "Batch", cell_type = "Group", hvgs = rownames(data), dims = 10, METHOD = fastMNNMethod(), alt_out = TRUE)
  expect_s4_class(mnn_out, "AltOutput")
})

test_that("harmony integration method works", {
    data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  harmony_out <- scIntegration(data, batch = "Batch", METHOD = harmonyMethod())
  expect_s4_class(harmony_out, "SingleCellExperiment")

  harmony_out <- scIntegration(data, batch = "Batch", cell_type = "Group", METHOD = harmonyMethod(), alt_out = TRUE)
  expect_s4_class(harmony_out, "AltOutput")
})

test_that("scanorama integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  ll <- lapply(unique(data$Batch), function(i) data[, data$Batch == i])

  batch <- data$Batch
  cell_type <-data$Group

  assaylist <- list()
  genelist <- list()
  for(i in seq_along(ll)) {
    assaylist[[i]] <- t(as.matrix(logcounts(ll[[i]])))
    genelist[[i]] <- rownames(ll[[i]])
  }

  scanorama_out <- scIntegration(assaylist, genelist = genelist, hvgs = length(rownames(data)), dims = 10, METHOD = scanoramaMethod())
  expect_type(scanorama_out, "list")

  scanorama_out <- scIntegration(assaylist, genelist = genelist, hvgs = length(rownames(data)), cell_type = cell_type, batch = batch, dims = 10, METHOD = scanoramaMethod(), alt_out = TRUE)
  expect_s4_class(scanorama_out, "AltOutput")
})

test_that("BBKNN integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  bbknn_out <- scIntegration(data, batch = "Batch", reduction = "PCA", METHOD = bbknnMethod())
  expect_s4_class(bbknn_out, "SingleCellExperiment")

  bbknn_out <- scIntegration(data, batch = "Batch", reduction = "PCA", METHOD = bbknnMethod(), cell_type = "Group", alt_out = TRUE)
  expect_s4_class(bbknn_out, "AltOutput")
})

# test_that("scVI integration method works", {
#   data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)
#
#   scvi_out <- scIntegration(data, batch = "Batch", METHOD = scVIMethod(), assay = "counts", alt_out = TRUE, cell_type = "Group")
#   expect_s4_class(scvi_out, "AltOutput")
# })

test_that("scMerge integration method works", {
  data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

  scmerge_out <- scIntegration(data, batch = "Batch", METHOD = scMergeMethod(), genelist = rownames(data), alt_out = TRUE, cell_type = "Group")
  expect_s4_class(scmerge_out, "AltOutput")
})
