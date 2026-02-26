test_that("Suggested method", {
  suppressMessages({
    sim <- simulate_data(
      n_genes = 1000, batch_cells = c(130, 110),
      group_prob = c(0.3, 0.2, 0.5), n_hvgs = 1000,
      compute_pca = FALSE, output_format = "Seurat"
    )

    pred <- suggested_method(input = sim, batch = "Batch")
    expect_true(inherits(pred, "ggplot"))
  })
})
