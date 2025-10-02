#' scVI method
#'
#' scVI is a deep learning-based method.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param batch A string specifying the batch for each cell.
#' @param assay_type A string specifying the assay.
#' @param layer A string specifying the key in adata.layers for raw count data.
#' @param labels_key A string specifying the key in adata.obs for
#' label information.
#' @param size_factor_key A string specifying th key in adata.obs for
#' size factor information.
#' @param categorical_covariate_keys A string specifying the keys in adata.obs
#' that correspond to categorical data.
#' @param continuous_covariate_keys A string specifying the keys in adata.obs
#' that correspond to continuous data.
#' @param n_hidden Number of nodes per hidden layer.
#' @param n_latent Dimensionality of the latent space.
#' @param n_layers Number of hidden layers used for encoder and decoder NNs.
#' @param dropout_rate Dropout rate for neural networks.
#' @param dispersion Dispersion parameter. One of the following: 'gene' -
#' dispersion parameter of NB is constant per gene across cells. 'gene-batch' -
#' dispersion can differ between different batches. 'gene-label' - dispersion
#' can differ between different labels. 'gene-cell' - dispersion can differ
#' for every gene in every cell
#' @param gene_likelihood Gene likelihood. 'nb' - Negative binomial distribution.
#' 'zinb' (default) - Zero-inflated negative binomial distribution.
#' 'poisson' - Poisson distribution. 'normal' - EXPERIMENTAL Normal distribution
#' @param latent_distribution Distribution of latent space. Values:
#' 'normal' (default) - Normal distribution
#' 'ln' - Logistic normal distribution (Normal(0, I) transformed by softmax)
#' @param max_epochs The maximum number of epochs to train the model
#' @param accelerator Supports passing different accelerator types
#' (“cpu”, “gpu”, “tpu”, “ipu”, “hpu”, “mps, “auto”) as well as custom
#' accelerator instances.
#' @param devices The devices to use.
#' @param train_size Size of training set in the range [0.0, 1.0]. Default is NULL.
#' @param validation_size Size of the test set
#' @param shuffle_set_split Boolean (default: TRUE).
#' Whether to shuffle indices before splitting.
#' @param load_sparse_tensor Boolean value (default: FALSE). If TRUE, loads data
#' with sparse CSR or CSC layout as a Tensor with the same layout.
#' @param batch_size Minibatch size to use during training.
#' @param early_stopping Boolean value (default: FALSE). Perform early stopping.
#' @param datasplitter_kwargs Additional keyword arguments passed into DataSplitter.
#' @param plan_kwargs Additional keyword arguments passed into TrainingPlan.
#' @param datamodule A LightningDataModule instance to use for training
#' in place of the default DataSplitter.
#' @param indices Indices of observations in adata to use (default: NULL)
#' @param give_mean Boolean value (default: TRUE). If TRUE, returns
#' the mean of the latent distribution. If FALSE, returns an estimate of
#' the mean using mc_samples Monte Carlo samples.
#' @param mc_samples Number of Monte Carlo samples
#' @param return_dist Boolean value (default: FALSE). If TRUE, returns
#' the mean and variance of the latent distribution.
#' Otherwise, returns the mean of the latent distribution.
#' @param dataloader An iterator over minibatches of data on which to
#' compute the metric.
#' @param transform_batch Batch to condition on.  If transform_batch is:
#' - NULL, then real observed batch is used.
#' - int, then batch transform_batch is used.
#' - Otherwise based on string
#' @param gene_list Return frequencies of expression for a subset of genes.
#' @param library_size Scale the expression frequencies to a common library size.
#' @param n_samples Number of posterior samples to use for estimation.
#' @param n_samples_overall Number of posterior samples to use for estimation.
#' Overrides n_samples.
#' @param weights Weights to use for sampling. If None, defaults to “uniform”.
#' @param return_mean Whether to return the mean of the samples
#' @param return_numpy Return a ndarray instead of a DataFrame.
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom reticulate r_to_py
#'
#' @return A list that contains the corrected gene expression matrix and the
#' corrected low-dimensional space.
#'
#' @examples
#' sim <- simulate_data(n_genes = 500, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' scvi <- scVIRun(input = sim, batch = "Batch", max_epochs = 1)
#'
scVIRun <- function(input, batch, assay_type = "counts", layer = NULL, labels_key = NULL,
                    size_factor_key = NULL, categorical_covariate_keys = NULL,
                    continuous_covariate_keys = NULL, n_hidden = 128,
                    n_latent = 10, n_layers = 1, dropout_rate = 0.1,
                    dispersion = "gene", gene_likelihood = "zinb",
                    latent_distribution = "normal", max_epochs = NULL,
                    accelerator = "auto", devices = 1, train_size = 0.25,
                    validation_size = NULL, shuffle_set_split = TRUE,
                    load_sparse_tensor = FALSE, batch_size = 128,
                    early_stopping = FALSE, datasplitter_kwargs = NULL,
                    plan_kwargs = NULL, datamodule = NULL,
                    indices = NULL, give_mean = TRUE, mc_samples = 5000,
                    return_dist = FALSE, dataloader = NULL,
                    transform_batch = NULL, gene_list = NULL, library_size = 1,
                    n_samples = 1, n_samples_overall = NULL, weights = NULL,
                    return_mean = TRUE, return_numpy = NULL) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, batch, assay_type = "counts",
                                                 layer = NULL, labels_key = NULL,
                                                 size_factor_key = NULL,
                                                 categorical_covariate_keys = NULL,
                                                 continuous_covariate_keys = NULL,
                                                 n_hidden = 128, n_latent = 10,
                                                 n_layers = 1, dropout_rate = 0.1,
                                                 dispersion = "gene",
                                                 gene_likelihood = "zinb",
                                                 latent_distribution = "normal",
                                                 max_epochs = 400,
                                                 accelerator = "auto",
                                                 devices = 1, train_size = 0.25,
                                                 validation_size = NULL,
                                                 shuffle_set_split = TRUE,
                                                 load_sparse_tensor = FALSE,
                                                 batch_size = 128,
                                                 early_stopping = FALSE,
                                                 datasplitter_kwargs = NULL,
                                                 plan_kwargs = NULL,
                                                 datamodule = NULL,
                                                 indices = NULL, give_mean = TRUE,
                                                 mc_samples = 5000,
                                                 return_dist = FALSE,
                                                 dataloader = NULL,
                                                 transform_batch = NULL,
                                                 gene_list = NULL,
                                                 library_size = 1, n_samples = 1,
                                                 n_samples_overall = NULL,
                                                 weights = NULL,
                                                 return_mean = TRUE,
                                                 return_numpy = NULL) {
    scvi <- reticulate::import("scvi")

    if (is(input, "SingleCellExperiment")) {
      andata <- zellkonverter::SCE2AnnData(input)
      andata$layers[assay_type] <- andata$X
    }
    else if (is(input, "Seurat")) {
      sce <- as.SingleCellExperiment(input)
      andata <- zellkonverter::SCE2AnnData(sce)
      andata$layers[assay_type] <- andata$X
    }
    else {
      andata <- reticulate::r_to_py(input, convert = TRUE)
    }

    arg <- c(list(adata = andata,layer = layer, batch_key = batch,
                  labels_key = labels_key, size_factor_key = size_factor_key,
                  categorical_covariate_keys = categorical_covariate_keys,
                  continuous_covariate_keys = continuous_covariate_keys))
    do.call(scvi$model$SCVI$setup_anndata, arg)

    arg <- c(list(adata = andata, n_hidden = as.integer(n_hidden),
                  n_latent = as.integer(n_latent), n_layers = as.integer(n_layers),
                  dropout_rate = dropout_rate, dispersion = dispersion,
                  gene_likelihood = gene_likelihood,
                  latent_distribution = latent_distribution))
    model_scvi <- do.call(scvi$model$SCVI, arg)

    arg <- c(list(max_epochs = as.integer(max_epochs), accelerator = accelerator,
                  devices = as.integer(devices), train_size = train_size,
                  validation_size = validation_size,
                  shuffle_set_split = shuffle_set_split,
                  load_sparse_tensor = load_sparse_tensor,
                  batch_size = as.integer(batch_size), early_stopping = early_stopping,
                  datasplitter_kwargs = datasplitter_kwargs,
                  plan_kwargs = plan_kwargs, datamodule = datamodule))
    do.call(model_scvi$train, arg)

    latent <- model_scvi$get_latent_representation(adata = andata,
                                                   indices = indices,
                                                   batch_size = as.integer(batch_size),
                                                   give_mean = give_mean,
                                                   mc_samples = as.integer(mc_samples),
                                                   return_dist = return_dist,
                                                   dataloader = dataloader)

    if (is(input, "AnnDataR6")) {
      rownames(latent) <- andata$obs_names
    }
    else {
      rownames(latent) <- colnames(input)
    }
    colnames(latent) <- paste0("scvi_", 1:n_latent)

    corrected <- model_scvi$get_normalized_expression(adata = andata,
                                                      indices = indices,
                                                      transform_batch = transform_batch,
                                                      gene_list = gene_list,
                                                      library_size = as.integer(library_size),
                                                      n_samples = as.integer(n_samples),
                                                      n_samples_overall = n_samples_overall,
                                                      weights = weights,
                                                      batch_size = as.integer(batch_size),
                                                      return_mean = return_mean,
                                                      return_numpy = return_numpy)

    return(list(corrected, latent))
  }, input = input, batch = batch, layer = layer, labels_key = labels_key,
  size_factor_key = size_factor_key,
  categorical_covariate_keys = categorical_covariate_keys,
  continuous_covariate_keys = continuous_covariate_keys,
  n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers,
  dropout_rate = dropout_rate,
  dispersion = dispersion, gene_likelihood = gene_likelihood,
  latent_distribution = latent_distribution,
  max_epochs = max_epochs, accelerator = accelerator, devices = devices,
  train_size = train_size, validation_size = validation_size,
  shuffle_set_split = shuffle_set_split,
  load_sparse_tensor = load_sparse_tensor, batch_size = batch_size,
  early_stopping = early_stopping, datasplitter_kwargs = datasplitter_kwargs,
  plan_kwargs = plan_kwargs, datamodule = datamodule,
  indices = indices, give_mean = give_mean, mc_samples = mc_samples,
  return_dist = return_dist, dataloader = dataloader,
  transform_batch = transform_batch, gene_list = gene_list,
  library_size = as.integer(library_size), n_samples = as.integer(n_samples),
  n_samples_overall = n_samples_overall, weights = weights,
  return_mean = return_mean, return_numpy = return_numpy)
}

