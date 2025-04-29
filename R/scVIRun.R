#' scVI method
#'
#' @param input input
#' @param batch batch
#' @param layer layer
#' @param labels_key labels_key
#' @param size_factor_key size_factor_key
#' @param categorical_covariate_keys categorical_covariate_keys
#' @param continuous_covariate_keys continuous_covariate_keys
#' @param n_hidden n_hidden
#' @param n_latent n_latent
#' @param n_layers n_layers
#' @param dropout_rate dropout_rate
#' @param dispersion dispersion
#' @param gene_likelihood gene_likelihood
#' @param latent_distribution latent_distribution
#' @param max_epochs max_epochs
#' @param accelerator accelerator
#' @param devices devices
#' @param train_size train_size
#' @param validation_size validation_size
#' @param shuffle_set_split shuffle_set_split
#' @param load_sparse_tensor load_sparse_tensor
#' @param batch_size batch_size
#' @param early_stopping early_stopping
#' @param datasplitter_kwargs datasplitter_kwargs
#' @param plan_kwargs plan_kwargs
#' @param datamodule datamodule
#' @param indices indices
#' @param give_mean give_mean
#' @param mc_samples mc_samples
#' @param return_dist return_dist
#' @param dataloader dataloader
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom zellkonverter SCE2AnnData
#'
scVIRun <- function(input, batch, layer = NULL, labels_key = NULL, size_factor_key = NULL, categorical_covariate_keys = NULL, continuous_covariate_keys = NULL,
                    n_hidden = 128, n_latent = 10, n_layers = 1, dropout_rate = 0.1,
                    dispersion = "gene", gene_likelihood = "zinb", latent_distribution = "normal",
                    max_epochs = NULL, accelerator = "auto", devices = 1, train_size = 0.25, validation_size = NULL, shuffle_set_split = TRUE,
                    load_sparse_tensor = FALSE, batch_size = 128, early_stopping = FALSE, datasplitter_kwargs = NULL, plan_kwargs = NULL, datamodule = NULL,
                    indices = NULL, give_mean = TRUE, mc_samples = 5000, return_dist = FALSE, dataloader = NULL) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(input, batch, layer = NULL, labels_key = NULL, size_factor_key = NULL, categorical_covariate_keys = NULL, continuous_covariate_keys = NULL,
                                                 n_hidden = 128, n_latent = 10, n_layers = 1, dropout_rate = 0.1,
                                                 dispersion = "gene", gene_likelihood = "zinb", latent_distribution = "normal",
                                                 max_epochs = 400, accelerator = "auto", devices = 1, train_size = 0.25, validation_size = NULL, shuffle_set_split = TRUE,
                                                 load_sparse_tensor = FALSE, batch_size = 128, early_stopping = FALSE, datasplitter_kwargs = NULL, plan_kwargs = NULL, datamodule = NULL,
                                                 indices = NULL, give_mean = TRUE, mc_samples = 5000, return_dist = FALSE, dataloader = NULL) {
    scvi <- import("scvi")

    andata <- SCE2AnnData(input)
    andata$layers["counts"] <- andata$X

    arg <- c(list(adata = andata,layer = layer, batch_key = batch, labels_key = labels_key, size_factor_key = size_factor_key,
                  categorical_covariate_keys = categorical_covariate_keys, continuous_covariate_keys = continuous_covariate_keys))
    do.call(scvi$model$SCVI$setup_anndata, arg)

    arg <- c(list(adata = andata, n_hidden = as.integer(n_hidden), n_latent = as.integer(n_latent), n_layers = as.integer(n_layers),
                  dropout_rate = dropout_rate, dispersion = dispersion, gene_likelihood = gene_likelihood, latent_distribution = latent_distribution))
    model_scvi <- do.call(scvi$model$SCVI, arg)

    arg <- c(list(max_epochs = as.integer(max_epochs), accelerator = accelerator, devices = as.integer(devices), train_size = train_size, validation_size = validation_size,
                  shuffle_set_split = shuffle_set_split, load_sparse_tensor = load_sparse_tensor, batch_size = as.integer(batch_size), early_stopping = early_stopping,
                  datasplitter_kwargs = datasplitter_kwargs, plan_kwargs = plan_kwargs, datamodule = datamodule))
    do.call(model_scvi$train, arg)

    andata$obsm["X_scVI"] <- model_scvi$get_latent_representation(adata = andata, indices = indices, give_mean = give_mean, mc_samples = as.integer(mc_samples), return_dist = return_dist, dataloader = dataloader)
    return(andata)
  }, input = input, batch = batch, layer = layer, labels_key = labels_key, size_factor_key = size_factor_key, categorical_covariate_keys = categorical_covariate_keys, continuous_covariate_keys = continuous_covariate_keys,
  n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers, dropout_rate = dropout_rate,
  dispersion = dispersion, gene_likelihood = gene_likelihood, latent_distribution = latent_distribution,
  max_epochs = max_epochs, accelerator = accelerator, devices = devices, train_size = train_size, validation_size = validation_size, shuffle_set_split = shuffle_set_split,
  load_sparse_tensor = load_sparse_tensor, batch_size = batch_size, early_stopping = early_stopping, datasplitter_kwargs = datasplitter_kwargs, plan_kwargs = plan_kwargs, datamodule = datamodule,
  indices = indices, give_mean = give_mean, mc_samples = mc_samples, return_dist = return_dist, dataloader = dataloader)
}
