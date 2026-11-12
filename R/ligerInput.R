#' Convert to a LIGER compatible object
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param features Vector of features to use.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @import methods
#' @return A LIGER compatible object
#' @rdname ligerInput
#'
setGeneric("ligerInput", function(input, batch, features, ...) {
  standardGeneric("ligerInput")
}, signature = c("input"))

#' @rdname ligerInput
#' @import methods
#' @importFrom Seurat SplitObject
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,Seurat,Seurat-method
#'
setMethod("ligerInput", "Seurat", function(input, batch, features, ...) {
  stopifnot(batch %in% colnames(input[[]]))
  so_ll <- SplitObject(input, split.by = batch)
  names(so_ll) <- unique(input[[batch]][, 1])

  args <- c(list(rawData = so_ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo)
  lo <- scaleNotCenter(object = lo, features = features)
})

#' @rdname ligerInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("ligerInput", "SingleCellExperiment", function(input, batch,
                                                         features, ...) {
  stopifnot(batch %in% colnames(colData(input)))
  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  names(ll) <- unique(colData(input)[, batch])

  args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo)
  lo <- scaleNotCenter(
    object = lo,
    features = features
  )
})

#' @rdname ligerInput
#' @import methods
#' @importFrom SingleCellExperiment colData reducedDim<- reducedDim
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @importFrom reticulate r_to_py
#' @aliases ligerInput,InMemoryAnnData,InMemoryAnnData-method
#'
setMethod("ligerInput", "ANY", function(input, batch, features, ...) {
  stopifnot(batch %in% colnames(input$obs))

  if (is(input, "InMemoryAnnData")) {
    input <- input$as_SingleCellExperiment(x_mapping = "counts",
                                           assays_mapping = c(logcounts = "normalized"),
                                           reducedDims_mapping = list(PCA = c(sampleFactors = "X_pca", featureLoadings = "PCs")))
    reducedDim(input, "PCA") <- as.matrix(reducedDim(input, "PCA"))
    ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
    names(ll) <- unique(colData(input)[, batch])

    args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
    lo <- do.call(createLiger, args)

    lo <- normalize(object = lo)
    lo <- scaleNotCenter(object = lo, features = features)
  } else {
    message("The input object must be of class 'InMemoryAnnData'.")
  }
})
