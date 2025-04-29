#' Title
#'
#' @param input input
#' @param batch batch
#' @param features features
#' @param useDatasets useDatasets
#' @param verbose verbose
#' @param format.type format.type
#' @param remove.missing remove.missing
#' @param ... ...
#'
#' @import methods
#' @rdname ligerInput
#'
setGeneric("ligerInput", function(input, batch, features, useDatasets = NULL,
                                  verbose = TRUE, format.type = NULL, remove.missing = NULL, ...)
  standardGeneric("ligerInput"), signature = c("input"))

#' @rdname ligerInput
#' @import methods
#' @importFrom Seurat SplitObject
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,Seurat,Seurat-method
#'
setMethod("ligerInput", "Seurat",  function(input, batch, features, useDatasets = NULL,
                                            verbose = TRUE, format.type = NULL, remove.missing = NULL, ...) {
  stopifnot(batch %in% colnames(input[[]]))
  so_ll <- SplitObject(input, split.by = batch)
  names(so_ll) <- unique(input[[batch]][, 1])

  args <- c(list(rawData = so_ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo, useDatasets = useDatasets, verbose = verbose, format.type = format.type, remove.missing = remove.missing)
  lo <- scaleNotCenter(object = lo, useDatasets = useDatasets, features = features, verbose = verbose, remove.missing = remove.missing)
})

#' @rdname ligerInput
#' @import methods
#' @importFrom SingleCellExperiment colData
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("ligerInput", "SingleCellExperiment",  function(input, batch, features, useDatasets = NULL,
                                                          verbose = TRUE, format.type = NULL, remove.missing = NULL, ...) {
  stopifnot(batch %in% colnames(colData(input)))
  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  names(ll) <- unique(colData(input)[, batch])

  args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo, useDatasets = useDatasets, verbose = verbose, format.type = format.type, remove.missing = remove.missing)
  lo <- scaleNotCenter(object = lo, useDatasets = useDatasets, features = features, verbose = verbose, remove.missing = remove.missing)
})

#' @rdname ligerInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom SingleCellExperiment colData
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,AnnDataR6,AnnDataR6-method
#'
setMethod("ligerInput", "AnnDataR6",  function(input, batch, features, useDatasets = NULL,
                                               verbose = TRUE, format.type = NULL, remove.missing = NULL, ...) {
  stopifnot(batch %in% colnames(input$obs))
  input <- AnnData2SCE(input)

  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  names(ll) <- unique(colData(input)[, batch])

  args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo, useDatasets = useDatasets, verbose = verbose, format.type = format.type, remove.missing = remove.missing)
  lo <- scaleNotCenter(object = lo, useDatasets = useDatasets, features = features, verbose = verbose, remove.missing = remove.missing)
})
