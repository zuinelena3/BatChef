#' Convert to a LIGER compatible object
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param batch A string specifying the batch variable.
#' @param features Vector of features to use.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @import methods
#' @return A LIGER compatible object
#' @rdname ligerInput
#'
setGeneric("ligerInput", function(input, batch, features, ...)
  standardGeneric("ligerInput"), signature = c("input"))

#' @rdname ligerInput
#' @import methods
#' @importFrom Seurat SplitObject
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,Seurat,Seurat-method
#'
setMethod("ligerInput", "Seurat",  function(input, batch, features, ...) {
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
setMethod("ligerInput", "SingleCellExperiment",  function(input, batch,
                                                          features, ...) {
  stopifnot(batch %in% colnames(colData(input)))
  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  names(ll) <- unique(colData(input)[, batch])

  args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo)
  lo <- scaleNotCenter(object = lo,
                       features = features)
})

#' @rdname ligerInput
#' @import methods
#' @importFrom zellkonverter AnnData2SCE
#' @importFrom SingleCellExperiment colData
#' @importFrom rliger createLiger normalize scaleNotCenter
#' @aliases ligerInput,AnnDataR6,AnnDataR6-method
#'
setMethod("ligerInput", "AnnDataR6",  function(input, batch, features, ...) {
  stopifnot(batch %in% colnames(input$obs))
  input <- AnnData2SCE(input)

  ll <- lapply(unique(colData(input)[, batch]), function(i) input[, colData(input)[, batch] == i])
  names(ll) <- unique(colData(input)[, batch])

  args <- c(list(rawData = ll), capture_params(createLiger, list(...)))
  lo <- do.call(createLiger, args)

  lo <- normalize(object = lo)
  lo <- scaleNotCenter(object = lo, features = features)
})
