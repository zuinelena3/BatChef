#' Convert to a Scanorama compatible object
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param batch A string specifying the batch for each cell.
#' @param assay_type A string specifying the assay.
#'
#' @import methods
#' @return List that contains the expression matrices for each batches and
#'  the genes.
#' @rdname scanoramaInput
#'
setGeneric("scanoramaInput", function(input, batch, assay_type = NULL) {
  standardGeneric("scanoramaInput")
}, signature = c("input"))

#' @rdname scanoramaInput
#' @import methods
#' @importFrom Seurat DefaultAssay SplitObject
#' @importFrom SeuratObject LayerData
#' @aliases anndatasInput,Seurat,Seurat-method
#'
setMethod("scanoramaInput", "Seurat", function(input, batch, assay_type) {
  stopifnot(batch %in% colnames(input[[]]))
  stopifnot("Error: assay_type has to be specified!" = !is.null(assay_type))

  ll <- SplitObject(input, split.by = batch)
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(ll)) {
    assaylist[[i]] <- t(as.matrix(LayerData(
      object = ll[[i]],
      assay = DefaultAssay(ll[[i]]),
      layer = assay_type
    )))
    genelist[[i]] <- rownames(ll[[i]])
  }
  return(c(assaylist, genelist))
})

#' @rdname scanoramaInput
#' @import methods
#' @importFrom SummarizedExperiment colData assay
#' @aliases anndatasInput,SingleCellExperiment,SingleCellExperiment-method
#'
setMethod("scanoramaInput", "SingleCellExperiment", function(input, batch,
                                                             assay_type) {
  stopifnot(batch %in% colnames(colData(input)))
  stopifnot("Error: assay_type has to be specified!" = !is.null(assay_type))

  ll <- lapply(unique(colData(input)[, batch]), function(x) input[, colData(input)[, batch] == x])
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(ll)) {
    assaylist[[i]] <- t(as.matrix(assay(ll[[i]], assay_type)))
    genelist[[i]] <- rownames(ll[[i]])
  }
  return(c(assaylist, genelist))
})

#' @rdname scanoramaInput
#' @import methods
#' @aliases anndatasInput,AnnDataR6,AnnDataR6-method
#'
setMethod("scanoramaInput", "AnnDataR6", function(input, batch, assay_type) {
  stopifnot(batch %in% colnames(input$obs))
  stopifnot("Error: assay_type has to be specified!" = !is.null(assay_type))

  batches <- unique(input$obs[batch])[, 1]
  ll <- lapply(batches, function(i) input[input$obs[batch] == i, ])

  assaylist <- list()
  genelist <- list()

  for (i in seq_along(ll)) {
    assaylist[[i]] <- ll[[i]]$layers[assay_type]
    genelist[[i]] <- ll[[i]]$var_names
  }

  return(c(assaylist, genelist))
})
