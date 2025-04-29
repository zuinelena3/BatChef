#' Title
#'
#' @param input input
#' @param output output
#' @param batch batch
#' @param return_dimred return_dimred
#' @param method method
#'
#' @import methods
#' @rdname scanoramaPost
#'
setGeneric("scanoramaPost", function(input, output, batch, return_dimred, method)
  standardGeneric("scanoramaPost"), signature = c("input"))

#' @rdname scanoramaPost
#' @aliases scanoramaPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject CreateDimReducObject DefaultAssay
#'
setMethod("scanoramaPost", "Seurat",  function(input, output, batch, return_dimred, method) {
  ll <- SplitObject(input, split.by = batch)
  if (return_dimred == TRUE) {
    corrected_mat <- t(do.call(rbind, output[[2]]))
    colnames(corrected_mat) <- unlist(lapply(ll, function(x) colnames(x)))
    rownames(corrected_mat) <- output[[3]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]

    input[[paste0(method, "_mat")]] <- CreateAssayObject(data = as.matrix(corrected_mat))
    input[[method]] <- CreateDimReducObject(embeddings = embedding,
                                            key = "scanorama_", assay = DefaultAssay(input))
    return(input)
  }
  else {
    corrected_mat <- t(do.call(rbind, output[[1]]))
    colnames(corrected_mat) <- unlist(lapply(ll, function(x) colnames(x)))
    rownames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    input[[paste0(method, "_mat")]] <- CreateAssayObject(data = as.matrix(corrected_mat))
    return(input)
  }
})

#' @rdname scanoramaPost
#' @aliases scanoramaPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SummarizedExperiment assay<- assay
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("scanoramaPost", "SingleCellExperiment",  function(input, output, batch, return_dimred, method) {
  ll <- lapply(unique(colData(input)[, batch]), function(x) input[, colData(input)[, batch] == x])
  if (return_dimred == TRUE) {
    corrected_mat <- t(do.call(rbind, output[[2]]))
    colnames(corrected_mat) <- unlist(lapply(ll, function(x) colnames(x)))
    rownames(corrected_mat) <- output[[3]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]

    assay(input, method) <- as.matrix(corrected_mat)
    reducedDim(input, method) <- embedding
    return(input)
  }
  else {
    corrected_mat <- t(do.call(rbind, output[[1]]))
    colnames(corrected_mat) <- unlist(lapply(ll, function(x) colnames(x)))
    rownames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    assay(input, method) <- as.matrix(corrected_mat)
    return(input)
  }
})

#' @rdname scanoramaPost
#' @aliases scanoramaPost,AnnDataR6,AnnDataR6-method
#'
setMethod("scanoramaPost", "AnnDataR6",  function(input, output, batch, return_dimred, method) {
  batches <- unique(input$obs[batch])[, 1]
  ll <- lapply(batches, function(i) input[input$obs[batch] == i, ])

  if (return_dimred == TRUE) {
    corrected_mat <- do.call(rbind, output[[2]])
    rownames(corrected_mat) <- unlist(lapply(ll, function(x) x$obs_names))
    colnames(corrected_mat) <- output[[3]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- rownames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), rownames(input))), ]

    input$layers[method] <- corrected_mat
    input$obsm[[method]] <- embedding
    return(input)
  }
  else {
    corrected_mat <- do.call(rbind, output[[1]])
    rownames(corrected_mat) <- unlist(lapply(ll, function(x) x$obs_names))
    colnames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat), colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat), rownames(input))), ]

    input$layers[method] <- corrected_mat
    return(input)
  }
})
