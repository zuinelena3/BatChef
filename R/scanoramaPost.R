#' Convert the Scanorama output
#'
#' Convert the Scanorama output into a \linkS4class{SingleCellExperiment},
#' \linkS4class{Seurat} or `AnnData` object.
#'
#' @param input A \linkS4class{SingleCellExperiment}, \linkS4class{Seurat} or
#' `AnnData` object can be supplied.
#' @param list List version of input.
#' @param output Scanorama output: a list which contains the transposed corrected
#' gene expression matrix and the corrected low-dimensional space.
#' @param return_dimred A logical to returning integrated low-dimesional embeddings.
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return A \linkS4class{SingleCellExperiment},
#' \linkS4class{Seurat} or `AnnData` object.
#' @rdname scanoramaPost
#'
setGeneric("scanoramaPost", function(input, list, output, return_dimred, method)
  standardGeneric("scanoramaPost"), signature = c("input"))

#' @rdname scanoramaPost
#' @aliases scanoramaPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject CreateDimReducObject DefaultAssay
#'
setMethod("scanoramaPost", "Seurat", function(input, list, output,
                                              return_dimred, method) {
  n <- length(list)
  n_batch <- n/2
  if (return_dimred == TRUE) {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[2]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[3]]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]

    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

    input[[paste0(method, "_mat")]] <- CreateAssayObject(data = as.matrix(corrected_mat))
    input[[method]] <- CreateDimReducObject(embeddings = embedding,
                                            key = "scanorama_",
                                            assay = DefaultAssay(input))
    return(input)
  }
  else {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[1]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

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
setMethod("scanoramaPost", "SingleCellExperiment",  function(input, list, output,
                                                             return_dimred,
                                                             method) {
  n <- length(list)
  n_batch <- n/2
  if (return_dimred == TRUE) {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[2]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[3]]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]

    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

    assay(input, method) <- as.matrix(corrected_mat)
    reducedDim(input, method) <- embedding
    return(input)
  }
  else {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[1]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

    assay(input, method) <- as.matrix(corrected_mat)
    return(input)
  }
})

#' @rdname scanoramaPost
#' @aliases scanoramaPost,AnnDataR6,AnnDataR6-method
#'
setMethod("scanoramaPost", "AnnDataR6",  function(input, list, output,
                                                  return_dimred, method) {
  n <- length(list)
  n_batch <- n/2

  if (return_dimred == TRUE) {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[2]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[3]]

    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(corrected_mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]

    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

    input$layers[method] <- corrected_mat
    input$obsm[[method]] <- embedding
    return(input)
  }
  else {
    corrected_mat <- t(as.matrix(do.call(rbind, output[[1]])))
    colnames(corrected_mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
    rownames(corrected_mat) <- output[[2]]
    corrected_mat <- corrected_mat[, order(match(colnames(corrected_mat),
                                                 colnames(input)))]
    corrected_mat <- corrected_mat[order(match(rownames(corrected_mat),
                                               rownames(input))), ]

    input$layers[method] <- corrected_mat
    return(input)
  }
})
