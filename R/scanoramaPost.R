#' Convert the Scanorama output
#'
#' Convert the Scanorama output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param list List version of input.
#' @param output Scanorama output: a list which contains the transposed corrected
#' gene expression matrix and the corrected low-dimensional space.
#' @param return_dimred A logical to returning integrated low-dimesional embeddings.
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname scanoramaPost
#'
setGeneric("scanoramaPost", function(input, list, output, return_dimred, method) {
  standardGeneric("scanoramaPost")
}, signature = c("input"))

#' @rdname scanoramaPost
#' @aliases scanoramaPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateAssayObject CreateDimReducObject DefaultAssay
#'
setMethod("scanoramaPost", "Seurat", function(input, list, output,
                                              return_dimred, method) {
  n <- length(list)
  n_batch <- n / 2

  indx <- length(output)
  mat <- t(as.matrix(do.call(rbind, output[[indx - 1]])))
  colnames(mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
  rownames(mat) <- output[[indx]]

  mat <- mat[
    order(match(rownames(mat), rownames(input))),
    order(match(colnames(mat), colnames(input)))
  ]

  input[[paste0(method, "_mat")]] <- CreateAssayObject(counts = mat)

  if (return_dimred == TRUE) {
    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]
    input[[method]] <- CreateDimReducObject(
      embeddings = embedding,
      key = "scanorama_",
      assay = DefaultAssay(input)
    )
    return(input)
  } else {
    return(input)
  }
})

#' @rdname scanoramaPost
#' @aliases scanoramaPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SummarizedExperiment assay<- assay
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("scanoramaPost", "SingleCellExperiment", function(input, list, output,
                                                            return_dimred,
                                                            method) {
  n <- length(list)
  n_batch <- n / 2

  indx <- length(output)
  mat <- t(as.matrix(do.call(rbind, output[[indx - 1]])))
  colnames(mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
  rownames(mat) <- output[[indx]]

  mat <- mat[
    order(match(rownames(mat), rownames(input))),
    order(match(colnames(mat), colnames(input)))
  ]

  assay(input, method) <- mat

  if (return_dimred == TRUE) {
    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- colnames(mat)
    embedding <- embedding[order(match(rownames(embedding), colnames(input))), ]
    reducedDim(input, method) <- embedding
    return(input)
  } else {
    return(input)
  }
})

#' @rdname scanoramaPost
#' @aliases scanoramaPost,AnnDataR6,AnnDataR6-method
#'
setMethod("scanoramaPost", "AnnDataR6", function(input, list, output,
                                                 return_dimred, method) {
  n <- length(list)
  n_batch <- n / 2

  indx <- length(output)
  mat <- as.matrix(do.call(rbind, output[[indx - 1]]))
  rownames(mat) <- unlist(lapply(list[1:n_batch], function(x) rownames(x)))
  colnames(mat) <- output[[indx]]

  mat <- mat[
    order(match(rownames(mat), rownames(input))),
    order(match(colnames(mat), colnames(input)))
  ]

  input$layers[method] <- mat

  if (return_dimred == TRUE) {
    embedding <- do.call(rbind, output[[1]])
    colnames(embedding) <- paste0("scanorama_", seq(1, dim(embedding)[2]))
    rownames(embedding) <- rownames(mat)
    embedding <- embedding[order(match(rownames(embedding), rownames(input))), ]
    input$obsm[[method]] <- embedding
    return(input)
  } else {
    return(input)
  }
})
