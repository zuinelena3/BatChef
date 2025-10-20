#' Convert the LIGER output
#'
#' Convert the LIGER output into a \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#'
#' @param input  A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param output LIGER output: liger object
#' @param method A string specifying the correction method
#'
#' @import methods
#' @return  A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object.
#' @rdname ligerPost
#'
setGeneric("ligerPost", function(input, output, method)
  standardGeneric("ligerPost"), signature = c("input"))

#' @rdname ligerPost
#' @aliases ligerPost,Seurat,Seurat-method
#' @import methods
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("ligerPost", "Seurat",  function(input, output, method) {
  embedding <- as.data.frame(output@H.norm)
  rownames(embedding) <- colnames(input)

  input[[method]] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                          key = "liger_",
                                          assay = DefaultAssay(input))
  return(input)
})

#' @rdname ligerPost
#' @aliases ligerPost,SingleCellExperiment,SingleCellExperiment-method
#' @import methods
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#'
setMethod("ligerPost", "SingleCellExperiment",  function(input, output, method) {
  embedding <- as.data.frame(output@H.norm)
  rownames(embedding) <- colnames(input)
  reducedDim(input, method) <- embedding
  return(input)
})

#' @rdname ligerPost
#' @aliases ligerPost,AnnDataR6,AnnDataR6-method
#' @importFrom SingleCellExperiment reducedDim
#'
setMethod("ligerPost", "AnnDataR6",  function(input, output, method) {
  embedding <- as.data.frame(output@H.norm)
  rownames(embedding) <- input$obs_names
  input$obsm[[method]] <- as.matrix(embedding)
  return(input)
})

