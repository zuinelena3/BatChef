#' Batch effect correction methods
#'
#' A common interface for single-cell batch correction methods.
#'
#' @param obj a \linkS4class{SingleCellExperiment} object or list containing single-cell gene expression matrices
#' @param batch a string specifying the batch
#' @param assay a string specifying the assay to use for correction
#' @param hvgs a vector specifying which features to use for correction
#' @param dims number of dimensions to use for dimension reduction
#' @param reduction a string specifying the dimension reduction to use for correction
#' @param anchor a string specifying the anchor integration type for Seurat methods
#' @param k_anchor number of anchors
#' @param genelist list containing the genes in each batch
#' @param cell_type string specifying the cell-type labels
#' @param METHOD a \code{MethodParam} object specifying the batch correction method
#' @param alt_out alternative output: a \code{\link{AltOutput}} class
#'
#' @import methods
#' @rdname scIntegration
#' @importFrom limma removeBatchEffect
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData assays
#'
setMethod("scIntegration", "limmaMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                    dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                    genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  if(!("logcounts" %in% names(assays(data)))) {
    print("Error: the object has to contain log-expression values!")
  }

  out <- removeBatchEffect(x = logcounts(obj), batch = colData(obj)[, batch])

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = out,
               embedding = NULL,
               meta = data.frame(cell_id = colnames(out),
                                 batch = colData(obj)[, batch],
                                 cell_type = colData(obj)[, cell_type]))
  }
  else return(out)
})

#' @rdname scIntegration
#' @importFrom SummarizedExperiment assays colData
#' @importFrom singleCellTK runComBatSeq
#'
setMethod("scIntegration", "combatMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                    dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                    genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  out <- runComBatSeq(inSCE = obj, useAssay = assay, batch = batch, covariates = cell_type)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = assays(out)$ComBatSeq,
               embedding = NULL,
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = colData(obj)[, batch],
                                       cell_type = colData(obj)[, cell_type])))
    return(res)
  }

  else return(out)
})

#' @rdname scIntegration
#' @param obj A SingleCellExperiment object
#' @param anchor a string specifying the anchors finding type (cca, rpca, jpca, rlsi)
#' @param k_anchor number of anchors (default: k_anchor = 5)
#'
#' @importFrom Seurat FindIntegrationAnchors IntegrateData GetAssayData SplitObject VariableFeatures<- CreateDimReducObject
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SingleCellExperiment counts logcounts colData reducedDim
#' @importFrom scater runPCA
#'
setMethod("scIntegration", "seuratv3Method", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                      dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                      genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  options(Seurat.object.assay.version = "v3", future.globals.maxSize = 1000 * 1024^2)

  so <- CreateSeuratObject(counts = counts(obj))
  so@assays$RNA$data <- logcounts(obj)
  so@meta.data <- cbind(so@meta.data, as.data.frame(colData(obj)[, c(batch, cell_type)]))
  VariableFeatures(so) <- rownames(so)

  if(is.null(reduction)) {
    set.seed(333)
    obj <- runPCA(obj, subset_row = hvgs, ncomponent = dims)

    so@reductions[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(obj, "PCA"),
                                                                loadings = attr(reducedDim(obj, "PCA"), "rotation"),
                                                                key = "pca_", assay = 'RNA')
  }

  so@reductions[[tolower(reduction)]] <- CreateDimReducObject(embeddings = reducedDim(obj, reduction),
                                                              loadings = attr(reducedDim(obj, reduction), "rotation"),
                                                              key = paste0(tolower(reduction), '_'), assay = 'RNA')

  so_ll <- SplitObject(so, split.by = batch)

  anchorset <- FindIntegrationAnchors(object.list = so_ll, reduction = anchor, anchor.features = hvgs, dims = 1:dims,  k.anchor = k_anchor, verbose = FALSE)
  out <- IntegrateData(anchorset = anchorset, verbose = FALSE)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = GetAssayData(object = out, assay = "integrated"),
               embedding = NULL,
               meta = data.frame(cbind(cell_id = colnames(out),
                                       batch = out@meta.data[, batch],
                                       cell_type = out@meta.data[, cell_type])))

    return(res)
  }
  else return(out)
})

#' @rdname scIntegration
#' @param anchor a string specifying the anchors finding type (CCAIntegration, RPCAIntegration, HarmonyIntegration, JointPCAIntegration)
#' @param k_anchor number of anchors (default: k_anchor = 5)
#'
#' @importFrom Seurat ScaleData CreateDimReducObject IntegrateLayers
#' @importFrom SeuratObject CreateAssay5Object JoinLayers CreateSeuratObject
#' @importFrom SingleCellExperiment reducedDim counts logcounts
#' @importFrom SummarizedExperiment colData
#'
setMethod("scIntegration", "seuratv5Method", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                      dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                      genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  options(Seurat.object.assay.version = "v5", future.globals.maxSize = 1000 * 1024^2)

  out <- CreateAssay5Object(counts = counts(obj), data = logcounts(obj))
  out <- CreateSeuratObject(out)
  out@meta.data <- cbind(out@meta.data, as.data.frame(colData(obj)[, c(batch, cell_type)]))

  if(is.null(reduction)) {
    set.seed(333)
    obj <- runPCA(obj, subset_row = hvgs, ncomponent = dims)

    so@reductions[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(obj, "PCA"),
                                                   loadings = attr(reducedDim(obj, "PCA"), "rotation"),
                                                   key = "pca_", assay = 'RNA')
  }

  out[["RNA"]] <- split(out[["RNA"]], f = out@meta.data[, batch])
  out@reductions[[tolower(reduction)]] <- CreateDimReducObject(embeddings = reducedDim(obj, reduction),
                                                               loadings = attr(reducedDim(obj, reduction), "rotation"),
                                                               key = paste0(reduction , "_"), assay = 'RNA')
  out <- ScaleData(out, verbose = FALSE)

  out <- IntegrateLayers(object = out, method = anchor, orig.reduction = tolower(reduction),
                         new.reduction = "integrated", features = hvgs, k.anchor = k_anchor, verbose = FALSE)

  out[["RNA"]] <- JoinLayers(out[["RNA"]])

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = NULL,
               embedding = out@reductions$integrated@cell.embeddings,
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = out@meta.data[, batch],
                                       cell_type = out@meta.data[, cell_type])))

    return(res)
  }
  else return(out)
})
