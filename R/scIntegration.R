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
  if(!("logcounts" %in% names(assays(obj)))) {
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
#' @importFrom sva ComBat_seq
#'
setMethod("scIntegration", "combatMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                    dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                    genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  out <- ComBat_seq(counts = as.matrix(counts(obj)), batch =  colData(obj)[, batch])

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = out,
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
  options(Seurat.object.assay.version = "v3")

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
  options(Seurat.object.assay.version = "v5")

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

#' @rdname scIntegration
#' @param obj a SingleCellExperiment object  or a list containing SingleCellExperiment objects
#' @param batch a string specifying the batch. batch = NULL when obj is a list
#'
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assays colData
#' @importFrom Seurat as.sparse
#'
setMethod("scIntegration", "fastMNNMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  out <- fastMNN(obj = obj, batch = colData(obj)[batch][, 1], subset.row = hvgs, d = dims)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = as.sparse(as.matrix(assays(out)$reconstructed)),
               embedding = reducedDim(out, "corrected"),
               meta = data.frame(cell_id = NA,
                                 batch = NA,
                                 cell_type = NA))

    if (is(obj, "list")) {
      res@meta <- data.frame(cbind(cell_id = unlist(lapply(obj, function(x) colnames(x))),
                                   batch = unlist(lapply(obj, function(x) colData(x)[, batch])),
                                   cell_type = unlist(lapply(obj, function(x) colData(x)[, cell_type]))))
    }

    else {
      res@meta <- data.frame(cell_id = colnames(obj),
                             batch = colData(obj)[, batch],
                             cell_type = colData(obj)[, cell_type])
    }
    return(res)
  }

  else return(out)
})

#' @rdname scIntegration
#' @importFrom harmony RunHarmony
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assays colData
#' @importFrom Seurat as.sparse
#'
setMethod("scIntegration", "harmonyMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  out <- RunHarmony(obj, batch)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = NULL,
               embedding = reducedDim(out, "HARMONY"),
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = out$batch,
                                       cell_type = colData(obj)[, cell_type])))
    return(res)
  }

  else return(out)
})

#' @rdname scIntegration
#' @param obj a list of "matrix" "array" objects
#' @param genelist a list of genes
#' @param hvgs number of highly variable genes
#'
#' @importFrom reticulate import
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#'
setMethod("scIntegration", "scanoramaMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                       dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                       genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))
  out <- basiliskRun(proc = proc, fun = function(obj, genelist, hvgs) {
    scanorama <- import("scanorama")
    method <- scanorama$correct(obj, genelist, return_dimred = TRUE, return_dense = TRUE, verbose = FALSE, hvg = hvgs, dimred = as.integer(dims))
  }, obj = obj, genelist = genelist, hvgs = hvgs)

  if (alt_out == TRUE) {
    corrected_counts <- t(do.call(rbind, out[[2]]))
    colnames(corrected_counts) <- unlist(lapply(obj, function(x) rownames(x)))
    rownames(corrected_counts) <- out[[3]]

    embedding <- do.call(rbind, out[[1]])
    colnames(embedding) <- paste0("Scanorama_", seq(1, dims))
    rownames(embedding) <- colnames(corrected_counts)

    res <- new("AltOutput", corrected = corrected_counts,
               embedding = embedding,
               meta = data.frame(cbind(cell_id = colnames(corrected_counts),
                                       batch = batch,
                                       cell_type = cell_type)))
    return(res)
  }

  else return(out)
})

#' @rdname scIntegration
#' @param obj A SingleCellExperiment
#'
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom Rtsne Rtsne_neighbors
#'
setMethod("scIntegration", "bbknnMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                   dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                   genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  proc <- basiliskStart(py_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, batch, reduction) {
    bbknn <- import("bbknn")
    sc <- import("scanpy")
    anndata <- import("anndata")

    adata <- SCE2AnnData(obj)
    adata$X <- adata$layers["logcounts"]
    adata$obs$batch <- colData(obj)[, batch]
    adata$obsm["X_pca"] <- reducedDim(obj, reduction)

    if (adata$n_obs > 100000) {
      neighbors_within_batch = 25
    } else neighbors_within_batch = 3

    bbknn$bbknn(adata, batch_key = batch, neighbors_within_batch = as.integer(neighbors_within_batch))
    sc$tl$umap(adata)
    bbknn_umap <- adata$obsm[["X_umap"]]
    reducedDim(obj, "UMAP_bbknn") <- bbknn_umap

    bbknn <- knn_index_dist(dist = adata$obsp[["distances"]])
    perplexity <- ncol(x = bbknn$idx) - 1
    tsne <- Rtsne_neighbors(
      index = bbknn$idx,
      distance = bbknn$dist,
      perplexity = perplexity
    )$Y
    colnames(x = tsne) <- paste0("tSNE_", c(1, 2))
    rownames(x = tsne) <- colnames(obj)
    reducedDim(obj, "TSNE_bbknn") <- tsne
    return(obj)
  }, obj = obj, batch = batch, reduction = reduction)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = NULL,
               embedding = list(tSNE = reducedDim(out, "TSNE_bbknn"), UMAP = reducedDim(out, "UMAP_bbknn")),
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = colData(obj)[, batch],
                                       cell_type = colData(obj)[, cell_type])))
    return(res)
  }
  else return(out)
})

#' @rdname scIntegration
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom zellkonverter SCE2AnnData
#' @importFrom reticulate import
#'
setMethod("scIntegration", "scVIMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                  dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                  genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  proc <- basiliskStart(scvi_env)
  on.exit(basiliskStop(proc))

  out <- basiliskRun(proc = proc, fun = function(obj, batch, assay, dims) {
    scvi <- import("scvi")

    andata <- SCE2AnnData(obj)
    andata$layers["counts"] <- andata$X
    scvi$model$SCVI$setup_anndata(andata, layer = assay, batch_key = batch)
    model_scvi <- scvi$model$SCVI(andata, n_latent = as.integer(dims))
    model_scvi$train()

    andata$obsm["X_scVI"] <- model_scvi$get_latent_representation()
    return(andata)
  }, obj = obj, batch = batch, assay = assay, dims = dims)

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = NULL,
               embedding = out$obsm["X_scVI"],
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = colData(obj)[, batch],
                                       cell_type = colData(obj)[, cell_type])))
    return(res)
  }
  else return(out)
})

#' @rdname scIntegration
#' @param genelist negative controls
#'
#' @importFrom scMerge scMerge2
#' @importFrom SingleCellExperiment logcounts colData
#'
setMethod("scIntegration", "scMergeMethod", function(obj, batch = NULL, assay = NULL, hvgs = NULL,
                                                     dims = NULL, reduction = NULL, anchor = NULL, k_anchor = NULL,
                                                     genelist = NULL, cell_type = NULL, METHOD, alt_out = FALSE) {
  if(is.null(cell_type)){
    out <- scMerge2(exprsMat = logcounts(obj), batch = colData(obj)[, batch], ctl = genelist,
                    verbose = FALSE)
  }
  else {
    out <- scMerge2(exprsMat = logcounts(obj), batch = colData(obj)[, batch], ctl = genelist,
                    cellTypes = colData(obj)[, cell_type], verbose = FALSE)
  }

  if (alt_out == TRUE) {
    res <- new("AltOutput", corrected = out$newY,
               embedding = NULL,
               meta = data.frame(cbind(cell_id = colnames(obj),
                                       batch = colData(obj)[, batch],
                                       cell_type = colData(obj)[, cell_type])))
    return(res)
  }
  else return(out)
})
