#' BatChefParams class
#'
#' @export
#' @importClassesFrom S4Vectors SimpleList
#' @import methods
#'
#' @rdname BatChefParams
setClass("BatChefParams", contains = "SimpleList")

#' @export
#' @rdname BatChefParams
setClass("LimmaParams",
  contains = "BatChefParams",
  slots = c(assay_type = "character", extra = "list")
)

#' @export
#' @rdname BatChefParams
setClass("CombatParams",
  contains = "BatChefParams",
  slots = c(assay_type = "character", extra = "list")
)

#' @export
#' @rdname BatChefParams
setClass("SeuratV3Params",
  contains = "BatChefParams",
  slots = c(
    features = "character", pca_name = "ANY", assay = "ANY",
    reference = "ANY", anchor_features = "numeric",
    scale = "logical", normalization_method = "character",
    sct_clip_range = "ANY", reduction = "character",
    l2_norm = "logical", dims = "numeric", k_anchor = "numeric",
    k_filter = "numeric", k_score = "numeric",
    max_features = "numeric", nn_method = "character",
    n_trees = "numeric", eps = "numeric",
    verbose = "logical", new_assay_name = "character",
    features_to_integrate = "ANY", k_weight = "numeric",
    weight_reduction = "ANY", sd_weight = "numeric",
    sample_tree = "ANY", preserve_order = "logical"
  )
)

#' @export
#' @rdname BatChefParams
setClass("SeuratV5Params",
  contains = "BatChefParams",
  slots = c(
    pca_name = "ANY", method = "character",
    orig_reduction = "character", assay = "ANY",
    features = "ANY", layers = "ANY", scale_layer = "character",
    new_reduction = "character", reference = "ANY",
    normalization_method = "character", dims = "numeric",
    k_filter = "ANY", dims_to_integrate = "ANY",
    k_weight = "numeric", weight_reduction = "ANY",
    sd_weight = "numeric", sample_tree = "ANY",
    preserve_order = "logical", verbose = "logical",
    l2_norm = "logical", k_anchor = "numeric", k_score = "numeric",
    max_features = "numeric", nn_method = "character",
    n_trees = "numeric", eps = "numeric"
  )
)

#' @export
#' @rdname BatChefParams
setClass("FastMNNParams", contains = "BatChefParams", slots = c(extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("HarmonyParams", contains = "BatChefParams", slots = c(extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("ScanoramaParams",
  contains = "BatChefParams",
  slots = c(
    assay_type = "ANY",
    return_dimred = "logical", extra = "list"
  )
)

#' @export
#' @rdname BatChefParams
setClass("ScMerge2Params",
  contains = "BatChefParams",
  slots = c(assay_type = "character", extra = "list")
)

#' @export
#' @rdname BatChefParams
setClass("BBKNNParams",
  contains = "BatChefParams",
  slots = c(reduction = "character", extra = "list")
)

#' @export
#' @rdname BatChefParams
setClass("LigerParams",
  contains = "BatChefParams",
  slots = c(
    features = "character",
    method = "character", extra = "list"
  )
)
