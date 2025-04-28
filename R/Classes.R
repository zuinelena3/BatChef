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
setClass("LimmaParams", contains = "BatChefParams", slots = c(assay_type = "character", extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("CombatParams", contains = "BatChefParams", slots = c(assay_type = "character", extra = "list"))

#' @export
#' @rdname BatChefParams
setClass("SeuratV3Params", contains = "BatChefParams", slots = c(features = "character", pca_name = "ANY", assay = "ANY", reference = "ANY",
                                                                 anchor.features = "numeric", scale = "logical",
                                                                 normalization.method = "character",
                                                                 sct.clip.range = "ANY", reduction = "character",
                                                                 l2.norm = "logical", dims = "numeric", k.anchor = "numeric",
                                                                 k.filter = "numeric", k.score = "numeric",
                                                                 max.features = "numeric", nn.method = "character",
                                                                 n.trees = "numeric", eps = "numeric", verbose = "logical",
                                                                 new.assay.name = "character",
                                                                 features.to.integrate = "ANY", k.weight = "numeric",
                                                                 weight.reduction = "ANY", sd.weight = "numeric", sample.tree = "ANY",
                                                                 preserve.order = "logical"))

#' @export
#' @rdname BatChefParams
setClass("SeuratV5Params", contains = "BatChefParams", slots = c(pca_name = "ANY", method = "character", orig.reduction = "character",
                                                                 assay = "ANY", features = "ANY", layers = "ANY", scale.layer = "character",
                                                                 new.reduction = "character", reference = "ANY",
                                                                 normalization.method = "character", dims = "numeric",
                                                                 k.filter = "ANY", dims.to.integrate = "ANY", k.weight = "numeric",
                                                                 weight.reduction = "ANY", sd.weight = "numeric", sample.tree = "ANY",
                                                                 preserve.order = "logical", verbose = "logical",
                                                                 l2.norm = "logical", k.anchor = "numeric", k.score = "numeric",
                                                                 max.features = "numeric", nn.method = "character",
                                                                 n.trees = "numeric", eps = "numeric"))
