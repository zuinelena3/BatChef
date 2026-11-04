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
setClass("ScanoramaParams",
         contains = "BatChefParams",
         slots = c(
           assay_type = "ANY",
           return_dimred = "logical", extra = "list"
           )
)

#' @export
#' @rdname BatChefParams
setClass("BBKNNParams",
         contains = "BatChefParams",
         slots = c(reduction = "character", extra = "list")
)
