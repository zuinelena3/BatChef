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
