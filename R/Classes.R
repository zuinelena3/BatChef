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
