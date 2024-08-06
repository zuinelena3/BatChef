#' MethodParam class
#'
#' @export
#' @importClassesFrom S4Vectors SimpleList
#' @import methods
#'
#' @rdname MethodParam
setClass("MethodParam", contains = "SimpleList")

#' @export
#' @rdname MethodParam
setClass("limmaMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("combatMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("seuratv3Method", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("seuratv5Method", contains = "MethodParam")
