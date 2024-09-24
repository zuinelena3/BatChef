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

#' @export
#' @rdname MethodParam
setClass("fastMNNMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("harmonyMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("scanoramaMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("bbknnMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("scVIMethod", contains = "MethodParam")

#' @export
#' @rdname MethodParam
setClass("scMergeMethod", contains = "MethodParam")
