#' MethodParam methods
#'
#' Constructors and methods for the method parameter classes.
#'
#' @param ...  Named arguments to pass to individual methods upon dispatch
#'
#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
limmaMethod <- function(...) {
  new("limmaMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
combatMethod <- function(...) {
  new("combatMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
seuratv3Method <- function(...) {
  new("seuratv3Method", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
seuratv5Method <- function(...) {
  new("seuratv5Method", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
fastMNNMethod <- function(...) {
  new("fastMNNMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
harmonyMethod <- function(...) {
  new("harmonyMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
scanoramaMethod <- function(...) {
  new("scanoramaMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
bbknnMethod <- function(...) {
  new("bbknnMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
scVIMethod <- function(...) {
  new("scVIMethod", SimpleList(list(...)))
}

#' @rdname MethodParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#'
scMergeMethod <- function(...) {
  new("scMergeMethod", SimpleList(list(...)))
}
