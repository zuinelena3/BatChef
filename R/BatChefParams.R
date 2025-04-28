#' BatChefParams methods
#'
#' Constructors and methods for the params parameter classes.
#' BatChefParams objects contain method specific parameters to pass to the batchCorrect generic.
#'
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
LimmaParams <- function(assay_type = "logcounts", ...) {
  new("LimmaParams", assay_type = assay_type, extra = list(...))
}

#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
CombatParams <- function(assay_type = "counts", ...) {
  new("CombatParams", assay_type = assay_type, extra = list(...))
}
