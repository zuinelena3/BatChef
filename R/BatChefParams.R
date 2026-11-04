#' BatChefParams methods
#'
#' Constructors and methods for the params parameter classes.
#' BatChefParams objects contain method specific parameters
#' to pass to the batchCorrect generic.
#'
#' @param assay_type A string specifying the assay.
#' @param ... Named arguments to pass to individual methods upon dispatch.
#'
#' @return
#' A BatChefParams object of the specified subclass, containing parameter
#' settings for the corresponding batch correction method.
#' @param assay_type A string specifying the assay.
#' @param return_dimred A logical to returning integrated low-dimesional
#' embeddings.
#' @param ... Named arguments to pass to individual methods upon dispatch
#'
#' @export
#' @rdname BatChefParams
#' @importFrom methods new
ScanoramaParams <- function(assay_type = NULL, return_dimred = FALSE, ...) {
  new("ScanoramaParams",
      assay_type = assay_type, return_dimred = return_dimred,
      extra = list(...)
  )
}

#' @export
#' @param reduction reduction
#' @rdname BatChefParams
#' @importFrom methods new
BBKNNParams <- function(reduction, ...) {
  new("BBKNNParams", reduction = reduction, extra = list(...))
}

