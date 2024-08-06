#' AltOutput class
#'
#' Class of an alternative output, that contains three slots: corrected matrix, embedding, and metadata.
#'
#' @slot corrected ANY.
#' @slot embedding matrix.
#' @slot meta data.frame.
#'
#' @importFrom methods new
#' @export
#'
AltOutput <- setClass("AltOutput",
                      slots = c(
                        corrected = "ANY",
                        embedding = "ANY",
                        meta = "data.frame"
                      ))
