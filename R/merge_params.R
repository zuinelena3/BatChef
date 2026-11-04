#' Parameters merging
#'
#' @param base base params
#' @param extra extra params
#' @param class class
#' @return Vector of strings of base and extra parameters.
#'
merge_params <- function(base, extra, class) {
  arg_base <- names(base)
  arg_extra <- names(extra)

  class_msg <- paste(class, "arguments are duplicated!")
  if (any(arg_base %in% arg_extra)) {
    message(class_msg)
  } else {
    return(c(base, extra))
  }
}
