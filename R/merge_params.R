#' Parameters merging
#'
#' @param base base params
#' @param extra extra params
#' @param class class
#'
merge_params <- function(base, extra, class) {
  arg_base <- names(base)
  arg_extra <- names(extra)

  if (any(arg_base %in% arg_extra)) {
    message(paste(class, "arguments are duplicated!"))
  }
  else return(c(base, extra))
}



