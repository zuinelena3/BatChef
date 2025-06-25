#' Parameters merging
#'
#' @param base base params
#' @param extra extra params
#' @param class class
#'
merge_params <- function(base, extra, class) {
  arg_base <- names(base)
  arg_extra <- names(extra)

  if (length(intersect(arg_base, arg_extra)) > 0) {
    print(paste(class, "argoments are duplicated!", sep = " "))
  }
  else return(c(base, extra))
}
