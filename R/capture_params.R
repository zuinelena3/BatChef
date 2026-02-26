#' Capture parameters of the correction methods
#'
#' @param target Target parameters
#' @param params Parameters
#'
#' @importFrom purrr keep
#' @return Parameters
capture_params <- function(target, params) {
  ns <- target |>
    formals() |>
    names() |>
    keep(function(n) !startsWith(n, ".") && n %in% names(params))
  params[ns]
}
