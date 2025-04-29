#' Title
#'
#' @param target target
#' @param params params
#'
#' @importFrom purrr keep
capture_params <- function(target, params) {
  ns <- target |> formals() |> names() |> keep(function(n) !startsWith(n, ".") && n %in% names(params))
  params[ns]
}
