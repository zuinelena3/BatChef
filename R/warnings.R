suppressWarningsByMsg <- function(msg, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      matches <- vapply(msg, \(m) grepl(m, conditionMessage(w)), FUN.VALUE = logical(1))
      if (any(matches)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
