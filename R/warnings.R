suppressWarningsByMsg <- function(msg, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      matches <- sapply(msg, \(m) grepl(m, conditionMessage(w)))
      if (any(matches)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
