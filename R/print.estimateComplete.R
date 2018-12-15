#' @export
print.estimateComplete <- function(x, digits = getOption("digits"), ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\n Estimates:\n")
  print(x$estimates, digits = digits)
  invisible(x)
}
