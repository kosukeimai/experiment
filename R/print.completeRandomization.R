#' @export
print.completeRandomization <- function(x, digits = getOption("digits"), ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("\n Treatment Assignment:\n")
  print(table(x$treatment, exclude = NULL, ...), digits = digits)
  
  cat("\nTotal number of observations:",
      length(na.omit(x$treatment)), "\n\n")
  invisible(x)
}
