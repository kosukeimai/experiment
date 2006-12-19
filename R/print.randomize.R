print.randomize <- function(x, digits = getOption("digits"), ...) {

  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (!is.null(block))
    print(ftable(object$block, object$treatment, exclude = NULL, ...),
          digits = digits)
  else
    print(table(object$treatment, exclude = NULL, ...), digits = digits)

  cat("\nTotal number of observations:",
      length(na.rm(object$treatment)), "\n\n")
  invisible(x)
}
