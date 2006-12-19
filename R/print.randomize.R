print.randomize <- function(x, digits = getOption("digits"), ...) {

  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (!is.null(x$block))
    print(ftable(x$block, x$treatment, exclude = NULL, ...),
          digits = digits)
  else
    print(table(x$treatment, exclude = NULL, ...), digits = digits)

  cat("\nTotal number of observations:",
      length(na.rm(x$treatment)), "\n\n")
  invisible(x)
}
