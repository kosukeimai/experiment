#' @export
summary.completeRandomization <- function(object, ...) {
  cat("\nCall:\n", deparse(object$call), "\n\n", sep= "")
  cat("\nTreatment Assignment:\n")
  print(table(object$treatment, exclude = NULL, ...))
  
  cat("\nTotal number of observations:",
      length(na.omit(object$treatment)), "\n\n")
}