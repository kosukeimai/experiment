#' @export
summary.estimateComplete <- function(object, ...) {
  ## print call
  cat("\nCall:\n", deparse(object$call), "\n", sep= "")
  
  cat("\nEstimates:\n")
  print(object$estimates)
  
  cat("\n")
  cat(object$df, "Degrees of Freedom \n")
  cat("R-squared:", object$r.squared,"\n")
  cat("Adjusted R-squared", object$adj.r.squared,"\n")
}