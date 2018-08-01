#' Complete Randomization of the Treatment Assignment for Conducting Experiments
#' 
#' This function implements the complete randomization of treatment assignment for
#' randomized experiments.
#'  
#' @aliases completeRandomization
#' @param data A data frame for which the treatments are randomly assigned.
#' @param n The number of observations from \code{data} to be randomly assigned 
#' to each treatment group. Defaults to \code{nrow(data)/2} for 2 treatment groups, 
#' and rounds down if \code{nrow(data)} is odd. If a single number is provided, \code{n} 
#' provides the number of units in the treatment group. A vector can be used for two groups 
#' as well (ex: \code{c(10,10)}). If more than two groups are necessary, provide a vector with the 
#' number of observations desired for each group (ex: \code{c(10,5,5)}). The total number 
#' of treatment options is taken to be the length of this vector unless a single number for \code{n}
#' is provided.
#' @return A list of class \code{completeRandomize} which contains the following items:
#' \item{call}{ the matched call. } \item{treatment}{ The vector of randomized
#' treatments.} \item{data}{The data frame that was used to conduct the
#' randomization.}
#' @author Kosuke Imai, Professor of Government and Statistics, Harvard University
#' \email{imai@@harvard.edu};
#' @author Tyler Simko, PhD Student, Harvard University
#' \email{tsimko@@g.harvard.edu};
#' 
#' @examples 
#' ## example data
#' data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' 
#' ## assigns half of the units to treatment, half to control
#' completeRandomization(data) 
#' 
#' ## assigns 10 units to treatment, remainder to control
#' completeRandomization(data, n = 10) 
#' 
#' ## 10 units to treatment, 10 to control
#' completeRandomization(data, n = c(10, 10))
#' 
#' ## assigns 10 units to first group, 5 to other two groups.
#' completeRandomization(data, n = c(10, 5, 5))
#' @export completeRandomization
#' 

completeRandomization <- function(data, n = nrow(data)/2) {  
  
  ## call
  call <- match.call()
  
  ## single treatment group, no vector provided
  if (length(n) == 1) {
    assign <- ifelse(1:nrow(data) %in% sample(1:nrow(data), n), 1, 0)
  }
  
  ## single treatment group, vector provided
  else if (length(n) == 2) {
    assign <- ifelse(1:nrow(data) %in% sample(1:nrow(data), n), 1, 0)
  }
  
  ## multiple treatment groups, vector provided
  else if (length(n) > 2){
    assign <- sample(rep(1:length(n) - 1, n))
  }
  
  if (sum(n) > nrow(data)) warning("More requested assignments than number of observations.")
  
  res <- list(call = call, data = data, treatment = assign)
  class(res) <- "completeRandomization"
  
  return(res)
}