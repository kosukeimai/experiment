#' Block Randomization of Treatment Assignment for Experiments
#' 
#' This function implements block randomization of treatment assignment for
#' randomized experiments. Complete randomization is conducted within blocks 
#' defined by a variable provided by the user.
#'  
#' @aliases blockRandomization
#' @param data A data frame for which the treatments are randomly assigned.
#' @param n Number of observations from each block to be randomized into the treatment group. 
#' Must be the same length as the number of unique values in \code{block}.
#' @param block A vector which denotes which block each observation is in. Observations will
#' be randomized within these blocks.
#' @param sort If \code{TRUE}, returns data in original order. If \code{FALSE}, returns data
#' grouped by the categories defined by block.
#' 
#' @return A list of class \code{blockRandomize} which contains the following items:
#' \item{call}{ the matched call. }
#' \item{data}{The data frame that was used to conduct the randomization.}
#' \item{treatment}{ The vector of randomized treatments.}
#' \item{block}{The block vector.}
#' 
#' @author Kosuke Imai, Professor of Government and Statistics, Harvard University
#' \email{imai@@harvard.edu}
#' @author Tyler Simko, PhD Student, Harvard University
#' \email{tsimko@@g.harvard.edu}
#' 
#' @examples 
#' ## create sample data
#' block <- rep(c("A", "B", "A", "C"), times = c(15,20,5,20))
#' data <- data.frame(x1 = rnorm(60), x2 = rnorm(60), block = block)
#' 
#' ## assigns 10, 5, and 5 to treatment in groups A, B, C
#' blockRandomization(data, n = c(10,5,5), block = data$block)
#' 
#' @export blockRandomization
#' 

blockRandomization <- function(data, n = NULL, block = NULL, sort = TRUE) {  
  
  if (is.null(n))     stop("n (number of observations to be assigned to treatment) must be provided.")
  if (is.null(block)) stop("block indicator must be provided.")
  
  if (length(n) != length(unique(block))) 
    stop("length of n must match number of unique categories in block")
  
  ## call
  call <- match.call()
  
  ## for ordering
  data$order <- 1:nrow(data)
  
  ## for requesting blocks
  names(n) <- unique(data$block)
  
  ## split by block categories
  splitUp <- split(data, data$block)
  
  ## apply complete randomization to each block
  mapplied <- mapply(function(x) {
      t <- completeRandomization(x, n[unique(as.character(x$block))])
    }, splitUp)
  
  dataToRet <- NULL
  dd <- apply(mapplied, 2, function(x) {
    dat <- x$data
    dat$treatment <- x$treatment
      rbind(dataToRet, dat)
    })
  
  (finalData <- do.call(rbind, dd))
  
  ## return data grouped by block if requested
  if (sort) {
    finalData <- finalData[order(finalData$order),]
  }
  
  res <- list(call = call, data = finalData, treatment = finalData$treatment,
              block = finalData$block)
  class(res) <- "blockRandomization"
  
  return(res)
}