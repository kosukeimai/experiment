###
### Calculate the ATE without covariates
###
#' Estimation of the Average Treatment Effect in Randomized Experiments
#' 
#' This function computes the standard ``difference-in-means'' estimate of the
#' average treatment effect in randomized experiments without using
#' pre-treatment covariates. The treatment variable is assumed to be binary.
#' Currently, the two designs are allowed: complete randomized design and
#' matched-pair design.
#'
#' @useDynLib experiment 
#' @importFrom stats coef complete.cases cov fitted ftable lm mahalanobis model.frame model.matrix model.response na.fail na.omit printCoefmat qnorm quantile rbinom rnorm runif pnorm uniroot sd terms var vcov weighted.mean
#' @importFrom utils packageDescription
#' @importFrom MASS mvrnorm
#' @importFrom boot boot
#' 
#' @param Y The outcome variable of interest.
#' @param Z The (randomized) treatment variable. This variable should be
#' binary.
#' @param data A data frame containing the relevant variables.
#' @param match A variable indicating matched-pairs. The two units in the same
#' matched-pair should have the same value.
#' @return A list of class \code{ATEnocov} which contains the following items:
#' \item{call}{ The matched call.  } \item{Y}{ The outcome variable.  }
#' \item{Z}{ The treatment variable.  } \item{match}{ The matched-pair
#' indicator variable.  } \item{ATEest}{ The estimated average treatment
#' effect.  } \item{ATE.var}{ The estimated variance of the average treatment
#' effect estimator.  } \item{diff}{ Within-pair differences if the
#' matched-pair design is analyzed.  }
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu};
#' @references Imai, Kosuke, (2007). \dQuote{Randomization-based Inference and
#' Efficiency Analysis in Experiments under the Matched-Pair Design}, Technical
#' Report. Department of Politics, Princeton University.
#' @keywords design
#' @export ATEnocov
ATEnocov <- function(Y, Z, data = parent.frame(), match = NULL){

  ## an internal function that checks match and returns diff
  match.check <- function(Y, Z, match) { 
    n <- length(Y)
    if ((n %% 2) != 0)
      stop("pair randomization requires the even number of observations")
    if (length(unique(table(match))) > 1)
      stop("invalid input for `match'")
    if (unique(table(match)) != 2)
      stop("invalid input for `match'")
    umatch <- sort(unique(match))
    diff <- rep(NA, n/2)
    for (i in 1:length(umatch))
      diff[i] <- Y[(Z == 1) & (match == umatch[i])] -
        Y[(Z == 0) & (match == umatch[i])] 
    return(diff)
  }
  
  ## getting the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  match <- eval(call$match, envir = data)

  ## checking data
  if (sum(sort(unique(Z)) == c(0,1)) != 2)
    stop("`Z' should be binary taking the value of 0 or 1")
  if (length(Y) != length(Z))
    stop("`Y' and `Z' have different numbers of observations")
  if (!is.null(match))
    if (length(match) != length(Y))
      stop("`match' and `Y' have different numbers of observations")
    
  ## ATE for unit randomization
  res$ATE.est <- mean(Y[Z==1])-mean(Y[Z==0])
  res <- list(call = call, Y = Y, Z = Z, match = match)
  if (is.null(match)) { # without matching
    res$ATE.var <- var(Y[Z==1])/sum(Z==1)+var(Y[Z==0])/sum(Z==0)
  } else { # with matching
    res$diff <- diff <- match.check(Y, Z, match)
    res$ATE.var <- var(diff)/length(diff)
  }
  class(res) <- "ATEnocov"
  return(res)
}

