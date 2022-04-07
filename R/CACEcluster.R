#' Estimation of the Complier Average Causal Effects in Cluster-Randomized
#' Experiments with Unit-level Noncompliance
#' 
#' This function estimates various complier average causal effect in
#' cluster-randomized experiments without using pre-treatment covariates when
#' unit-level noncompliance exists. Both the encouragement and treatment
#' variables are assumed to be binary. Currently, only the matched-pair design
#' is allowed. The details of the methods for this design are given in Imai,
#' King, and Nall (2007).
#' 
#' 
#' @param Y The outcome variable of interest.
#' @param D The unit-level treatment receipt variable. This variable should be
#' binary but can differ across units within each cluster.
#' @param Z The (randomized) cluster-level encouragement variable. This
#' variable should be binary. Two units in the same cluster should have the
#' same value.
#' @param grp A variable indicating clusters of units. Two units in the same
#' cluster should have the same value.
#' @param data A data frame containing the relevant variables.
#' @param match A variable indicating matched-pairs of clusters. Two units in
#' the same matched-pair of clusters should have the same value. The default is
#' \code{NULL} (i.e., no matching).
#' @param weights A variable indicating the population cluster sizes, which
#' will be used to construct weights for each pair of clusters. Two units in
#' the same cluster should have the same value. The default is \code{NULL}, in
#' which case sample cluster sizes will be used for constructing weights.
#' @param ...  Optional arguments passed to \code{ATEcluster}, which is called
#' internally.
#' @return A list of class \code{CACEcluster} which contains the following
#' items: \item{call}{ The matched call.  } \item{ITTY}{ The output object from
#' \code{ATEcluster} which is used to estimate the ITT effect of the
#' encouragement on the outcome variable.  } \item{ITTD}{ The output object
#' from \code{ATEcluster} which is used to estimate the ITT effect of the
#' encouragement on the treatment receipt variable.  } \item{n1}{ The total
#' number of units in the treatment group.  } \item{n0}{ The total number of
#' units in the control group.  } \item{Z}{ The treatment variable.  }
#' \item{est}{ The estimated complier average causal effect.  } \item{var}{ The
#' estimated variance of the complier average causal effect estimator.  }
#' \item{cov}{ The estimated covariance between two ITT estimator.  } \item{m}{
#' The number of pairs in the matched-pair design.  } \item{N1}{ The population
#' cluster sizes for the treatment group.  } \item{N0}{ The population cluster
#' sizes for the control group.  } \item{w}{ Pair-specific normalized
#' arithmetic mean weights. These weights sum up to the total number of units
#' in the sample, i.e., \code{n}.  }
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' @references Imai, Kosuke, Gary King, and Clayton Nall (2007). \dQuote{The
#' Essential Role of Pair Matching in Cluster-Randomized Experiments, with
#' Application to the Mexican Universal Health Insurance Evaluation}, Technical
#' Report. Department of Politics, Princeton University.
#' @keywords design
#' @export CACEcluster
CACEcluster <- function(Y, D, Z, grp, data = parent.frame(),
                        match = NULL, weights = NULL, ...) {

  cov.internal <- function(x){
    return(cov(x[,1],x[,2]))
  }
  
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)
  
  ITTY <- ATEcluster(Y = Y, Z = Z, grp = grp, match = match,
                     weights = weights, ...)
  ITTD <- ATEcluster(Y = D, Z = Z, grp = grp, match = match,
                     weights = weights, ...)

  ## point estimate
  n <- ITTY$n
  ITTY.est <- ITTY$est
  ITTD.est <- ITTD$est
  CACE.est <- ITTY.est/ITTD.est

  ## outputs
  res <- list(est = CACE.est, ITTY = ITTY, ITTD = ITTD)

  ## calculation
  if (is.null(match))
    stop("the estimator is not yet available.")
  else {
    res$n1 <- n1 <- ITTY$n1
    res$n0 <- n0 <- ITTY$n0
    res$N1 <- N1 <- ITTY$N1
    res$N0 <- N0 <- ITTY$N0
    res$m <- m <- ITTY$m
    res$w <- w <- ITTY$w
    diffY <- ITTY$diff
    diffD <- ITTD$diff
    res$cov <- Cov <- m*sum((diffY*w - n*ITTY.est/m) *
                            (diffD*w - n*ITTD.est/m))/((m-1)*(n^2))
  }
  
  res$var <- CACE.var <-
    (ITTY$var*(ITTD.est^2) + ITTD$var*(ITTY.est^2) -
     2*Cov*ITTY.est*ITTD.est)/(ITTD.est^4)
  
  class(res) <- "CACEcluster"
  return(res)
}
