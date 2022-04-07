#' Estimation of the Average Treatment Effects in Cluster-Randomized
#' Experiments
#' 
#' This function estimates various average treatment effect in
#' cluster-randomized experiments without using pre-treatment covariates. The
#' treatment variable is assumed to be binary. Currently, only the matched-pair
#' design is allowed. The details of the methods for this design are given in
#' Imai, King, and Nall (2007).
#' 
#' 
#' @param Y The outcome variable of interest.
#' @param Z The (randomized) cluster-level treatment variable. This variable
#' should be binary. Two units in the same cluster should have the same value.
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
#' @param fpc A logical variable indicating whether or not finite population
#' correction should be used for estimating the lower bound of CACE variance.
#' This is relevant only when \code{weights} are specified.
#' @return A list of class \code{ATEcluster} which contains the following
#' items: \item{call}{ The matched call.  } \item{n}{ The total number of
#' units.  } \item{n1}{ The total number of units in the treatment group.  }
#' \item{n0}{ The total number of units in the control group.  } \item{Y}{ The
#' outcome variable.  } \item{Y1bar}{ The cluster-specific (unweighted) average
#' value of the observed outcome for the treatment group.  } \item{Y0bar}{ The
#' cluster-specific (unweighted) average value of the observed outcome for the
#' treatment group.  } \item{Y1var}{ The cluster-specific sample variance of
#' the observed outcome for the treatment group.  } \item{Y0var}{ The
#' cluster-specific sample variance of the observed outcome for the control
#' group.  } \item{Z}{ The treatment variable.  } \item{grp}{ The
#' cluster-indicator variable.  } \item{match}{ The matched-pair indicator
#' variable.  } \item{weights}{ The weight variable in its original form.  }
#' \item{est}{ The estimated average treatment effect based on the arithmetic
#' mean weights.  } \item{var}{ The estimated variance of the average treatment
#' effect estimator based on the arithmetic mean weights. This uses the
#' variance formula provided in Imai, King, and Nall (2007).  } \item{var.lb}{
#' The estimated sharp lower bound of the cluster average treatment effect
#' estimator using the arithmetic mean weights.  } \item{est.dk}{ The estimated
#' average treatment effect based on the harmonic mean weights.  }
#' \item{var.dk}{ The estimated variance of the average treatment effect
#' estimator based on the harmonic mean weights. This uses the variance formula
#' provided in Donner and Klar (1993).  } \item{dkvar}{ The estimated variance
#' of the average treatment effect estimator based on the harmonic mean
#' weights. This uses the variance formula provided in Imai, King, and Nall
#' (2007).  } \item{eff}{ The estimated relative efficiency of the matched-pair
#' design over the completely randomized design (the ratio of two estimated
#' variances).  } \item{m}{ The number of pairs in the matched-pair design.  }
#' \item{N1}{ The population cluster sizes for the treatment group.  }
#' \item{N0}{ The population cluster sizes for the control group.  } \item{w1}{
#' Cluster-specific weights for the treatment group.  } \item{w0}{
#' Cluster-specific weights for the control group.  } \item{w}{ Pair-specific
#' normalized arithmetic mean weights. These weights sum up to the total number
#' of units in the sample, i.e., \code{n}.  } \item{w.dk}{ Pair-specific
#' normalized harmonic mean weights. These weights sum up to the total number
#' of units in the sample, i.e., \code{n}.  } \item{diff}{ Within-pair
#' differences if the matched-pair design is analyzed. This equals the
#' difference between \code{Y1bar} and \code{Y0bar}.  }
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' @references Donner, A. and N. Klar (1993). \dQuote{Confidence interval
#' construction for effect measures arising from cluster randomized trials.}
#' Journal of Clinical Epidemiology. Vol. 46, No. 2, pp. 123-131.
#' 
#' Imai, Kosuke, Gary King, and Clayton Nall (2007). \dQuote{The Essential Role
#' of Pair Matching in Cluster-Randomized Experiments, with Application to the
#' Mexican Universal Health Insurance Evaluation}, Technical Report. Department
#' of Politics, Princeton University.
#' @keywords design
#' @export ATEcluster
ATEcluster <- function(Y, Z, grp, data = parent.frame(),
                       match = NULL, weights = NULL, fpc = TRUE) {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  n <- length(Y)
  res <- list(call = call, n = n, Y = Y, Z = Z, grp = grp,
              match = match, weights = weights) 
  if (is.null(match))
    stop("This option is not yet available.")
  else {
    res$m <- m <- length(unique(match))
    res$Y1bar <- Y1bar <- tapply(Y[Z==1], match[Z==1], mean)
    res$Y0bar <- Y0bar <- tapply(Y[Z==0], match[Z==0], mean)
    res$diff <- diff <- Y1bar-Y0bar
    res$n1 <- n1 <- tapply(rep(1, sum(Z==1)), match[Z==1], sum)
    res$n0 <- n0 <- tapply(rep(1, sum(Z==0)), match[Z==0], sum)
  }

  if (is.null(weights)) {
    ## variance for PATE1 (sampling of clusters)
    N1 <- w1 <- n1
    N0 <- w0 <- n0
  } else {
    ## variance for PATE2 (double sampling)
    w1 <- N1 <- tapply(weights[Z==1], match[Z==1], mean)
    w0 <- N0 <- tapply(weights[Z==0], match[Z==0], mean)
  }
  w <- w1 + w0
  w <- n*w/sum(w)
  ## estimates
  ATE.est <- weighted.mean(diff, w)
  ATE.var <- m*sum((w*diff-n*ATE.est/m)^2)/((m-1)*(n^2))
  ## donner&klar methods:
  w.dk <- w1*w0/(w1 + w0)
  w.dk <- n*w.dk/sum(w.dk)
  ATEdk.est <- weighted.mean(diff, w.dk)
  ATEdk.var <- sum(w.dk^2)*sum(w.dk*(diff-ATEdk.est)^2)/(n^3)
  ATE.dkvar <- sum(w^2)*sum(w*(diff-ATE.est)^2)/(n^3)
  ## lower bound for CATE variance
  if (!is.null(weights)) {
    Y1var <- tapply(Y[Z==1], match[Z==1], var)/n1
    Y0var <- tapply(Y[Z==0], match[Z==0], var)/n0
    if (fpc) {
      Y1var <- (1-n1/N1)*Y1var
      Y0var <- (1-n0/N0)*Y0var
      if ((sum(n0 > N0)+sum(n1 > N1))>0)
        stop("population size is smaller than sample size")
    }
    res$Y1var <- Y1var
    res$Y0var <- Y0var
    res$var.lb <- sum((w/n)^2*(Y1var+Y0var))
  }
  ## unbiased estimation
  ##Y1sum <- tapply(Y[Z==1], match[Z==1], sum)
  ##Y0sum <- tapply(Y[Z==0], match[Z==0], sum)
  ##ATE.estU <- 2*(sum(Y1sum)-sum(Y0sum))/n
  ##ATE.varU <- 4*m*var(Y1sum-Y0sum)/(n^2)
  
  ## return the resutls
  res$est <- ATE.est
  res$est.dk <- ATEdk.est
  res$var <- ATE.var
  res$dkvar <- ATE.dkvar
  res$var.dk <- ATEdk.var
  res$eff <- 1/(1-2*cov(w*Y1bar, w*Y0bar)/(var(w*Y1bar)+var(w*Y0bar)))
  res$w <- w
  res$w.dk <- w.dk
  if (!is.null(match)) {
    res$w1 <- w1
    res$w0 <- w0
    res$N0 <- N0
    res$N1 <- N1
  }
  class(res) <- "ATEcluster"
  return(res)
}
