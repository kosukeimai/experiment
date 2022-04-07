###
### Calculates the bounds for the ATE in the presence of missing
### response 
###



#' Bounding the Average Treatment Effect when some of the Outcome Data are
#' Missing
#' 
#' This function computes the sharp bounds on the average treatment effect when
#' some of the outcome data are missing. The confidence intervals for the
#' bounds are also computed.
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param formula A formula of the form \code{Y ~ X} where \code{Y} is the name
#' of the outcome variable and \code{X} is the name of the (randomized)
#' treatment variable. \code{X} should be a factor variable but its value can
#' take more than two levels. The missing values for \code{Y} should be coded
#' as \code{NA}.
#' @param data A data frame containing the relevant variables.
#' @param maxY A scalar. The maximum value of the outcome variable. The default
#' is the maximum sample value.
#' @param minY A scalar. The minimum value of the outcome variable. The default
#' is the minimum sample value.
#' @param alpha A positive scalar that is less than or equal to 0.5. This will
#' determine the (1-\code{alpha}) level of confidence intervals. The default is
#' \code{0.05}.
#' @param strata The variable name indicating strata. If this is specified, the
#' quantities of interest will be first calculated within each strata and then
#' aggregated. The default is \code{NULL}.
#' @param ratio A \eqn{J \times M} matrix of probabilities where \eqn{J} is the
#' number of strata and \eqn{M} is the number of treatment and control groups.
#' Each element of the matrix specifies the probability of a unit falling into
#' that category. The default is \code{NULL} in which case the sample estimates
#' of these probabilities are used for computation.
#' @param survey The variable name for survey weights. The default is
#' \code{NULL}.
#' @param n.reps A positive integer. The number of bootstrap replicates used
#' for the construction of confidence intervals via B-method of Berran (1988).
#' If it equals zero, the confidence intervals will not be constructed.
#' @param ...  The arguments passed to other functions.
#' @return A list of class \code{ATEbounds} which contains the following items:
#' \item{call}{ The matched call.  } \item{Y}{ The outcome variable.  }
#' \item{D}{ The treatment variable.  } \item{bounds}{ The point estimates of
#' the sharp bounds on the average treatment effect.  } \item{bounds.Y}{ The
#' point estimates of the sharp bounds on the outcome variable within each
#' treatment/control group.  } \item{bmethod.ci}{ The B-method confidence
#' interval of the bounds on the average treatment effect.  } \item{bonf.ci}{
#' The Bonferroni confidence interval of the bounds on the average treatment
#' effect.  } \item{bonf.ci.Y}{ The Bonferroni confidence interval of the
#' bounds on the outcome variable within each treatment/control group.  }
#' \item{bmethod.ci.Y}{ The B-method confidence interval of the bounds on the
#' outcome variable within each treatment/control group.  } \item{maxY}{ The
#' maximum value of the outcome variable used in the computation.  }
#' \item{minY}{ The minimum value of the outcome variable used in the
#' computation.  } \item{nobs}{ The number of observations.  } \item{nobs.Y}{
#' The number of observations within each treatment/control group.  }
#' \item{ratio}{ The probability of treatment assignment (within each strata if
#' \code{strata} is specified) used in the computation.  }
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' @references Horowitz, Joel L. and Charles F. Manski. (1998).
#' \dQuote{Censoring of Outcomes and Regressors due to Survey Nonresponse:
#' Identification and Estimation Using Weights and Imputations.} \emph{Journal
#' of Econometrics}, Vol. 84, pp.37-58.
#' 
#' Horowitz, Joel L. and Charles F. Manski. (2000).  \dQuote{Nonparametric
#' Analysis of Randomized Experiments With Missing Covariate and Outcome Data.}
#' \emph{Journal of the Americal Statistical Association}, Vol. 95, No. 449,
#' pp.77-84.
#' 
#' Harris-Lacewell, Melissa, Kosuke Imai, and Teppei Yamamoto. (2007).
#' \dQuote{Racial Gaps in the Responses to Hurricane Katrina: An Experimental
#' Study}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords design
#' @export ATEbounds
ATEbounds <- function(formula, data = parent.frame(), maxY = NULL,
                      minY = NULL, alpha = 0.05, n.reps = 0,
                      strata = NULL, ratio = NULL, survey = NULL, ...) {

  ## getting Y and D
  call <- match.call()
  tm <- terms(formula)
  attr(tm, "intercept") <- 0
  mf <- model.frame(tm, data = data, na.action = 'na.pass')
  D <- model.matrix(tm, data = mf)
  M <- ncol(D)
  if (max(D) > 1 || min(D) < 0)
    stop("the treatment variable should be a factor variable.")
  Y <- model.response(mf)
  if (is.null(maxY))
    maxY <- max(Y, na.rm = TRUE)
  if (is.null(minY))
    minY <- min(Y, na.rm = TRUE)
  if (!is.null(call$survey))
    survey <- eval(call$survey, data)
  else
    survey <- rep(1, length(Y))
  ### computing the bounds
  if (!is.null(call$strata)) {
    strata <- eval(call$strata, data)
    res <- boundsAggComp(cbind(Y, strata, D), rep(1, length(Y)), maxY,
                         minY, alpha = alpha, ratio = ratio, survey = survey)
  } else {
    res <- boundsComp(cbind(Y, D), rep(1, length(Y)), maxY, minY,
                      alpha = alpha, survey = survey)
  }
  
  ## CI based on B-method
  if (n.reps > 0) {
    if (!is.null(call$strata)) {
      breps <- boot(data = cbind(Y, strata, D), statistic = boundsAggComp, 
                    R = n.reps, maxY = maxY, minY = minY, alpha =
                    NULL, survey = survey)$t
      res$bmethod.ci <- res$bonf.ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
      counter <- 1
      for (i in 1:(M-1)) 
        for (j in (i+1):M) {
          tmp <- boundsCI(breps[,counter], breps[,counter+1],
                          res$bounds[(counter+1)/2,1],
                          res$bounds[(counter+1)/2,2], alpha)
          res$bmethod.ci[(counter+1)/2,] <- tmp$bmethod
          res$bonf.ci[(counter+1)/2,] <- tmp$bonferroni
          counter <- counter + 2
        }
    } else {
      breps <- boot(data = cbind(Y, D), statistic = boundsComp,
                    R = n.reps, maxY = maxY, minY = minY, alpha = NULL,
                    survey = survey)$t
      res$bmethod.ci.Y <- matrix(NA, ncol = 2, nrow = M)
      res$bmethod.ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
      for (i in 1:M) { 
        tmp <- boundsCI(breps[,(i-1)*2+1], breps[,i*2],
                        res$bounds.Y[i,1],
                        res$bounds.Y[i,2], alpha)
        res$bmethod.ci.Y[i,] <- tmp$bmethod
        res$bonf.ci.Y[i,] <- tmp$bonferroni
      }
      counter <- 1
      for (i in 1:(M-1)) 
        for (j in (i+1):M) {
          tmp <- boundsCI(breps[,2*M+counter], breps[,2*M+counter+1],
                          res$bounds[(counter+1)/2,1],
                          res$bounds[(counter+1)/2,2], alpha)
          res$bmethod.ci[(counter+1)/2,] <- tmp$bmethod
          res$bonf.ci[(counter+1)/2,] <- tmp$bonferroni
          counter <- counter + 2
        }
    }
  }
  
  ## dimnames
  tmp <- NULL
  for (i in 1:(M-1)) 
    for (j in (i+1):M) 
      tmp <- c(tmp, paste(colnames(D)[i], "-", colnames(D)[j]))
  if (is.null(call$strata)) {
    rownames(res$bounds.Y) <- rownames(res$bonf.ci.Y) <- colnames(D)
    rownames(res$bounds) <- rownames(res$bonf.ci) <- tmp
    colnames(res$bounds) <- colnames(res$bounds.Y) <- c("lower", "upper")
    colnames(res$bonf.ci) <- colnames(res$bonf.ci.Y) <-
      c(paste("lower ", alpha/2, "%CI", sep=""),
        paste("upper ", 1-alpha/2, "%CI", sep=""))
    if (n.reps > 0) {
      rownames(res$bmethod.ci.Y) <- colnames(D)
      rownames(res$bmethod.ci) <- tmp
      colnames(res$bmethod.ci) <- colnames(res$bmethod.ci.Y) <-
        c(paste("lower ", alpha/2, "%CI", sep=""),
          paste("upper ", 1-alpha/2, "%CI", sep=""))
    }
  } else {
    rownames(res$bounds) <- tmp
    colnames(res$bounds) <- c("lower", "upper")
    if (n.reps > 0) {
      rownames(res$bmethod.ci) <- rownames(res$bonf.ci) <- tmp
      colnames(res$bmethod.ci) <- colnames(res$bonf.ci) <-
        c(paste("lower ", alpha/2, "%CI", sep=""),
          paste("upper ", 1-alpha/2, "%CI", sep=""))
    }
  }
  res$Y <- Y
  res$D <- D
  res$call <- call
  class(res) <- "ATEbounds"
  return(res)
}

###
### An internal function which computes the bounds and bonferroni CI
### if alpha is specified (when alpha = NULL, then it returns a vector
### of bounds; this is used for bootstrap)
###

boundsComp <- function(data, weights, maxY, minY, alpha = NULL,
                       survey = NULL) {
  Y <- data[,1]
  D <- data[,-1]
  M <- ncol(D)
  bounds.Y <- ci.Y <- vars.Y <- matrix(NA, ncol = 2, nrow = M)
  nobs.Y <- NULL
  if (is.null(survey))
    survey <- rep(1, length(Y))
  for (i in 1:M) {
    Ysub <- Y[D[,i]==1]
    w <- weights[D[,i]==1]*survey[D[,i]==1]
    n <- length(Ysub)
    Ymax <- Ymin <- Ysub
    Ymax[is.na(Ysub)] <- maxY
    Ymin[is.na(Ysub)] <- minY
    ## point estimates of the bounds
    bounds.Y[i,] <- c(weighted.mean(Ymin, w), weighted.mean(Ymax, w))
    if (!is.null(alpha)) {
      ## variances
      vars.Y[i,] <- c(weighted.var(Ymin, w)*sum(w^2)/(sum(w)^2),
                      weighted.var(Ymax, w)*sum(w^2)/(sum(w)^2))
      ## Bonferroni bounds
      ci.Y[i,] <- c(bounds.Y[i,1] - qnorm(1-alpha/2)*sqrt(vars.Y[i,1]),
                    bounds.Y[i,2] + qnorm(1-alpha/2)*sqrt(vars.Y[i,2]))
    }
    nobs.Y <- c(nobs.Y, n)
  }
  
  ## Bounds for the ATE
  bounds <- ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
  counter <- 1
  nobs <- tmp <- NULL
  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      bounds[counter,] <- c(bounds.Y[i,1]-bounds.Y[j,2],
                            bounds.Y[i,2]-bounds.Y[j,1])
      if (!is.null(alpha))
        ci[counter,] <- c(bounds[counter,1] -
                          qnorm(1-alpha/2)*sqrt(vars.Y[i,1]+vars.Y[j,2]),
                          bounds[counter,2] +
                          qnorm(1-alpha/2)*sqrt(vars.Y[i,2]+vars.Y[j,1]))
      counter <- counter + 1
      nobs <- c(nobs, nobs.Y[i]+nobs.Y[j])
    }
  }
  
  if (is.null(alpha))
    return(c(t(rbind(bounds.Y, bounds))))
  else
    return(list(bounds.Y = bounds.Y, bounds = bounds, bonf.ci = ci,
                bonf.ci.Y = ci.Y, maxY = maxY, minY = minY,
                nobs = nobs, nobs.Y = nobs.Y))
}

###
### Aggregate bounds
###

boundsAggComp <- function(data, weights, maxY, minY, alpha = NULL,
                          ratio = NULL, survey = NULL) {
  Y <- data[,1]
  S <- data[,2]
  Svalue <- unique(S)
  J <- length(Svalue)
  D <- data[,3:ncol(data)]
  M <- ncol(D)
  ## compute bounds within each strata and weights across strata
  res.sub <- list()
  if (is.null(ratio))
    ratio.cal <- TRUE
  else
    ratio.cal <- FALSE
  if (ratio.cal)
    ratio <- matrix(NA, nrow = J, ncol = M)
  if (is.null(survey))
    survey <- rep(1, length(Y))
  for (i in 1:J) {
    sub <- (S == Svalue[i])
    res.sub[[i]] <- boundsComp(data[sub,-2], weights[sub], maxY, minY,
                               0.05, survey[sub])
    if (ratio.cal)
      for (j in 1:M)
        ratio[i,j] <- sum(weights[sub & (D[,j] == 1)])
  }
  if (ratio.cal)
    ratio <- ratio/sum(weights)
  omega <- matrix(NA, nrow = (M-1)*M/2, ncol = J)
  counter <- 1
  for (j in 1:(M-1)) {
    for (k in (j+1):M) {
      tmp <- 0
      for (i in 1:J) {
        omega[counter,i] <- ratio[i,j] + ratio[i,k]
        tmp <- tmp + omega[counter,i]
      }
      omega[counter,] <- omega[counter,]/tmp
      counter <- counter + 1
    }
  }
  
  ## aggregate the results
  bounds <- matrix(0, ncol = 2, nrow = choose(M, 2))
  counter <- 1
  for (j in 1:(M-1)) {
    for (k in (j+1):M) {
      for (i in 1:J) {
        bounds[counter,] <- bounds[counter,] +
          (res.sub[[i]]$bounds)[counter,]*omega[counter,i]
      }
      counter <- counter + 1
    }
  }
  if (is.null(alpha))
    return(c(t(bounds)))
  else
    return(list(bounds = bounds, maxY = maxY, minY = minY,
                ratio = ratio))
    #return(list(bounds = bounds, maxY = maxY, minY = minY,
    #            ratio = ratio, omega = omega))
}
