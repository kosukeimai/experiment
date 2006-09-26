###
### Calculate the cluster adjusted variance for the sample mean
### the calculation is based on the standard ANOVA formula
###

varCluster <- function(Y, grp) {

  n <- c(table(grp))
  ugrp <- unique(grp)
  k <- length(ugrp)
  N <- length(Y)
  Ybar <- mean(Y)
  Ygbar <- rep(NA, k)
  MSw <- 0
  ## within-group mean squares
  for (i in 1:k) {
    Ygbar[i] <- mean(Y[grp == ugrp[i]])
    MSw <- MSw + sum((Y[grp == ugrp[i]]-Ygbar[i])^2)
  }
  MSw <- MSw / (N-k)
  ## between-group mean squares
  MSb <- sum(n*(Ygbar-Ybar)^2)/(k-1)
  ## intraclass correlation coefficient
  n0 <- (N - sum(n^2)/N)/(k-1)
  rho <- (MSb - MSw)/(MSb + (n0-1)*MSw)
  ## cluster-adjusted variance
  res <- var(Y)*(1+(k-1)*rho)/N
  return(list(var = res, MSb = MSb, MSw = MSw, rho = rho))
}

###
### Calculate the cluster-adjusted covariance for two sample means 
### the calculation is based on the formula analogous to the one
### above.
###

covCluster <- function(Y1, Y2, grp) {
  
  n <- c(table(grp))
  ugrp <- unique(grp)
  k <- length(ugrp)
  N <- length(Y1)
  Y1bar <- mean(Y1)
  Y2bar <- mean(Y2)
  Y1gbar <- Y2gbar <- rep(NA, k)
  MSw <- 0
  ## within-group mean squares
  for (i in 1:k) {
    Y1gbar[i] <- mean(Y1[grp == ugrp[i]])
    Y2gbar[i] <- mean(Y2[grp == ugrp[i]])
    MSw <- MSw + sum((Y1[grp == ugrp[i]]-Y1gbar[i])*(Y2[grp == ugrp[i]]-Y2gbar[i]))
  }
  MSw <- MSw / (N-k)
  ## between-group mean squares
  MSb <- sum(n*(Y1gbar-Y1bar)*(Y2gbar-Y2bar))/(k-1)
  ## intraclass correlation coefficient
  n0 <- (N - sum(n^2)/N)/(k-1)
  rho <- (MSb - MSw)/(MSb + (n0-1)*MSw)
  ## cluster-adjusted covariance
  res <- var(Y1, Y2)*(1+(k-1)*rho)/N
  return(list(cov = res, MSb = MSb, MSw = MSw, rho = rho))
}
