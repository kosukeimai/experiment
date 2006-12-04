###
### Calculate the cluster adjusted variance for the sample mean based
### on the average method of Donner (Statistics in Medicine, 1992)
###
### the estimation of intraclass coefficient is based on the standard
### ANOVA formula. see Donner (International Statistical Review, 1986)
###
###

varCluster <- function(Y, grp) {

  ## number of obs within each group
  n <- c(table(grp))
  ## groups
  ugrp <- unique(grp)
  ## number of groups
  k <- length(ugrp)
  ## total number of obs
  N <- length(Y)

  ## total mean
  Ybar <- mean(Y)
  ## group mean
  Ygbar <- rep(NA, k)
  ## within-group mean squares
  MSw <- 0
  for (i in 1:k) {
    Ygbar[i] <- mean(Y[grp == ugrp[i]])
    MSw <- MSw + sum((Y[grp == ugrp[i]]-Ygbar[i])^2)
  }
  MSw <- MSw / (N-k)
  ## between-group mean squares
  MSb <- sum(n*(Ygbar-Ybar)^2)/(k-1)
  ## intraclass correlation coefficient estimate based on ANOVA
  n0 <- mean(n)-sum((n-mean(n))^2)/((k-1)*N)
  rho <- (MSb - MSw)/(MSb + (n0-1)*MSw)
  ## cluster-adjusted variance based on the average method
  res <- var(Y)*(1+(mean(n)-1)*rho)/N
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
  n0 <- mean(n)-sum((n-mean(n))^2)/((k-1)*N)
  rho <- (MSb - MSw)/(MSb + (n0-1)*MSw)
  ## cluster-adjusted covariance
  res <- var(Y1, Y2)*(1+(mean(n)-1)*rho)/N
  return(list(cov = res, MSb = MSb, MSw = MSw, rho = rho))
}
