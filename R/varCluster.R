###
### Calculate the cluster adjusted variance for the sample mean based
### on the average method of Donner (Statistics in Medicine, 1992)
###
### the estimation of intraclass coefficient is based on the standard
### ANOVA formula. see Donner (International Statistical Review, 1986)
###
### if Z is specified, then it will use pooled rho: see chapter 7 of
### Donner and Klar book (page 114).
###

varCluster <- function(Y, Z = NULL, grp) {

  if (is.null(Z)) { # unpooled
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
  } else {
    ## prep
    M <- length(Y)
    Y1 <- Y[Z==1]
    Y0 <- Y[Z==0]
    grp1 <- grp[Z==1]
    grp0 <- grp[Z==0]
    m1 <- c(table(grp1))
    m0 <- c(table(grp0))
    ugrp1 <- unique(grp1)
    ugrp0 <- unique(grp0)
    k1 <- length(ugrp1)
    k0 <- length(ugrp0)
    mbar.A1 <- sum(m1^2)/length(Y1)
    mbar.A0 <- sum(m0^2)/length(Y0)
    Y1gbar <- rep(NA, k1)
    Y0gbar <- rep(NA, k0)
    ## within-group mean squares
    MSw <- 0
    for (i in 1:k0) {
      Y0gbar[i] <- mean(Y0[grp0 == ugrp0[i]])
      MSw <- MSw + sum((Y0[grp0 == ugrp0[i]]-Y0gbar[i])^2)
    }
    for (i in 1:k1) {
      Y1gbar[i] <- mean(Y1[grp1 == ugrp1[i]])
      MSw <- MSw + sum((Y1[grp1 == ugrp1[i]]-Y1gbar[i])^2)
    }
    MSw <- MSw / (M - k0 - k1)
    ## between-group mean squares
    MSb <- (sum(m0*(Y0gbar-mean(Y0))^2) +
            sum(m1*(Y1gbar-mean(Y1))^2))/(k0 + k1 - 2)
    ## ICC
    m0 <- (M - mbar.A1 - mbar.A0)/(k0 + k1 - 2)
    rho <- (MSb - MSw)/(MSb + (m0 - 1)*MSw)
    ## variance
    C0 <- 1 + (mbar.A0 - 1)*rho
    C1 <- 1 + (mbar.A1 - 1)*rho
    Sp2 <- MSw + (MSb - MSw)/m0
    res <- Sp2*(C0/length(Y0)+C1/length(Y1))
  }
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
