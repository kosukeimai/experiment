###
### Calculate the cluster adjusted variance for the sample mean
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

