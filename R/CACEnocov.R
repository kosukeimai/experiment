###
### Calculate the CACE for individual and group randomized trials 
###

CACEnocov <- function(Y, D, Z, data = parent.frame(), grp = NULL,
                      match = NULL, grp.size = NULL, grp.method = "neyman"){

  ## get the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  grp.size <- eval(call$grp.size, envir = data)
  match <- eval(call$match, envir = data)
  
  ## point estimates and variances
  ITTY <- ATEnocov(Y = Y, Z = Z, grp = grp, match = match,
                   grp.size = grp.size, grp.method = grp.method) 
  ITTD <- ATEnocov(Y = D, Z = Z, grp = grp, match = match,
                   grp.size = grp.size, grp.method = grp.method) 
  CACEest <- ITTY$est/ITTD$est

  ## covariance calculation
  if (is.null(grp) && is.null(grp.size)) { # unit randomization
    if (is.null(match)) { # without matching
      Cov <- cov(Y[Z==0], D[Z==0])/sum(Z==0)+
        cov(Y[Z==1], D[Z==1])/sum(Z==1)
      CACEvar <- (ITTY$var*(ITTD$est^2) + ITTD$var*(ITTY$est^2) -
                  2*Cov*ITTY$est*ITTD$est)/(ITTD$est^4)
    } else { # with matching
      K <- length(ITTY$diff)
      Cov <- cov(ITTY$diff, ITTD$diff)/K
      CACEvar <- (var(ITTY$diff)*(ITTD$est^2)/K +
                  var(ITTD$diff)*(ITTY$est^2)/K -
                  2*Cov**ITTY$est*ITTD$est)/(ITTD$est^4)
    }
  } else if (grp.method %in% c("textbook", "unpooled")) {
    if (is.null(match)) { # without matching
      if (grp.method == "unpooled") {
        Cov <- covCluster(Y[Z==0], D[Z==0], grp[Z==0])$cov +
          covCluster(Y[Z==1], D[Z==1], grp[Z==1])$cov
      } else {
        Cov <- covCluster(Y, D, grp)$cov
      }
      CACEvar <- (ITTY$var*(ITTD$est^2) + ITTD$var*(ITTY$est^2) -
                  2*Cov*ITTY$est*ITTD$est)/(ITTD$est^4)
    } else { # with matching
      stop("this estimator is not available for matched-pair designs.")
    }
  } else if (grp.method == "standard") {
    if (is.null(match)) { # without matching
      Cov <- cov(ITTY$Ysum[ITTY$Z==0],
                 ITTD$Ysum[ITTD$Z==0])/sum(ITTY$Z==0) +
                   cov(ITTY$Ysum[ITTY$Z==1],
                       ITTD$Ysum[ITTD$Z==1])/sum(ITTY$Z==1)
      CACEvar <- (ITTY$var*(ITTD$est^2) + ITTD$var*(ITTY$est^2) -
                  2*Cov*ITTY$est*ITTD$est)/(ITTD$est^4)
    } else { # with matching
      stop("this estimator is not available for matched-pair designs.")
    }
  } else { # my method
    if (is.null(match)) { # without matching
      if (ITTD$M == (ITTD$m1*2)) { 
        Cov <- cov(ITTY$Ysum[ITTY$Z==0],
                   ITTD$Ysum[ITTD$Z==0])/sum(ITTY$Z==0) 
        + cov(ITTY$Ysum[ITTY$Z==1],
              ITTD$Ysum[ITTD$Z==1])/sum(ITTY$Z==1)
        CACEvar <- 2*ITTY$M*((var(ITTY$Ysum[ITTY$Z==1]) +
                              var(ITTY$Ysum[ITTY$Z==1]))*ITTD$est^2 +
                             (var(ITTD$Ysum[ITTD$Z==1]) +
                              var(ITTD$Ysum[ITTD$Z==1]))*ITTY$est^2 -
                             2*cov*ITTY$est*ITTD$est)/((ITTY$N^2)*(ITTD$est^4))
      } else {
        stop("the treatment and control groups must have the same number of clusters for this estimator")
      }  
    } else { # with matching
      K <- length(ITTY$diff)
      Cov <- cov(ITTY$diff, ITTD$diff)/K
      CACEvar <- (var(ITTY$diff)*(ITTD$est^2)/K +
                  var(ITTD$diff)*(ITTY$est^2)/K -
                  2*Cov**ITTY$est*ITTD$est)/(ITTD$est^4)
    }
  }
  return(list(est = CACEest, var = CACEvar, ITTd = ITTD, ITTy = ITTY,
              cov = Cov)) 
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
