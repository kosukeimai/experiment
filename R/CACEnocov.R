CACEnocov <- function(Y, D, Z, data = parent.frame(), 
                      match = NULL) {

  ## get the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  match <- eval(call$match, envir = data)
  
  ## point estimates and variances
  ITTY <- ATEnocov(Y = Y, Z = Z, grp = grp, match = match)
  ITTD <- ATEnocov(Y = D, Z = Z, grp = grp, match = match)
  CACEest <- ITTY$est/ITTD$est

  ## covariance calculation
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

  ## returning the results
  res <- list(est = CACEest, var = CACEvar, ITTd = ITTD, ITTy = ITTY,
              cov = Cov)
  class(res) <- "CACEnocov"
  return(res)
}
