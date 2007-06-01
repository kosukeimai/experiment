CACEcluster <- function(Y, D, Z, grp, data = parent.frame(),
                        match = NULL, weights = NULL,
                        method = "conditional") {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  ITTY <- ATEcluster(Y = Y, Z = Z, grp = grp, match = match, weights =
                     weights, method = method)
  ITTD <- ATEcluster(Y = D, Z = Z, grp = grp, match = match, weights =
                     weights, method = method)

  CACE.est <- ITTY$est/ITTD$est
  CACE.est1 <- ITTY$estU/ITTD$estU
  N <- length(Y)
  M <- length(unique(match))
  w <- ITTY$weights
  w1 <- ITTY$weights1
  w0 <- ITTY$weights0
  
  ## conditional ATE
  Y1 <- Y[Z==1]
  Y0 <- Y[Z==0]
  D1 <- D[Z==1]
  D0 <- D[Z==0]
  allmatch <- unique(match)
  Y1cov <- Y0cov <- NULL
  tmp1 <- match[Z==1]
  tmp0 <- match[Z==0]
  for (i in 1:M) {
    tmp <- allmatch[i]
    Y1cov <- c(Y1cov, cov(Y1[(tmp1 == tmp)],
                          D1[(tmp1 == tmp)])/length(Y1[(tmp1 == tmp)]))
    Y0cov <- c(Y0cov, cov(Y0[tmp0 == tmp],
                          D0[tmp0 == tmp])/length(Y0[(tmp0 == tmp)]))
  }
               
  CACE.cov <- sum(w1^2*Y1cov+w0^2*Y0cov)/N^2
  CACE.cov1 <- sum(w^2*(Y1cov+Y0cov))/N^2
  CACE.var <- (ITTD$est^2*ITTY$var + ITTY$est^2*ITTD$var -
               2*ITTD$est*ITTY$est*CACE.cov)/ITTD$est^4
  CACE.var1 <- (ITTD$estU^2*ITTY$varU + ITTY$estU^2*ITTD$varU -
                2*ITTD$estU*ITTY$estU*CACE.cov1)/ITTD$estU^4
  
  return(list(est = CACE.est, estU = CACE.est1, var = CACE.var,
              varU = CACE.var1, cov = CACE.cov, covU = CACE.cov1,
              N = N, M = M, ITTY = ITTY, ITTD = ITTD,
              weights = w, weights1 = w1, weights0 = w0))
}
