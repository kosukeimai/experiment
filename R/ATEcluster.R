ATEcluster <- function(Y, Z, grp, data = parent.frame(),
                       match = NULL, weights = NULL,
                       method = "conditional") {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  N <- length(Y)
  M <- length(unique(match))
  
  ###
  ### conditional ATE
  ###
  Y1bar <- tapply(Y[Z==1], match[Z==1], mean)
  Y0bar <- tapply(Y[Z==0], match[Z==0], mean)
  Y1var <- tapply(Y[Z==1], match[Z==1], var)/tapply(rep(1, sum(Z==1)),
                      match[Z==1], sum)
  Y0var <- tapply(Y[Z==0], match[Z==0], var)/tapply(rep(1, sum(Z==0)),
                      match[Z==0], sum)
  if (is.null(weights))
    weights <- rep(1, N)

  ## unbiased estimation
  w1 <- tapply(weights[Z==1], match[Z==1], mean)
  w0 <- tapply(weights[Z==0], match[Z==0], mean)
  tmp <- sum(c(w1, w0))
  w1 <- N*w1/tmp
  w0 <- N*w0/tmp
  ATE.est1 <- sum(w1*Y1bar - w0*Y0bar)/N
  ATE.var1 <- sum(w1^2*Y1var+w0^2*Y0var)/N^2

  ## weighted estimation
  w <- w1 + w0
  w <- N*w/sum(w)
  ATE.est <- sum(w*(Y1bar-Y0bar))/N
  ATE.var <- sum(w^2*(Y1var+Y0var))/N^2
  
  return(list(est = ATE.est, estU = ATE.est1,
              var = ATE.var, varU = ATE.var1, N = N, M = M,
              Y1bar = Y1bar, Y0bar = Y0bar, Y1var = Y1var,
              Y0var = Y0var, weights = w, weights1 = w1, weights0 = w0))
}
