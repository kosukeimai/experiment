ATEcluster <- function(Y, Z, grp, data, match = NULL,
                       weights = NULL, method = "conditional") {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  N <- length(Y)
  M <- length(unique(match))
  Y1bar <- tapply(Y[Z==1], match[Z==1], mean)
  Y0bar <- tapply(Y[Z==0], match[Z==0], mean)

  ## conditional ATE
  Y1var <- tapply(Y[Z==1], match[Z==1], var)/tapply(rep(1, sum(Z==1)),
                      match[Z==1], sum)
  Y0var <- tapply(Y[Z==0], match[Z==0], var)/tapply(rep(1, sum(Z==0)),
                      match[Z==0], sum)
  if (is.null(weights))
    weights <- rep(1, N)
  w <- tapply(weights, match, mean)
  w <- N*w/sum(w)
  ATE.est <- sum(w*(Y1bar-Y0bar))/N
  ATE.var <- sum(w^2*(Y1var+Y0var))/N^2

  return(list(est = ATE.est, var = ATE.var, N = N, M = M, weights = w))
}
