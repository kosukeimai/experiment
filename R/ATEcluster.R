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
  Y1bar <- tapply(Y[Z==1], match, mean)
  Y0bar <- tapply(Y[Z==0], match, mean)
  Y1var <- tapply(Y[Z==1], match, var)/tapply(rep(1, sum(Z==1)),
                      match, sum)
  Y0var <- tapply(Y[Z==0], match, var)/tapply(rep(1, sum(Z==0)),
                      match, sum)
  w <- tapply(weights, match, mean)
  ATE.est <- sum(w*(Y1bar-Y0bar))/N
  ATE.var <- sum(w^2*(Y1var+Y0var))/N^2

  return(list(est = ATE.est, var = ATE.var, N = N, M = M, weights = w))
}
