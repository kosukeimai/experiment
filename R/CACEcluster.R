CACEcluster <- function(Y, D, Z, grp, data = parent.frame(),
                        match = NULL, estimand = "CATE1") {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)

  if (estimand == "CATE1")
    ITTestimand <- "PATE1"
  else if (estimand == "CATE2")
    ITTestimand <- "PATE2"
  else
    stop("invalid estimand")
  
  ITTY <- ATEcluster(Y = Y, Z = Z, grp = grp, match = match,
                     estimand = ITTestimand)
  ITTD <- ATEcluster(Y = D, Z = Z, grp = grp, match = match,
                     estimand = ITTestimand)

  ITTY.est <- ITTY$est
  ITTD.est <- ITTD$est
  CACE.est <- ITTY.est/ITTD.est
  CACE.bias <- CACE.est - (ITTY.est - ITTY$bias)/(ITTD.est - ITTD$bias)
  N <- ITTY$N
  M <- ITTY$M
  diffY <- ITTY$Y1bar-ITTY$Y0bar
  diffD <- ITTD$Y1bar-ITTD$Y0bar
  w <- ITTY$weights

  if (is.null(match))
    stop("the estimator is not yet available.")
  else {
    Cov <- M*sum((diffY*w - N*ITTY.est/M) *
                 (diffD*w - N*ITTD$est/M))/((M-1)*(N^2))
    CACE.var <- (ITTY$var*(ITTD.est^2) + ITTD$var*(ITTY.est^2) -
                2*Cov*ITTY.est*ITTD.est)/(ITTD.est^4)
  }

  ## return the results
  res <- list(est = CACE.est, bias = CACE.bias, var = CACE.var,
              cov = Cov, N = N, M = M, ITTY = ITTY, ITTD = ITTD,
              weights = w)
  class(res) <- "CACEcluster"
  return(res)
}
