ATEcluster <- function(Y, Z, grp, data = parent.frame(),
                       match = NULL, weights = NULL,
                       estimand = "PATE1") {

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
  diff <- Y1bar-Y0bar

  if (estimand == "PATE1") {
    ## PATE without sampling within clusters
    if (is.null(weights)) {
      weights <- rep(1, N)
    } else {
      stop("For this estimand, the cluster sample sizes will be used as weights.")
    }
    if (is.null(match)) {
      stop("This option is not yet available.")
    } else {
      w1 <- tapply(weights[Z==1], match[Z==1], sum)
      w0 <- tapply(weights[Z==0], match[Z==0], sum)
      w <- w1 + w0
      ## donner&klar weights: w <- w1*w0/(w1 + w0)
      w <- N*w/sum(w)
      ## weighted estimation
      ATE.est <- weighted.mean(diff, w)
      ATE.var <- M*sum((w*diff-N*ATE.est/M)^2)/((M-1)*(N^2))
      ## donner&klar formula: ATE.var <- sum(wT^2)*sum(wT*(diff-ATE.est)^2)/N^3
      ## unbiased estimation
      Y1sum <- tapply(Y[Z==1], match[Z==1], sum)
      Y0sum <- tapply(Y[Z==0], match[Z==0], sum)
      ATE.estU <- 2*(sum(Y1sum)-sum(Y0sum))/N
      #ATE.varU <- 4*M*var(Y1sum-Y0sum)/(N^2)
      ATE.varlb <- ATE.varub <- NULL
      Y1var <- Y0var <- NULL
    }
  } else  if (estimand %in% c("CATE", "PATE2")) {
    ## CATE, PATE with two-stage sampling
    if (is.null(match)) {
      stop("This option is not yet available.")
    } else {
      Y1var <- tapply(Y[Z==1], match[Z==1], var)/tapply(rep(1, sum(Z==1)),
                                       match[Z==1], sum)
      Y0var <- tapply(Y[Z==0], match[Z==0], var)/tapply(rep(1, sum(Z==0)),
                                       match[Z==0], sum)
      if (is.null(weights)) {
        if (estimand == "CATE")
          stop("CATE requires the specification of cluster weights")
        else
          weights <- rep(1, N)
      }
      
      ## unbiased estimation
      w1 <- tapply(weights[Z==1], match[Z==1], mean)
      w0 <- tapply(weights[Z==0], match[Z==0], mean)
      tmp <- sum(c(w1, w0))
      w1 <- N*w1/tmp
      w0 <- N*w0/tmp
      ATE.estU <- sum(w1*Y1bar - w0*Y0bar)/N
      ##ATE.var1 <- sum(w1^2*Y1var+w0^2*Y0var)/(N^2)
      ##ATE.var1 <- sum(w1^2*Y1var+w0^2*Y0var+(w1*Y1bar-w0*Y0bar)^2)/(N^2)
      
      ## weighted estimation
      w <- w1 + w0
      w <- N*w/sum(w)
      ATE.est <- sum((w/N)*diff)
      if (estimand == "PATE2") {
        ## two-stage sampling estimate
        ATE.var <- M*(mean((w^2)*(Y1var+Y0var)) + var(w*diff))/(N^2)
        ATE.varlb <- ATE.varub <- NULL
      } else {
        ATE.var <- NULL
        ## lower bound
        ATE.varlb <- sum((w/N)^2*(Y1var+Y0var))
        ## upper bound
        ATE.varub <- sum((w/N)^2*(Y1var+Y0var+(diff-mean(diff))^2))
      }
    }
  } else {
    stop("invalid input for `estimand'")
  }

  ## return the resutls
  res <- list(est = ATE.est, bias = ATE.est-ATE.estU,
              var = ATE.var, var.lb = ATE.varlb,
              var.ub = ATE.varub, N = N, M = M,
              Y1bar = Y1bar, Y0bar = Y0bar, Y1var = Y1var,
              Y0var = Y0var, weights = w, weights1 = w1,
              weights0 = w0)
  class(res) <- "ATEcluster"
  return(res)
}
