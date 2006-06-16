###
### Bayesian probit with nonignorable missing outcomes
###

NIbprobit <- function(Y, D, X, data = parent.frame(),
                      n.draws = 5000, param = TRUE, mda = TRUE,
                      p.mean.o = 0, p.var.o = 1000,
                      p.mean.r = 0, p.var.r = 1000,
                      coef.start.o = 0, coef.start.r = 0,
                      burnin = 0, thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  tm <- terms(X)
  attr(tm, "intercept") <- 0
  X <- model.matrix(tm, data = data, na.action = 'na.pass')
  ind <- complete.cases(cbind(D, X))
  Y <- Y[ind]
  D <- D[ind]
  X <- X[ind,]
  R <- (!is.na(Y))*1
  Y[is.na(Y)] <- rbinom(sum(is.na(Y)), size = 1, prob = 0.5)
  Xo <- cbind(1-D, D, X)
  Xr <- cbind(1-Y, Y, X)
  
  res <- list(call = call, Y = Y, X = X, D = D, n.draws = n.draws)

  n <- length(Y)
  k <- ncol(Xo)
  ## starting values
  if(length(coef.start.o) != k)
    coef.start.o <- rep(coef.start.o, k)
  if(length(coef.start.r) != k)
    coef.start.r <- rep(coef.start.r, k)
 
  ## prior
  if(length(p.mean.o) != k)
    p.mean.o <- rep(p.mean.o, k)
  if(length(p.mean.r) != k)
    p.mean.r <- rep(p.mean.r, k)
  if(!is.matrix(p.var.o))
    p.var.o <- diag(p.var.o, k)
  if(!is.matrix(p.var.r))
    p.var.r <- diag(p.var.r, k)
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function to do MCMC
  par <- .C("NIbprobit",
            as.integer(Y), as.integer(R), as.integer(D), 
            as.double(Xo), as.double(Xr),
            as.double(coef.start.o), as.double(coef.start.r),
            as.integer(n), as.integer(k), 
            as.double(p.mean.o), as.double(p.mean.r),
            as.double(solve(p.var.o)), as.double(solve(p.var.r)), 
            as.integer(param), as.integer(mda),
            as.integer(n.draws), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coef.o = double(k*(ceiling((n.draws-burnin)/keep))),
            coef.r = double(k*(ceiling((n.draws-burnin)/keep))),
            ATE = double(3*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="are")
  if (param) {
    res$coef.o <- matrix(par$coef.o, byrow = TRUE, ncol = k)
    res$coef.r <- matrix(par$coef.r, byrow = TRUE, ncol = k)
  }
  res$ATE <- matrix(par$ATE, byrow = TRUE, ncol = 3)
  
  class(res) <- "Classical.bprobit"
  return(res)
}
