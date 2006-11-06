###
### Bayesian probit with nonignorable missing outcomes and
### multi-valued treatments
###

NIbprobitMixed <- function(formula, formulae.o, formulae.r, grp,
                           data = parent.frame(), n.draws = 5000,
                           insample = FALSE, param = TRUE, mda = TRUE,
                           p.mean.o = 0, p.prec.o = 0.01,
                           p.mean.r = 0, p.prec.r = 0.01,
                           p.df.o = 1, p.df.r = 1,
                           p.scale.o = 1, p.scale.r = 1,
                           coef.start.o = 0, coef.start.r = 0,
                           Psi.start.o = 1, Psi.start.r = 1,
                           burnin = 0, thin = 0, verbose = TRUE) {

  ## getting Y and D
  call <- match.call()
  tm <- terms(formula)
  attr(tm, "intercept") <- 0
  mf <- model.frame(tm, data = data, na.action = 'na.pass')
  D <- model.matrix(tm, data = mf)
  if (max(D) > 1 || min(D) < 0)
    stop("the treatment variable should be a factor variable.")
  Y <- model.response(mf)
  m <- ncol(D) # number of treatment levels including control

  ## getting fixed effects
  tm <- terms(formulae.o[[1]])
  attr(tm, "intercept") <- 1
  Xo <- model.matrix(tm, data = data, na.action = 'na.pass')
  Xo <- Xo[,(colnames(Xo) != "(Intercept)")]
  tm <- terms(formulae.r[[1]])
  attr(tm, "intercept") <- 1
  Xr <- model.matrix(tm, data = data, na.action = 'na.pass')
  Xr <- Xr[,(colnames(Xr) != "(Intercept)")]

  ## getting random effects
  tm <- terms(formulae.o[[2]])
  Zo <- model.matrix(tm, data = data, na.action = 'na.pass')
  Zo <- Zo[,(colnames(Zo) != "(Intercept)")]
  tm <- terms(formulae.r[[2]])
  Zr <- model.matrix(tm, data = data, na.action = 'na.pass')
  Zr <- Zr[,(colnames(Zr) != "(Intercept)")]
    
  ## taking care of NA's in D, X, and Z
  ind <- complete.cases(cbind(D, Xo, Xr, Zo, Zr))
  Y <- Y[ind]
  D <- D[ind,]
  Xo <- Xo[ind,]
  Xr <- Xr[ind,]
  Zo <- Zo[ind,]
  Zr <- Zr[ind,]
  R <- (!is.na(Y))*1
  Y[is.na(Y)] <- rbinom(sum(is.na(Y)), size = 1, prob = 0.5)
  cnameso <- c(colnames(D), colnames(Xo))
  cnamesr <- c("1-Y", "Y", colnames(Xr))
  Xo <- cbind(D, Xo)
  colnames(Xo) <- cnameso
  Xr <- cbind(1-Y, Y, Xr)
  colnames(Xr) <- cnamesr

  ## group indicators
  grp <- eval(call$grp, envir = data)
  grp <- as.integer(factor(grp))-1
  ngrp <- length(table(grp))

  ## output
  res <- list(call = call, Y = Y, Xo = Xo, Xr = Xr, Zo = Zo, Zr = Zr,
              n.draws = n.draws, grp = grp)

  ## some important numbers
  n <- length(Y)
  ncovo <- ncol(Xo)
  ncovr <- ncol(Xr)
  ncovo1 <- ncol(Zo)
  ncovr1 <- ncol(Zr)
  
  ## starting values
  if(length(coef.start.o) != ncovo)
    coef.start.o <- rep(coef.start.o, ncovo)
  if(length(coef.start.r) != ncovr)
    coef.start.r <- rep(coef.start.r, ncovr)
  if(is.matrix(Psi.start.o)) {
    if (dim(Psi.start.o) != ncovo1)
      stop(paste("the dimension of Psi.start.o should be", ncovo1))    
  } else if (length(Psi.start.o) == 1)
    Psi.start.o <- diag(Psi.start.o, ncovo1)
  else if (length(Psi.start.o) == ncovo1)
    Psi.start.o <- diag(Psi.start.o)
  else
    stop("Incorrect input for Psi.start.o")

  if(is.matrix(Psi.start.r)) {
    if (dim(Psi.start.r) != ncovr1)
      stop(paste("the dimension of Psi.start.r should be", ncovr1))    
  } else if (length(Psi.start.r) == 1)
    Psi.start.r <- diag(Psi.start.r, ncovr1)
  else if (length(Psi.start.r) == ncovr1)
    Psi.start.r <- diag(Psi.start.r)
  else
    stop("Incorrect input for Psi.start.r")
  
  ## prior
  if(length(p.mean.o) != ncovo)
    p.mean.o <- rep(p.mean.o, ncovo)
  if(length(p.mean.r) != ncovr)
    p.mean.r <- rep(p.mean.r, ncovr)
  if(!is.matrix(p.prec.o))
    p.prec.o <- diag(p.prec.o, ncovo)
  if(!is.matrix(p.prec.r))
    p.prec.r <- diag(p.prec.r, ncovr)

  if(is.matrix(p.scale.o)) {
    if (dim(p.scale.o) != ncovo1)
      stop(paste("the dimension of p.scale.o should be", ncovo1))    
  } else if (length(p.scale.o) == 1)
    p.scale.o <- diag(p.scale.o, ncovo1)
  else if (length(p.scale.o) == ncovo1)
    p.scale.o <- diag(p.scale.o)
  else
    stop("Incorrect input for p.scale.o")

  if(is.matrix(p.scale.r)) {
    if (dim(p.scale.r) != ncovr1)
      stop(paste("the dimension of p.scale.r should be", ncovr1))    
  } else if (length(p.scale.r) == 1)
    p.scale.r <- diag(p.scale.r, ncovr1)
  else if (length(p.scale.r) == ncovr1)
    p.scale.r <- diag(p.scale.r)
  else
    stop("Incorrect input for p.scale.r")
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function to do MCMC
  par <- .C("NIbprobitMixed",
            as.integer(Y), as.integer(R), as.integer(grp),
            as.integer(max(grp)+1), as.integer(max(table(grp))),
            as.double(Xo), as.double(Xr), as.double(Zo), as.double(Zr), 
            as.double(coef.start.o), as.double(coef.start.r),
            as.double(Psi.start.o), as.double(Psi.start.r),
            as.integer(n), as.integer(ncovo), as.integer(ncovr),
            as.integer(ncovo1), as.integer(ncovr1), as.integer(m), 
            as.double(p.mean.o), as.double(p.mean.r),
            as.double(p.prec.o), as.double(p.prec.r),
            as.integer(p.df.o), as.integer(p.df.r),
            as.double(p.scale.o), as.double(p.scale.r),
            as.integer(insample), as.integer(param), as.integer(mda),
            as.integer(n.draws), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coef.o = double(ncovo*(ceiling((n.draws-burnin)/keep))),
            coef.r = double(ncovr*(ceiling((n.draws-burnin)/keep))),
            sPsiO = double(ncovo1*(ncovo1+1)*(ceiling((n.draws-burnin)/keep))/2),
            sPsiR = double(ncovr1*(ncovr1+1)*(ceiling((n.draws-burnin)/keep))/2),
            ATE = double((m-1)*(ceiling((n.draws-burnin)/keep))),
            BASE = double(m*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="experiment")
  
  if (param) {
    res$coef.o <- matrix(par$coef.o, byrow = TRUE, ncol = ncovo)
    colnames(res$coef.o) <- colnames(Xo)
    res$coef.r <- matrix(par$coef.r, byrow = TRUE, ncol = ncovr)
    colnames(res$coef.r) <- colnames(Xr)
  }
  res$ATE <- matrix(par$ATE, byrow = TRUE, ncol = m-1)
  res$base <- matrix(par$BASE, byrow = TRUE, ncol = m)
  
  class(res) <- "NIbprobit"
  return(res)
}
