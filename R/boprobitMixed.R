boprobitMixed <- function(formula, formulaR, grp, data = parent.frame(),
                          coef.start = 0, tau.start = NULL,
                          Phi.start = 1, p.mean = 0, p.prec = 0.01,
                          p.df = 1, p.scale = 1, prop = 0.01, mh = TRUE,
                          n.draws = 5000) {

  call <- match.call()
  mf <- model.frame(formula, data = data)
  X <- model.matrix(formula, data = mf)
  Y <- model.response(mf)
  grp <- as.integer(factor(eval(call$grp, data)))-1
  n.beta <- ncol(X)
  n.gamma <- ncol(Z)

  ## prior parameters
  if (length(p.mean) == 1)
    p.mean <- rep(p.mean, n.beta)
  else if (length(p.mean) != n.beta)
    stop(paste("the length of p.mean should be", n.beta))        

  if(is.matrix(p.prec)) {
    if (sum(dim(p.prec) == rep(n.beta, 2)) < 2)
      stop(paste("the dimension of p.prec should be",
                 rep(n.beta, 2)))    
  } else if (length(p.prec) == 1){
    p.prec <- diag(p.prec, n.beta)
  } else {
    stop("Incorrect input for p.prec")
  }

  if(is.matrix(p.scale)) {
    if (sum(dim(p.scale) == rep(n.gamma, 2)) < 2)
      stop(paste("the dimension of p.scale should be",
                 rep(n.gamma, 2)))    
  } else if (length(p.scale) == 1){
    p.scale <- diag(p.scale, n.gamma)
  } else {
    stop("Incorrect input for p.scale")
  }

  if (length(prop) == 1)
    prop <- rep(prop, max(Y)-1)
  else if (length(prop) != 2) 
    stop("Incorrect input for prop")
    
  ## starting values
  if (length(coef.start) == 1)
    coef.start <- rep(coef.start, n.beta)
  else if (length(coef.start) != n.beta)
    stop(paste("the length of coef.start should be", n.beta))        
  
  if(is.matrix(Phi.start)) {
    if (sum(dim(Phi.start) == rep(n.gamma, 2)) < 2)
      stop(paste("the dimension of Phi.start should be",
                 rep(n.gamma, 2)))    
  } else if (length(Phi.start) == 1){
    Phi.start <- diag(Phi.start, n.gamma)
  } else {
    stop("Incorrect input for Phi.start")
  }
  
  if (is.null(tau.start))
    tau.start <- c(seq(from = 0, to = max(Y)-1, by = 1), 1000)
  else if (length(tau.start) != (max(Y)+1))
    stop("Incorrect input for tau.start")
  
  ## this code assumes the equal number of obs within each group
  accept <- 0
  res <- .C("R2boprobitMixedMCMC", as.integer(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(coef.start),
            as.double(tau.start),
            as.double(Psi.start), as.integer(nrow(X)),
            as.integer(ncol(X)), as.integer(ncol(Z)),
            as.integer(length(table(grp))), as.integer(max(table(grp))),
            as.integer(max(Y)+1),
            as.double(p.mean), as.double(p.prec),
            as.integer(df), as.double(T0), as.integer(mh),
            as.double(prop), accept = as.integer(accept), 
            as.integer(n.draws), betaStore = double(n.draws*ncol(X)),
            gammaStore = double(n.draws*ncol(Z)*ngrp),
            tauStore = double(n.draws*max(Y)),
            PsiStore = double(n.draws*ncol(Z)*(ncol(Z)+1)/2),
            PACKAGE = "experiment")

  return(list(beta = matrix(res$betaStore, byrow = TRUE, ncol = ncol(X)),
              gamma = array(res$gammaStore, dim = c(ncol(Z), ngrp, n.draws)),
              Psi = matrix(res$PsiStore, byrow = TRUE, nrow = n.draws),
              tau = matrix(res$tauStore, byrow = TRUE, ncol = max(Y)),
              accept = res$accept/n.draws))
}
