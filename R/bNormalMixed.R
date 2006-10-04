bNormalMixed <- function(formula, formulaR, grp, data = parent.frame(),
                         coef.start = 0, sig2.start = 1, Phi.start = 1,
                         p.mean = 0, p.prec = 0.01, nu0 = 1, s0 = 1,
                         p.df = 1, p.scale = 1, p.improper = TRUE, n.draws = 5000) {

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
 
  gamma.start <- rnorm(n.gamma*length(unique(grp)))
    
  ## this code assumes the equal number of obs within each group
  res <- .C("R2bNormalMixedGibbs", as.double(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(beta.start),
            as.double(gamma.start), as.double(sig2.start),
            as.double(solve(Psi.start)),
            as.integer(nrow(X)), as.integer(n.beta), as.integer(n.gamma),
            as.integer(length(table(grp))), 
            as.integer(max(table(grp))), as.double(p.mean), as.double(p.prec),
            as.integer(p.improper), as.integer(nu0), as.integer(s0),
            as.integer(p.df), as.double(p.scale), as.integer(n.draws),
            betaStore = double(n.draws*n.beta),
            gammaStore = double(n.draws*n.gamma*ngrp),
            sig2Store = double(n.draws),
            PsiStore = double(n.draws*n.gamma*(n.gamma+1)/2),
            PACKAGE = "experiment")

  return(list(beta = matrix(res$betaStore, byrow = TRUE, ncol = n.beta),
              sig2 = res$sig2Store,
              gamma = array(res$gammaStore, dim = c(n.gamma, ngrp, n.draws)),
              Psi = matrix(res$PsiStore, byrow = TRUE, nrow = n.draws)))
}
