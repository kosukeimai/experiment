bNegBin <- function(Y, X, beta.start, sig2.start, sims, beta0, A0,
                    a0, b0, tune.beta, tune.sig2) {
  
  ncov <- ncol(X)
  counter <- rep(0, 2)
  tmp <- .C("R2bNegBin", as.integer(Y), as.double(X),
            as.double(beta.start), as.double(sig2.start),
            as.integer(nrow(X)), as.integer(ncov), as.integer(sims), 
            as.double(beta0), as.double(A0), as.double(a0),
            as.double(b0), as.double(tune.beta), 
            as.double(tune.sig2), 
            betaStore = double(sims*ncov), sig2Store = double(sims),
            counter = as.integer(counter), PACKAGE = "experiment")

  res <- list(beta = matrix(tmp$betaStore, byrow = TRUE, ncol = ncov),
              sig2 = tmp$sig2Store, accept = tmp$counter/sims)
  
  return(res)
}
