bNegBin <- function(Y, X, beta.start, sig2.start, sims, beta0, A0,
                    a0, b0, tune.beta, tune.sig2) {
  
  counter <- rep(0, 2)
  tmp <- .C("R2bNegBin", as.integer(Y), as.double(X),
            as.double(beta.start), as.double(sig2.start),
            as.integer(nrow(X)), as.integer(ncol(X)), as.integer(sims), 
            as.double(beta0), as.double(A0), as.double(tune),
            as.integer(sims), counter = as.integer(counter),
            store = double(sims*ncol(X)*ndim),
            PACKAGE = "experiment")

  res <- list(beta = matrix(tmp$store, byrow = TRUE, ncol =
                ndim*ncov), accept = tmp$counter/sims)
  
  return(res)
}
