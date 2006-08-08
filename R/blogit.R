blogit <- function(Y, X, beta.start, sims, beta0, A0, tune) {

  ncov <- ncol(X)
  ndim <- length(unique(Y))-1
  counter <- rep(0, ncov*ndim)

  tmp <- .C("R2logitMetro", as.integer(Y), as.double(X),
            as.double(beta.start), 
            as.integer(nrow(X)), as.integer(ndim), as.integer(ncov), 
            as.double(beta0), as.double(A0), as.double(tune),
            as.integer(sims), counter = as.integer(counter),
            store = double(sims*ncol(X)*ndim),
            PACKAGE = "experiment")

  res <- list(beta = matrix(tmp$store, byrow = TRUE, ncol =
                ndim*ncov), accept = tmp$counter/sims)
  
  return(res)
}
