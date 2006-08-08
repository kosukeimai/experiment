boprobit <- function(Y, X, beta.start, tau.start, sims, beta0, A0) {

  tmp <- .C("R2boprobit", as.integer(Y), as.double(X),
            as.double(beta.start), as.double(tau.start),
            as.integer(nrow(X)), as.integer(ncol(X)),
            as.integer(max(Y)+1), as.double(beta0), as.double(A0),
            as.integer(sims), 
            betaStore = double(sims*ncol(X)),
            tauStore = double(sims*max(Y)),
            PACKAGE = "experiment")

  res <- list(beta = matrix(tmp$betaStore, byrow = TRUE, ncol = ncol(X)),
              tau = matrix(tmp$tauStore, byrow = TRUE, ncol = max(Y)))
  
  return(res)
}
