bNormalReg <- function(Y, X, beta.start, sig2.start, sims, beta0, A0,
                       pbeta = TRUE, psig2 = FALSE, sig2.fixed = FALSE) {

  tmp <- .C("R2bNormalReg", as.double(Y), as.double(X),
            as.double(beta.start), as.double(sig2.start),
            as.integer(nrow(X)), as.integer(ncol(X)), as.integer(sims), 
            as.integer(pbeta), as.double(beta0), as.double(A0),
            as.integer(psig2), as.double(1), as.integer(ncol(X)),
            as.integer(sig2.fixed),
            betaStore = double(sims*ncol(X)),
            sig2Store = double(sims),
            PACKAGE = "experiment")

  res <- list(beta = matrix(tmp$betaStore, byrow = TRUE, ncol = ncol(X)),
              sig2 = tmp$sig2Store)
  
  return(res)
}
