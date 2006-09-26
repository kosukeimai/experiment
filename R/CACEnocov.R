###
### Calculate the CACE with optional clustering
###

CACEnocov <- function(Y, D, Z, grp = NULL){
  
  N0 <- sum(Z==0)
  N1 <- sum(Z==1)
  ITTD <- ITTnocov(D, Z, grp)
  ITTY <- ITTnocov(Y, Z, grp)
  if (is.null(grp)) {
    Cov0 <- cov(Y[Z==0], D[Z==0])/N0
    Cov1 <- cov(Y[Z==1], D[Z==1])/N1
  } else {
    Cov0 <- covCluster(Y[Z==0], D[Z==0], grp[Z==0])
    Cov1 <- covCluster(Y[Z==1], D[Z==1], grp[Z==1])
  }
  
  est <- ITTY$est/ITTD$est
  var <- (ITTY$var*(ITTD$est^2) + ITTD$var*(ITTY$est^2) -
          2*(Cov0$cov + Cov1$cov)*ITTY$est*ITTD$est)/(ITTD$est^4)

  return(list(est = est, var = var, ITTD = ITTD, ITTY = ITTY,
              cov = Cov0$cov + Cov1$cov)) 
}
