###
### Calculate the ITT estimate with optional clustering
###


ITTnocov <- function(Y, Z, grp = NULL){
  
  ITT <- mean(Y[Z==1])-mean(Y[Z==0])
  if (is.null(grp))
    ITT.var <- var(Y[Z==1])/sum(Z==1) + var(Y[Z==0])/sum(Z==0)
  else
    ITT.var <- varCluster(Y[Z==1], grp[Z==1])$var +
      varCluster(Y[Z==0], grp[Z==0])$var

  return(list(est = ITT, var = ITT.var))
}
