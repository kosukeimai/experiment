###
### Calculate the ITT estimate with optional clustering or matching
###

ITTnocov <- function(Y, Z, grp = NULL, match = NULL){
  
  ITT <- mean(Y[Z==1])-mean(Y[Z==0])
  if (is.null(match)) {
    if (is.null(grp))
      ITT.var <- var(Y[Z==1])/sum(Z==1) + var(Y[Z==0])/sum(Z==0)
    else
      ITT.var <- varCluster(Y[Z==1], grp[Z==1])$var +
        varCluster(Y[Z==0], grp[Z==0])$var
  } else {
    if (is.null(grp)) {
      if (sum(Z==1) != sum(Z==0))
        stop("invalid input for `match'.")
      Y1 <- Y[Z==1]
      Y0 <- Y[Z==0]
      ind0 <- sort(match[Z==0], index.return = TRUE)
      ind1 <- sort(match[Z==1], index.return = TRUE)
      if (sum(ind0$x == ind1$x) != sum(Z==1))
        stop("invalid input for `match'.")
      D <- Y1[ind1$ix] - Y0[ind0$ix]
      ITT.var <- var(D)/length(D)
    } else {
      stop("invalid input for `match'.")
    }
  }
  return(list(est = ITT, var = ITT.var))
}

