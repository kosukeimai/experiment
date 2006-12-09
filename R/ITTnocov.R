###
### Calculate the ITT estimate with optional clustering or matching
###

ITTnocov <- function(Y, Z, grp = NULL, match = NULL, size = NULL){

  weighted.var <- function(x, w) {
    w <- w/sum(w)*length(w)
    sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)
  }
  
  Y1 <- Y[Z==1]
  Y0 <- Y[Z==0]
  if (is.null(size))
    size <- rep(1, length(Y))
  else if (!is.null(grp))
    stop("`grp' should be NULL when `size' is specified")
  w1 <- size[Z==1]
  w0 <- size[Z==0]
  ITT <- weighted.mean(Y1, w1) - weighted.mean(Y0, w0)
  if (is.null(match)) {
    if (is.null(grp))
      ITT.var <- weighted.var(Y1, w1)/sum(Z==1) +
        weighted.var(Y0, w0)/sum(Z==0)
    else
      ITT.var <- varCluster(Y1, grp[Z==1])$var +
        varCluster(Y0, grp[Z==0])$var
  } else {
    if (is.null(grp)) {
      if (sum(Z==1) != sum(Z==0))
        stop("invalid input for `match'.")
      ind0 <- sort(match[Z==0], index.return = TRUE)
      ind1 <- sort(match[Z==1], index.return = TRUE)
      if (sum(ind0$x == ind1$x) != sum(Z==1))
        stop("invalid input for `match'.")
      D <- Y1[ind1$ix] - Y0[ind0$ix]
      w <- w1[ind1$ix] + w0[ind0$ix]
      ITT.var <- weighted.var(D, w)/length(D)
    } else {
      stop("invalid input for `match'.")
    }
  }
  
  return(list(est = ITT, var = ITT.var))
}

