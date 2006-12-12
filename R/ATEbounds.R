###
### Calculates the bounds for the ATE in the presence of missing
### response 
###

ATEbounds <- function(formula, data = parent.frame(), maxY = NULL,
                      minY = NULL, alpha = 0.05, n.reps = 0, ...) {

  ## getting Y and D
  call <- match.call()
  tm <- terms(formula)
  attr(tm, "intercept") <- 0
  mf <- model.frame(tm, data = data, na.action = 'na.pass')
  D <- model.matrix(tm, data = mf)
  M <- ncol(D)
  if (max(D) > 1 || min(D) < 0)
    stop("the treatment variable should be a factor variable.")
  Y <- model.response(mf)
  if (is.null(maxY))
    maxY <- max(Y, na.rm = TRUE)
  if (is.null(minY))
    minY <- min(Y, na.rm = TRUE)

  ### computing the bounds
  res <- boundsComp(cbind(Y, D), rep(1, length(Y)), maxY, minY, alpha = alpha)
  
  ## CI based on B-method
  if (n.reps > 0) {
    breps <- boot(data = cbind(Y, D), statistic = boundsComp,
                  R = n.reps, maxY = maxY, minY = minY, alpha = NULL)$t
    res$bmethod.ci.Y <- matrix(NA, ncol = 2, nrow = M)
    res$bmethod.ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
    for (i in 1:M) { 
      tmp <- boundsCI(breps[,(i-1)*M+1], breps[,i*M],
                      res$bounds.Y[i,1],
                      res$bounds.Y[i,2], alpha)
      res$bmethod.ci.Y[i,] <- tmp$bmethod
      res$bonf.ci.Y[i,] <- tmp$bonferroni
    }
    counter <- 1
    for (i in 1:(M-1)) 
      for (j in (i+1):M) {
        tmp <- boundsCI(breps[,2*M+counter], breps[,2*M+counter+1],
                        res$bounds[(counter+1)/2,1],
                        res$bounds[(counter+1)/2,2], alpha)
        res$bmethod.ci[(counter+1)/2,] <- tmp$bmethod
        res$bonf.ci[(counter+1)/2,] <- tmp$bonferroni
        counter <- counter + 2
      }
  }
  
  ## dimnames
  rownames(res$bounds.Y) <- rownames(res$bonf.ci.Y) <-
    rownames(res$bmethod.ci.Y) <- colnames(D)
  tmp <- NULL
  for (i in 1:(M-1)) 
    for (j in (i+1):M) 
      tmp <- c(tmp, paste(colnames(D)[i], "-", colnames(D)[j]))
  rownames(res$bounds) <- rownames(res$bonf.ci) <-
    rownames(res$bmethod.ci) <- tmp
  colnames(res$bounds) <- colnames(res$bounds.Y) <- c("lower", "upper")
  colnames(res$bonf.ci) <- colnames(res$bonf.ci.Y) <-
    colnames(res$bmethod.ci) <- colnames(res$bmethod.ci.Y) <-
      c(paste("lower ", alpha/2, "%CI", sep=""),
        paste("upper ", 1-alpha/2, "%CI", sep=""))
  
  return(res)
}

###
### An internal function which computes the bounds and bonferroni CI
### if alpha is specified (when alpha = NULL, then it returns a vector
### of bounds; this is used for bootstrap)
###

boundsComp <- function(data, weights, maxY, minY, alpha = NULL) {
  Y <- data[,1]
  D <- data[,-1]
  M <- ncol(D)
  bounds.Y <- ci.Y <- vars.Y <- matrix(NA, ncol = 2, nrow = M)
  for (i in 1:M) {
    Ysub <- Y[D[,i]==1]
    w <- weights[D[,i]==1]
    n <- length(Ysub)
    Ymax <- Ymin <- Ysub
    Ymax[is.na(Ysub)] <- maxY
    Ymin[is.na(Ysub)] <- minY
    ## point estimates of the bounds
    bounds.Y[i,] <- c(weighted.mean(Ymin, w), weighted.mean(Ymax, w))
    if (!is.null(alpha)) {
      ## variances
      vars.Y[i,] <- c(weighted.var(Ymin, w)*sum(w^2)/(sum(w)^2),
                      weighted.var(Ymax, w)*sum(w^2)/(sum(w)^2))
      ## Bonferroni bounds
      ci.Y[i,] <- c(mean(Ymin) - qnorm(1-alpha/2)*sqrt(vars.Y[i,1]),
                    mean(Ymax) + qnorm(1-alpha/2)*sqrt(vars.Y[i,2]))
    }
  }
  
  ## Bounds for the ATE
  bounds <- ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
  counter <- 1
  tmp <- NULL
  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      bounds[counter,] <- c(bounds.Y[i,1]-bounds.Y[j,2],
                            bounds.Y[i,2]-bounds.Y[j,1])
      if (!is.null(alpha))
        ci[counter,] <- c(bounds[counter,1] -
                          qnorm(1-alpha/2)*sqrt(vars.Y[i,1]+vars.Y[j,2]),
                          bounds[counter,2] +
                          qnorm(1-alpha/2)*sqrt(vars.Y[i,2]+vars.Y[j,1]))
      counter <- counter + 1
      tmp <- c(tmp, paste(colnames(D)[i], "-", colnames(D)[j]))
    }
  }
  if (is.null(alpha))
    return(c(t(rbind(bounds.Y, bounds))))
  else
    return(list(bounds.Y = bounds.Y, bounds = bounds, bonf.ci = ci,
                bonf.ci.Y = ci.Y, maxY = maxY, minY = minY))
}

