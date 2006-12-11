ATEbounds <- function(formula, data = parent.frame(), alpha = 0.05,
                      maxY = NULL, minY = NULL) {
  
  ## getting Y and D
  call <- match.call()
  tm <- terms(formula)
  attr(tm, "intercept") <- 0
  mf <- model.frame(tm, data = data, na.action = 'na.pass')
  D <- model.matrix(tm, data = mf)
  if (max(D) > 1 || min(D) < 0)
    stop("the treatment variable should be a factor variable.")
  M <- ncol(D)
  Y <- model.response(mf)
  if (is.null(maxY))
    maxY <- max(Y, na.rm = TRUE)
  if (is.null(minY))
    minY <- min(Y, na.rm = TRUE)
  
  ## Bounds for average Ys
  bounds.Y <- ci.Y <- vars.Y <- matrix(NA, ncol = 2, nrow = M) 
  for (i in 1:M) {
    Ysub <- Y[D[,i]==1]
    n <- length(Ysub)
    Ymax <- Ymin <- Ysub
    Ymax[is.na(Ysub)] <- maxY
    Ymin[is.na(Ysub)] <- minY
    ## point estimates of the bounds
    bounds.Y[i,] <- c(mean(Ymin), mean(Ymax))
    ## variances
    vars.Y[i,] <- c(var(Ymin)/n, var(Ymax)/n)
    ## Bonferroni bounds
    ci.Y[i,] <- c(mean(Ymin) - qnorm(1-alpha/2)*sqrt(var(Ymin)/n),
                  mean(Ymax) + qnorm(1-alpha/2)*sqrt(var(Ymax)/n))
  }

  ## Bounds for the ATE
  bounds <- ci <- matrix(NA, ncol = 2, nrow = choose(M, 2))
  counter <- 1
  tmp <- NULL
  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      bounds[counter,] <- c(bounds.Y[i,1]-bounds.Y[j,2],
                            bounds.Y[i,2]-bounds.Y[j,1])
      ci[counter,] <- c(bounds[counter,1] -
                        qnorm(1-alpha/2)*sqrt(vars.Y[i,1]+vars.Y[j,2]),
                        bounds[counter,2] +
                        qnorm(1-alpha/2)*sqrt(vars.Y[i,2]+vars.Y[j,1]))
      counter <- counter + 1
      tmp <- c(tmp, paste(colnames(D)[i], "-", colnames(D)[j]))
    }
  }

  ## dimnames
  rownames(bounds) <- rownames(ci) <- tmp
  rownames(bounds.Y) <- rownames(ci.Y) <- colnames(D)
  colnames(bounds) <- colnames(bounds.Y) <- c("lower", "upper")
  colnames(ci) <- colnames(ci.Y) <-
    c(paste("lower ", alpha/2, "%CI", sep=""),
      paste("upper ", 1-alpha/2, "%CI", sep=""))
  
  return(list(bounds = bounds, ci = ci,
              bounds.Y = bounds.Y,
              ci.Y = ci.Y), maxY = maxY, minY = minY)
}
