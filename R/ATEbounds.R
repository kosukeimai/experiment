ATEbounds <- function(Y, D, data = parent.frame(), alpha = 0.95) {
  
  ## prep
  Y <- eval(Y, data)
  D <- eval(D, data)
  Y1 <- Y[D==1]
  Y0 <- Y[D==0]
  n1 <- length(Y1)
  n0 <- length(Y0)
  maxY <- max(Y, na.rm = TRUE)
  minY <- min(Y, na.rm = TRUE)
  Y1max <- Y1min <- Y1
  Y1max[is.na(Y1)] <- maxY
  Y1min[is.na(Y1)] <- minY
  Y0max <- Y0min <- Y0
  Y0max[is.na(Y0)] <- maxY
  Y0min[is.na(Y0)] <- minY
  
  ## point estimates of the bounds
  ub <- mean(Y1max)-mean(Y0min)
  lb <- mean(Y1min)-mean(Y0max)

  ## Bonferroni bounds
  ub.ci <- ub + qnorm(1-alpha/2)*sqrt(var(Y1max)/n1+var(Y0min)/n0)
  lb.ci <- lb - qnorm(1-alpha/2)*sqrt(var(Y1min)/n1+var(Y0max)/n0)
  
  return(list(upper = ub, lower = lb, upper.ci = ub.ci, lower.ci = lb.ci))
}
