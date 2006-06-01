AllorNone <- function(Y, D, Z, data = parent.frame()) {
  
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  N <- length(Y)
  N1 <- length(Y[Z==1])
  N0 <- length(Y[Z==0])
  R <- (!is.na(Y))*1
  if (sum(R)>0) {
    ## Frangakis and Rubin (1999) Biometrika
    U <- sum(D*Z)/sum(Z)
    R01 <- sum(R*(1-D)*Z)/sum((1-D)*Z) 
    R0 <- sum(R*(1-Z))/sum(1-Z)
    Y0 <- mean(Y[R==1 & Z==0])
    Y01 <- sum(na.omit(Y*R*(1-D)*Z))/sum(R*(1-D)*Z)
    Y10 <- (Y0*R0-Y01*R01*(1-U))/(R0-R01*(1-U))
    Y11 <- sum(na.omit(Y*R*D*Z))/sum(R*D*Z)
    q <- var(Y[R==1 & D==1 & Z==1])/mean(Z*D*R)
    v <- delta <- rep(0, 5)
    v[1] <- U*(1-U)/mean(Z==1)
    v[2] <- var(Y[R==1 & D==0 & Z == 1])/mean(Z*(1-D)*R==1) 
    v[3] <- R01*(1-R01)/mean(Z*(1-D)==1)
    v[4] <- R0*(1-R0)/mean(Z==0)
    v[5] <- var(Y[R==1 & Z==0])/mean(R*(1-Z)==1)
    w <- 1/(R0-R01*(1-U))
    delta[1] <- -R0*R01*(Y0-Y01)*(w^2)
    delta[2] <- -R01*(1-U)*w
    delta[3] <- R0*(Y0-Y01)*(1-U)*(w^2)
    delta[4] <- -R01*(Y0-Y01)*(1-U)*(w^2)
    delta[5] <- R0*w
    
    IVest <- Y11-Y10
    IVvar <- (q + sum(v*(delta^2)))/N
    ITTest <- U*(Y11-Y10)
    ITTvar <-
      ((U^2)*q+v[1]*(Y11-Y10-U*delta[1])^2+(U^2)*sum(v[2:5]*(delta[2:5]^2)))/N
  }
  else {
    ## Imbens and Rubin (1997) Annals of Statistics
    Y1bar <- mean(Y[Z==1])
    Y0bar <- mean(Y[Z==0])
    D1bar <- mean(D[Z==1])
    D0bar <- mean(D[Z==0])
    ITTest <- Y1bar - Y0bar
    ITTestD <- D1bar - D0bar
    IVest <- ITTest/ITTestD
    Y1var <- sum(Z*(Y-Y1bar))/(N1^2)
    Y0var <- sum((1-Z)*(Y-Y0bar))/(N0^2)
    D1var <- sum(Z*(D-D1bar))/(N1^2)
    D0var <- sum((1-Z)*(D-D0bar))/(N0^2)
    ITTvar <- Y1var + Y0var
    ITTvarD <- D1var + D0var
    ITTcov <- sum(Z*(Y-Y1bar)*(D-D1bar))/(N1^2) +
      sum((1-Z)*(Y-Y0bar)*(D-D0bar))/(N0^2)
    IVvar <- (ITTvar*(ITTestD^2)+ITTvarD*(ITTest^2)-2*ITTcov*ITTest*ITTestD)/(ITTestD^4)
  }
  return(list(CACEest=IVest, CACEse=sqrt(IVvar), ITTest=ITTest, ITTse=sqrt(ITTvar)))
}
