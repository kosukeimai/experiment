Noncomp.bprobitMixed <- function(formulae, Z, D, grp, data = parent.frame(),
                                 n.draws = 5000, param = TRUE,
                                 in.sample = FALSE, model.c = "probit", 
                                 tune.c = 1, 
                                 p.mean.c = 0, p.prec.c = 1000, p.mean.o = 0,
                                 p.prec.o = 1000, p.mean.r = 0, p.prec.r = 1000,
                                 mda = TRUE,
                                 coef.start.c = 0, coef.start.o = 0,
                                 coef.start.r = 0, Psi.start.c = 1,
                                 Psi.start.o = 1, Psi.start.r = 1,
                                 p.df.c = 1, p.df.o = 1, p.df.r = 1,
                                 p.scale.c = 1, p.scale.o = 1,
                                 p.scale.r = 1,
                                 burnin = 0, thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()
  ## outcome model: fixed effects
  mf <- model.frame(formulae[[1]], data=data, na.action='na.pass')
  Xo <- model.matrix(formulae[[1]], data=mf)
  if (sum(is.na(Xo)) > 0)
    stop("missing values not allowed in covariates")
  Y <- as.integer(model.response(mf))
  ## outcome model: random effects
  mf <- model.frame(formulae[[2]], data=data, na.action='na.pass')
  Wo <- model.matrix(formulae[[2]], data=mf)

  ## compliance model
  mf <- model.frame(formulae[[3]], data=data, na.action='na.fail')
  Xc <- model.matrix(formulae[[3]], data=mf)
  mf <- model.frame(formulae[[4]], data=data, na.action='na.fail')
  Wc <- model.matrix(formulae[[4]], data=mf)
  
  ## response model
  mf <- model.frame(formulae[[5]], data=data, na.action='na.pass')
  Xr <- model.matrix(formulae[[5]], data=mf)
  mf <- model.frame(formulae[[6]], data=data, na.action='na.pass')
  Wr <- model.matrix(formulae[[6]], data=mf)

  ## other variables 
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)
  grp <- eval(call$grp, envir = data)
  ngrp <- length(table(grp))
  if (sum(is.na(Z)) > 0)
    stop("missing values not allowed in the encouragement variable")
  if (sum(is.na(D)) > 0)
    stop("missing values not allowed in the treatment variable")
  if (sum(is.na(grp)) > 0)
    stop("missing values not allowed in the group variable")
  if (model.c == "logit") 
    logit.c <- TRUE
  else
    logit.c <- FALSE
  
  ## Random starting values for missing Y using Bernoulli(0.5)
  R <- (!is.na(Y))*1
  NR <- is.na(Y)
  Ymiss <- sum(NR)
  if (Ymiss > 0) 
    Y[NR] <- (runif(Ymiss) > 0.5)*1
  
  ## Compliance status: 0 = noncomplier, 1 = complier
  C <- rep(NA, N)
  C[Z == 1 &  D == 0] <- 0 # never-takers

  ## Always-takers: 0 = never-takers/compliers, 1 = always-takers
  if (sum(Z == 0 & D == 1)>0) { # some always-takers
    AT <- TRUE
    A <- rep(NA, N)
    if (logit.c)
      C[Z == 0 & D == 1] <- 2 # always-takers
    else
      C[Z == 0 & D == 1] <- 0 # always-takers
    A[Z == 0 & D == 1] <- 1  
    A[Z == 1 & D == 0] <- 0 # never-takers
  } else { # no always-takers
    A <- rep(0, N)
    AT <- FALSE
    C[Z == 1 & D == 1] <- 1 # compliers
  }
  res <- list(call = call, Y = Y, R = R, Xo = Xo, Xc = Xc, Xr = Xr,
              Zo = Wo, Zc = Wc, Zr = Wr, A = A, C = C, D = D, Z = Z,
              grp = grp, n.draws = n.draws)
  
  ## Random starting values for missing compliance status
  if (AT) {
    if (logit.c)
      C[is.na(C) & Z == 1] <- 1 + (runif(sum(is.na(C) & Z == 1)) > 0.5)*1
    else
      C[is.na(C) & Z == 1] <- (runif(sum(is.na(C) & Z == 1)) > 0.5)*1
    C[is.na(C) & Z == 0] <- (runif(sum(is.na(C) & Z == 0)) > 0.5)*1
    if (logit.c)
      A[is.na(A)] <- (C[is.na(A)] == 2)*1
    else {
      A[is.na(A) & Z == 1] <- (C[is.na(A) & Z == 1] == 0)*1
      A[is.na(A)] <- 0
    }
  }
  else
    C[is.na(C)] <- (runif(sum(is.na(C))) > 0.5)*1
  
  ## Completing the outcome model matrix 
  ## The default category is never-takers
  X <- Xo
  W <- Wo
  X1 <- Xr
  W1 <- Wr
  if (AT) { # when some always-takers exist
  ## Xo = [c1 c0 a X] where c1 for compliers with encouragement
  ##                        c0 for compliers without encouragement
  ##                        a for always-takers with/without encouragement  
    Xo <- cbind(0, 0, 0, X)
    Wo <- cbind(0, 0, 0, W)
    Xr <- cbind(0, 0, 0, X1)
    Wr <- cbind(0, 0, 0, W1)
    Xo[A == 1, 3] <- Wo[A == 1, 3] <- 1
    Xr[A == 1, 3] <- Wr[A == 1, 3] <- 1
    colnames(Xo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X))
    colnames(Wo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(W))
    colnames(Xr) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X1))
    colnames(Wr) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(W1))
  } else { # when always-takers do not exist
  ## Xo = [c1 c0 X] where c1 for compliers with encouragement
  ##                      c0 for compliers without encouragement
    Xo <- cbind(0, 0, X)
    Wo <- cbind(0, 0, W)
    Xr <- cbind(0, 0, X1)
    Wr <- cbind(0, 0, W1)
    colnames(Xo) <- c("Complier1", "Complier0", colnames(X))
    colnames(Wo) <- c("Complier1", "Complier0", colnames(W))
    colnames(Xr) <- c("Complier1", "Complier0", colnames(X1))
    colnames(Wr) <- c("Complier1", "Complier0", colnames(W1))
  }
  Xo[C == 1 & Z == 1, 1] <- Wo[C == 1 & Z == 1, 1] <- 1
  Xo[C == 1 & Z == 0, 2] <- Wo[C == 1 & Z == 0, 2] <- 1
  Xr[C == 1 & Z == 1, 1] <- Wr[C == 1 & Z == 1, 1] <- 1
  Xr[C == 1 & Z == 0, 2] <- Wr[C == 1 & Z == 0, 2] <- 1
  
  ## dimensions
  nfixedC <- ncol(Xc)
  nfixedO <- ncol(Xo)
  nfixedR <- ncol(Xr)
  nrandomC <- ncol(Wc)
  nrandomO <- ncol(Wo)
  nrandomR <- ncol(Wr)
  if (AT)
    nqoi <- 8
  else
    nqoi <- 7
  
  ## checking starting values and prior for fixed effects
  if (logit.c & AT) {
    if(length(p.mean.c) != nfixedC*2)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, nfixedC*2)
      else
        stop(paste("the length of p.mean.c should be", nfixedC*2))    
    if(length(coef.start.c) != nfixedC*2)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, nfixedC*2)
      else
        stop(paste("the length of coef.start.c should be", nfixedC*2))        
  } else {
    if(length(p.mean.c) != nfixedC)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, nfixedC)
      else
        stop(paste("the length of p.mean.c should be", nfixedC))        
    if(length(coef.start.c) != nfixedC)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, nfixedC)
      else
        stop(paste("the length of coef.start.c should be", nfixedC))    
  }

  if(length(coef.start.o) != nfixedO)
    if (length(coef.start.o) == 1)
      coef.start.o <- rep(coef.start.o, nfixedO)
    else
      stop(paste("the length of coef.start.o should be", nfixedO))      
  if(length(p.mean.o) != nfixedO)
    if (length(p.mean.o) == 1)
      p.mean.o <- rep(p.mean.o, nfixedO)
    else
      stop(paste("the length of p.mean.o should be", nfixedO))    

  if(length(coef.start.r) != nfixedR)
    if (length(coef.start.r) == 1)
      coef.start.r <- rep(coef.start.r, nfixedR)
    else
      stop(paste("the length of coef.start.r should be", nfixedR))    
  if(length(p.mean.r) != nfixedR)
    if (length(p.mean.r) == 1)
      p.mean.r <- rep(p.mean.r, nfixedR)
    else
      stop(paste("the length of p.mean.r should be", nfixedR))    

  if(is.matrix(p.prec.c)) {
    if (dim(p.prec.c) != rep(nfixedC*2, 2))
        stop(paste("the dimension of p.prec.c should be",
                   rep(nfixedC*2, 2)))    
  } else if (length(p.prec.c) == 1){
    if (logit.c & AT)
      p.prec.c <- diag(p.prec.c, nfixedC*2)
    else
      p.prec.c <- diag(p.prec.c, nfixedC)
  } else {
    stop("Incorrect input for p.prec.c")
  }

  if(is.matrix(p.prec.o)) {
    if (dim(p.prec.o) != rep(nfixedO, 2))
      stop(paste("the dimension of p.prec.o should be",
                 rep(nfixedO, 2)))    
  } else if (length(p.prec.o) == 1){
    p.prec.o <- diag(p.prec.o, nfixedO)
  } else {
    stop("Incorrect input for p.prec.o")
  }

  if(is.matrix(p.prec.r)) {
    if (dim(p.prec.r) != rep(nfixedR, 2))
      stop(paste("the dimension of p.prec.r should be",
                 rep(nfixedR, 2)))    
  } else if (length(p.prec.r) == 1){
    p.prec.r <- diag(p.prec.r, nfixedR)
  } else {
    stop("Incorrect input for p.prec.r")
  }

  ## starting values for Psi
  if(is.matrix(Psi.start.c)) {
    if (dim(Psi.start.c) != rep(nrandomC, 2))
      stop(paste("the dimension of Psi.start.c should be",
                 rep(nrandomC, 2)))    
  } else if (length(Psi.start.c) == 1){
    Psi.start.c <- diag(Psi.start.c, nrandomC)
  } else {
    stop("Incorrect input for Psi.start.c")
  }

  if(is.matrix(Psi.start.o)) {
    if (dim(Psi.start.o) != rep(nrandomO, 2))
      stop(paste("the dimension of Psi.start.o should be",
                 rep(nrandomO, 2)))    
  } else if (length(Psi.start.o) == 1){
    Psi.start.o <- diag(Psi.start.o, nrandomO)
  } else {
    stop("Incorrect input for Psi.start.o")
  }

  if(is.matrix(Psi.start.r)) {
    if (dim(Psi.start.r) != rep(nrandomR, 2))
      stop(paste("the dimension of Psi.start.r should be",
                 rep(nrandomR, 2)))    
  } else if (length(Psi.start.r) == 1){
    Psi.start.r <- diag(Psi.start.r, nrandomR)
  } else {
    stop("Incorrect input for Psi.start.r")
  }
  
  ## starting values for random effects
  xiC <- mvrnorm(ngrp, mu = rep(0, nrandomC),
                 Sigma = solve(Psi.start.c)) 
  xiO <- mvrnorm(ngrp, mu = rep(0, nrandomO),
                 Sigma = solve(Psi.start.r)) 
  xiR <- mvrnorm(ngrp, mu = rep(0, nrandomR),
                 Sigma = solve(Psi.start.r)) 
  
  ## proposal variance for logits
  if (AT) {
    if (length(tune.c) != nfixedC*2)
      if (length(tune.c) == 1)
        tune.c <- rep(tune.c, nfixedC*2)
      else
        stop(paste("the length of tune.c should be", nfixedC*2))
  } else {
    if (length(tune.c) != nfixedC)
      if (length(tune.c) == 1)
        tune.c <- rep(tune.c, nfixedC)
      else
        stop(paste("the length of tune.c should be", nfixedC))
  }
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function
  out <- .C("LIbprobitMixed",
            as.integer(Y), as.integer(R), as.integer(Z),
            as.integer(D), as.integer(C), as.integer(A),
            as.integer(grp), as.integer(Ymiss), as.integer(AT),
            as.integer(in.sample), as.double(Xc), as.double(Wc),
            as.double(Xo), as.double(Wo), as.double(Xr),
            as.double(Wr), as.double(coef.start.c),
            as.double(coef.start.c), as.double(xiC),
            as.double(xiA), as.double(xiO), as.double(xiR),
            as.double(coef.start.o), as.double(coef.start.r),
            as.integer(N), as.integer(n.draws), as.integer(ngrp),
            as.integer(table(grp)), as.integer(max(table(grp))),
            as.integer(nfixedC), as.integer(nfixedO), as.integer(nfixedR),
            as.integer(nrandomC), as.integer(nrandomO), as.integer(nrandomR),
            as.double(Psi.start.c), as.double(Psi.start.c),
            as.double(Psi.start.o), as.double(Psi.start.r),
            as.double(p.mean.c), as.double(p.mean.o),
            as.double(p.mean.r),
            as.double(p.prec.c), as.double(p.prec.o),
            as.double(p.prec.r),
            as.integer(p.df.c), as.integer(p.df.c),
            as.integer(p.df.o), as.integer(p.df.r),
            as.double(p.scale.c), as.double(p.scale.c),
            as.double(p.scale.o), as.double(p.scale.r),
            as.double(tune.c), as.integer(logit.c),
            as.integer(param), as.integer(mda), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
            coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
            coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
            coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
            QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
            PACKAGE = "are")

  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = nfixedC)
    colnames(res$coefC) <- colnames(Xc)
    if (AT) {
      res$coefA <- matrix(out$coefA, byrow = TRUE, ncol = nfixedC)
      colnames(res$coefA) <- colnames(Xc)
    }
    res$coefO <- matrix(out$coefO, byrow = TRUE, ncol = nfixedO)
    colnames(res$coefO) <- colnames(Xo)
    if (Ymiss > 0) {
      res$coefR <- matrix(out$coefR, byrow = TRUE, ncol = nfixedR)
      colnames(res$coefR) <- colnames(Xr)
    }
  }
  QoI <- matrix(out$QoI, byrow = TRUE, ncol = nqoi)
  res$ITT <- QoI[,1]
  res$CACE <- QoI[,2]
  res$pC <- QoI[,3]
  res$pN <- QoI[,4]
  if (AT)
    res$pA <- 1-QoI[,3]-QoI[,4]
  res$Y1barC <- QoI[,5]
  res$Y0barC <- QoI[,6]
  res$YbarN <- QoI[,7]
  if (AT) 
    res$YbarA <- QoI[,8]

  return(res)
}
