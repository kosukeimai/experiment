Noncomp.bayesMixed <- function(formulae, Z, D, grp, data = parent.frame(),
                               n.draws = 5000, param = TRUE, in.sample = FALSE,
                               model.c = "probit", model.o = "probit",
                               random = FALSE, tune.tau = 0.01,
                               tune.var = 0.01, tune.fixed.c = 0.01,
                               tune.random.c = 0.01, tune.fixed.o = 0.01,
                               tune.random.o = 0.01,
                               p.mean.c = 0, p.prec.c = 0.01, p.mean.o = 0,
                               p.prec.o = 0.01, p.mean.r = 0, p.prec.r = 0.01,
                               p.df.var = 10, p.shape.var = 1, p.scale.var = 1,
                               coef.start.c = 0, coef.start.o = 0,
                               tau.start.o = NULL, coef.start.r = 0,
                               var.start.o = 1, Psi.start.c = 1,
                               Psi.start.o = 1, Psi.start.r = 1,
                               p.df.c = 5, p.df.o = 5, p.df.r = 5,
                               p.scale.c = 1, p.scale.o = 1,
                               p.scale.r = 1, burnin = 0, thin = 0,
                               verbose = TRUE) {   
  
  ## models
  if (!(model.o %in% c("probit", "oprobit", "gaussian", "negbin", "twopart")))
    stop("no such model is supported as the outcome model")
  if (!(model.c %in% c("probit", "logit")))
    stop("no such model is supported as the compliance model")
  
  ## getting the data
  call <- match.call()
  ## outcome model: fixed effects
  mf <- model.frame(formulae[[1]], data=data, na.action='na.pass')
  Xo <- model.matrix(formulae[[1]], data=mf)
  if (sum(is.na(Xo)) > 0)
    stop("missing values not allowed in covariates")
  if (model.o %in% c("gaussian", "twopart")) {
    Y <- model.response(mf)
    if (model.o == "twopart") {
      Y1 <- Y
      Y[!is.na(Y)] <- (Y[!is.na(Y)] > 0)*1
    }
  } else if (model.o == "oprobit")
    Y <- as.integer(factor(model.response(mf)))-1
  else
    Y <- as.integer(model.response(mf))

  ## outcome model: random effects
  Wo <- model.matrix(formulae[[2]], data=data)

  ## compliance model
  Xc <- model.matrix(formulae[[3]], data=data)
  Wc <- model.matrix(formulae[[4]], data=data)
  
  ## response model
  Xr <- model.matrix(formulae[[5]], data=data)
  Wr <- model.matrix(formulae[[6]], data=data)

  ## other variables 
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)
  grp <- eval(call$grp, envir = data)
  grp <- as.integer(factor(grp))-1
  ngrp <- length(table(grp))
  if (sum(is.na(Z)) > 0)
    stop("missing values not allowed in the encouragement variable")
  if (sum(is.na(grp)) > 0)
    stop("missing values not allowed in the group variable")

  res <- list(call = call, Y = Y, Xo = Xo, Xc = Xc, Xr = Xr,
              Zo = Wo, Zc = Wc, Zr = Wr, D = D, Z = Z,
              grp = grp, n.draws = n.draws)
  if (model.o == "twopart") {
    res$Y1 <- Y1
    nsamp1 <- sum(Y1[!is.na(Y1)] > 0)
    Y1[is.na(Y1)] <- 0
  }
  
  ## Starting values for missing D
  RD <- (!is.na(D))*1
  NRD <- is.na(D)
  if (sum(NRD) > 0)
    D[NRD] <- Z[NRD]
  
  ## Random starting values for missing Y using Bernoulli(0.5)
  R <- (!is.na(Y))*1
  NR <- is.na(Y)
  Ymiss <- sum(NR)
  if (Ymiss > 0)
    if (model.o == "gaussian")
      Y[NR] <- rnorm(Ymiss)
    else
      Y[NR] <- (runif(Ymiss) > 0.5)*1
  if (model.o == "oprobit") {
    ncat <- max(Y, na.rm = TRUE) + 1
    if (is.null(tau.start.o))
      tau.start.o <- seq(from = 0, length = ncat-1)/10
    if (length(tau.start.o) != (ncat-1))
      stop("incorrect length for tau.start.o")
    if (!identical(sort(tau.start.o), tau.start.o))
      stop("incorrect input for tau.start.o")
    if (length(unique(tau.start.o)) != (ncat-1))
      stop("incorrect input for tau.start.o")
    tau.start.o <- c(tau.start.o, tau.start.o[ncat-1]+1000)
    if (length(tune.tau) != ncat-2)
      if (length(tune.tau) == 1)
        tune.tau <- rep(tune.tau, ncat-2)
      else
        stop(paste("the length of tune.tau should be", ncat-2))
  }
  
  ## Compliance status: 0 = noncomplier, 1 = complier
  C <- rep(NA, N)
  C[Z == 1 &  D == 0] <- 0 # never-takers

  ## Always-takers: 0 = never-takers/compliers, 1 = always-takers
  if (sum(Z == 0 & D == 1)>0) { # some always-takers
    AT <- TRUE
    A <- rep(NA, N)
    A[D == 0] <- 0            # never-takers or compliers
    A[Z == 0 & D == 1] <- 1  
    if (model.c == "logit")
      C[Z == 0 & D == 1] <- 2 # always-takers
    else
      C[Z == 0 & D == 1] <- 0 # always-takers
  } else { # no always-takers
    A <- rep(0, N)
    AT <- FALSE
    C[Z == 1 & D == 1] <- 1 # compliers
  }
  res$R <- R
  res$C <- C
  res$A <- A
  
  ## Random starting values for missing compliance status
  if (AT) {
    A[is.na(A)] <- (runif(sum(is.na(A))) > 0.5)*1
    if (model.c == "logit")
      C[A == 1] <- 2
    else
      C[A == 1] <- 0
  }
  C[is.na(C)] <- (runif(sum(is.na(C))) > 0.5)*1
  
  ## Completing the outcome model matrix 
  ## The default category is never-takers
  X <- Xo
  W <- Wo
  X1 <- Xr
  W1 <- Wr
  if (AT) { # when some always-takers exist
    ## X and W include an intercept
    ## Xo = [c1 c0 a X] where c1 for compliers with encouragement
    ##                        c0 for compliers without encouragement
    ##                        a for always-takers with/without encouragement  
    ## Wo = [c a W] - no Z for parsimony
    ##    = W if random == FALSE
    Xo <- cbind(0, 0, 0, X)
    Xr <- cbind(0, 0, 0, X1)
    Xo[A == 1, 3] <- Xr[A == 1, 3] <- 1
    colnames(Xo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X1))
    if (random) {
      Wo <- cbind(0, 0, W)
      Wr <- cbind(0, 0, W1)
      Wr[A == 1, 2] <- Wo[A == 1, 2] <- 1
      colnames(Wo) <- c("Complier", "AlwaysTaker", colnames(W))
      colnames(Wr) <- c("Complier", "AlwaysTaker", colnames(W1))
    }
  } else { # when always-takers do not exist
    ## Xo = [c1 c0 X] where c1 for compliers with encouragement
    ##                      c0 for compliers without encouragement
    ## Wo = [c W] for parsimony
    ##    = W if random == FALSE
    Xo <- cbind(0, 0, X)
    Xr <- cbind(0, 0, X1)
    colnames(Xo) <- c("Complier1", "Complier0", colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", colnames(X1))
    if (random) {
      Wo <- cbind(0, W)
      Wr <- cbind(0, W1)
      colnames(Wo) <- c("Complier", colnames(W))
      colnames(Wr) <- c("Complier", colnames(W1))
    }
  }
  if (random)
    Wo[C == 1, 1] <- Wr[C == 1, 1] <- 1
  Xo[C == 1 & Z == 1, 1] <- Xr[C == 1 & Z == 1, 1] <- 1
  Xo[C == 1 & Z == 0, 2] <- Xr[C == 1 & Z == 0, 2] <- 1
  
  ## dimensions
  nfixedC <- ncol(Xc)
  nfixedO <- ncol(Xo)
  nfixedR <- ncol(Xr)
  nrandomC <- ncol(Wc)
  nrandomO <- ncol(Wo)
  nrandomR <- ncol(Wr)
  if (model.o == "oprobit")
    if (AT)
      nqoi <- 2 + (ncat-1)*6
    else
      nqoi <- 2 + (ncat-1)*5
  else
    if (AT)
      nqoi <- 8
    else
      nqoi <- 7
  nqoi <- nqoi * (ngrp+1)
  
  ## checking starting values and prior for fixed effects
  if (model.c == "logit" & AT) {
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
    if (model.c == "logit" & AT)
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
  
  ## prior scale matrix for Psi's
  if(is.matrix(p.scale.c)) {
    if (dim(p.scale.c) != rep(nrandomC, 2))
        stop(paste("the dimension of p.scale.c should be",
                   rep(nrandomC, 2)))    
  } else if (length(p.scale.c) == 1)
    p.scale.c <- diag(p.scale.c, nrandomC)
  else 
    stop("Incorrect input for p.scale.c")
  
  if(is.matrix(p.scale.o)) {
    if (dim(p.scale.o) != rep(nrandomO, 2))
        stop(paste("the dimension of p.scale.o should be",
                   rep(nrandomO, 2)))    
  } else if (length(p.scale.o) == 1)
    p.scale.o <- diag(p.scale.o, nrandomO)
  else 
    stop("Incorrect input for p.scale.o")

  if(is.matrix(p.scale.r)) {
    if (dim(p.scale.r) != rep(nrandomR, 2))
        stop(paste("the dimension of p.scale.r should be",
                   rep(nrandomR, 2)))    
  } else if (length(p.scale.r) == 1)
    p.scale.r <- diag(p.scale.r, nrandomR)
  else 
   stop("Incorrect input for p.scale.r")
  
  
  ## proposal variance for logits
  if (AT) {
    if (length(tune.fixed.c) != nfixedC*2)
      if (length(tune.fixed.c) == 1)
        tune.fixed.c <- rep(tune.fixed.c, nfixedC*2)
      else
        stop(paste("the length of tune.fixed.c should be", nfixedC*2))
    if (length(tune.random.c) != 2)
      if (length(tune.random.c) == 1)
        tune.random.c <- rep(tune.random.c, 2)
      else
        stop("the length of tune.random.c should be 2")
  } else {
    if (length(tune.fixed.c) != nfixedC)
      if (length(tune.fixed.c) == 1)
        tune.fixed.c <- rep(tune.fixed.c, nfixedC)
      else
        stop(paste("the length of tune.fixed.c should be", nfixedC))
    if (length(tune.random.c) != 1)
      stop("the length of tune.random.c should be 1")
  }

  ## proposal variance for gaussian and negative binomial
  if (length(tune.var) != 1)
    stop("the length of tune.var should be 1")

  if (model.o == "negbin") {
    if (length(tune.fixed.o) != nfixedO)
      if (length(tune.fixed.o) == 1)
        tune.fixed.o <- rep(tune.fixed.o, nfixedO)
      else
        stop(paste("the length of tune.fixed.o should be", nfixedO))
    if (length(tune.random.o) != nrandomO)
      if (length(tune.random.o) == 1)
        tune.random.o <- rep(tune.random.o, nrandomO)
      else
        stop(paste("the length of tune.random.o should be", nrandomO))
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
  if (model.o == "probit") 
    out <- .C("LIbprobitMixed",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(grp), as.integer(Ymiss), as.integer(AT),
              as.integer(in.sample), as.integer(random),
              as.double(Xc), as.double(Wc), as.double(Xo),
              as.double(Wo), as.double(Xr), as.double(Wr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(coef.start.r),
              as.integer(N), as.integer(n.draws), as.integer(ngrp),
              as.integer(max(table(grp))), as.integer(c(nfixedC,nfixedO,nfixedR)),
              as.integer(c(nrandomC, nrandomO, nrandomR)),
              as.double(Psi.start.c), as.double(Psi.start.c),
              as.double(Psi.start.o), as.double(Psi.start.r),
              as.double(p.mean.c), as.double(p.mean.o), as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o), as.double(p.prec.r),
              as.integer(c(p.df.c,p.df.c,p.df.o,p.df.r)),
              as.double(p.scale.c), as.double(p.scale.c),
              as.double(p.scale.o), as.double(p.scale.r),
              as.double(tune.fixed.c), as.double(tune.random.c),
              as.integer(model.c == "logit"),
              as.integer(param), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
              sPsiC = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiA = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiR = double(nrandomR*(nrandomR+1)*(ceiling((n.draws-burnin)/keep))/2),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  if (model.o == "oprobit")
    out <- .C("LIboprobitMixed",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(grp), as.integer(Ymiss), as.integer(AT),
              as.integer(in.sample), as.integer(random), as.double(Xc), as.double(Wc),
              as.double(Xo), as.double(Wo), as.double(Xr),
              as.double(Wr), as.double(coef.start.c),
              as.double(coef.start.c), as.double(coef.start.o),
              as.double(tau.start.o), as.double(coef.start.r),
              as.integer(N), as.integer(ncat), as.integer(n.draws),
              as.integer(ngrp), as.integer(max(table(grp))),
              as.integer(c(nfixedC,nfixedO,nfixedR)),
              as.integer(c(nrandomC, nrandomO, nrandomR)),
              as.double(Psi.start.c), as.double(Psi.start.c),
              as.double(Psi.start.o), as.double(Psi.start.r),
              as.double(p.mean.c), as.double(p.mean.o), as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o), as.double(p.prec.r),
              as.integer(c(p.df.c,p.df.c,p.df.o,p.df.r)),
              as.double(p.scale.c), as.double(p.scale.c),
              as.double(p.scale.o), as.double(p.scale.r),
              as.double(tune.fixed.c), as.double(tune.random.c),
              as.double(tune.tau), as.integer(TRUE), as.integer(model.c == "logit"),
              as.integer(param), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
              tauO = double((ncat-1)*(ceiling((n.draws-burnin)/keep))),
              sPsiC = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiA = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiR = double(nrandomR*(nrandomR+1)*(ceiling((n.draws-burnin)/keep))/2),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  if (model.o == "negbin")
    out <- .C("LINegBinMixed",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(grp), as.integer(Ymiss), as.integer(AT),
              as.integer(in.sample), as.integer(random), as.double(Xc),
              as.double(Wc), as.double(Xo), as.double(Wo), as.double(Xr),
              as.double(Wr), as.double(coef.start.c),
              as.double(coef.start.c), as.double(coef.start.o),
              as.double(coef.start.r), as.double(var.start.o),
              as.integer(c(N, ngrp)), as.integer(n.draws), 
              as.integer(max(table(grp))), as.integer(c(nfixedC,nfixedO,nfixedR)),
              as.integer(c(nrandomC, nrandomO, nrandomR)),
              as.double(Psi.start.c), as.double(Psi.start.c),
              as.double(Psi.start.o), as.double(Psi.start.r),
              as.double(p.mean.c), as.double(p.mean.o), as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o), as.double(p.prec.r),
              as.double(p.shape.var), as.double(p.scale.var),
              as.integer(c(p.df.c,p.df.c,p.df.o,p.df.r)),
              as.double(p.scale.c), as.double(p.scale.c),
              as.double(p.scale.o), as.double(p.scale.r),
              as.double(tune.fixed.c), as.double(tune.random.c),
              as.double(tune.fixed.o), as.double(tune.random.o),
              as.double(tune.var), as.integer(model.c == "logit.c"),
              as.integer(param), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
              ssig2 = double(ceiling((n.draws-burnin)/keep)),
              sPsiC = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiA = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiR = double(nrandomR*(nrandomR+1)*(ceiling((n.draws-burnin)/keep))/2),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  if (model.o == "gaussian")
    out <- .C("LINormalMixed",
              as.double(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(grp), as.integer(Ymiss), as.integer(AT),
              as.integer(in.sample), as.integer(random), as.double(Xc),
              as.double(Wc),
              as.double(Xo), as.double(Wo), as.double(Xr),
              as.double(Wr), as.double(coef.start.c),
              as.double(coef.start.c), as.double(coef.start.o),
              as.double(coef.start.r), as.double(var.start.o),
              as.integer(N), as.integer(n.draws), as.integer(ngrp),
              as.integer(max(table(grp))),
              as.integer(c(nfixedC,nfixedO,nfixedR)),
              as.integer(c(nrandomC, nrandomO, nrandomR)),
              as.double(Psi.start.c), as.double(Psi.start.c),
              as.double(Psi.start.o), as.double(Psi.start.r),
              as.double(p.mean.c), as.double(p.mean.o), as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.integer(p.df.var), as.double(p.scale.var),
              as.integer(c(p.df.c,p.df.c,p.df.o,p.df.r)),
              as.double(p.scale.c), as.double(p.scale.c),
              as.double(p.scale.o), as.double(p.scale.r),
              as.double(tune.fixed.c), as.double(tune.random.c),
              as.integer(model.c == "logit.c"),
              as.integer(param), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
              ssig2 = double(ceiling((n.draws-burnin)/keep)),
              sPsiC = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiA = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiR = double(nrandomR*(nrandomR+1)*(ceiling((n.draws-burnin)/keep))/2),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  if (model.o == "twopart")
    out <- .C("LItwopartMixed",
              as.integer(Y), as.double(Y1), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(grp), as.integer(c(Ymiss,AT)),
              as.integer(in.sample), as.integer(random), as.double(Xc),
              as.double(Wc), as.double(Xo), as.double(Wo), as.double(Xr),
              as.double(Wr), as.double(coef.start.c),
              as.double(coef.start.c), as.double(coef.start.o), 
              as.double(coef.start.r), as.double(var.start.o),
              as.integer(c(N, nsamp1)), as.integer(n.draws), as.integer(ngrp),
              as.integer(max(table(grp))),
              as.integer(c(nfixedC,nfixedO,nfixedR)),
              as.integer(c(nrandomC, nrandomO, nrandomR)),
              as.double(Psi.start.c), as.double(Psi.start.c),
              as.double(Psi.start.o), as.double(Psi.start.r),
              as.double(p.mean.c), as.double(p.mean.o), as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.integer(p.df.var), as.double(p.scale.var),
              as.integer(c(p.df.c,p.df.c,p.df.o,p.df.r)),
              as.double(p.scale.c), as.double(p.scale.c),
              as.double(p.scale.o), as.double(p.scale.r),
              as.double(tune.fixed.c), as.double(tune.random.c),
              as.integer(model.c == "logit.c"),
              as.integer(param), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(nfixedC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefO1 = double(nfixedO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(nfixedR*(ceiling((n.draws-burnin)/keep))),
              ssig2 = double(ceiling((n.draws-burnin)/keep)),
              sPsiC = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiA = double(nrandomC*(nrandomC+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiO1 = double(nrandomO*(nrandomO+1)*(ceiling((n.draws-burnin)/keep))/2),
              sPsiR = double(nrandomR*(nrandomR+1)*(ceiling((n.draws-burnin)/keep))/2),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")

  
  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = nfixedC)
    colnames(res$coefC) <- colnames(Xc)
    res$PsiC <- matrix(out$sPsiC, byrow = TRUE, ncol = nrandomC*(nrandomC+1)/2)
    if (AT) {
      res$coefA <- matrix(out$coefA, byrow = TRUE, ncol = nfixedC)
      colnames(res$coefA) <- colnames(Xc)
      res$PsiA <- matrix(out$sPsiA, byrow = TRUE, ncol = nrandomC*(nrandomC+1)/2)
    }
    res$coefO <- matrix(out$coefO, byrow = TRUE, ncol = nfixedO)
    res$PsiO <- matrix(out$sPsiO, byrow = TRUE, ncol = nrandomO*(nrandomO+1)/2)
    colnames(res$coefO) <- colnames(Xo)
    if (model.o == "twopart") {
      res$coefO1 <- matrix(out$coefO1, byrow = TRUE, ncol = nfixedO)
      res$PsiO1 <- matrix(out$sPsiO1, byrow = TRUE, ncol = nrandomO*(nrandomO+1)/2)
      colnames(res$coefO1) <- colnames(Xo)
    }
    if (model.o == "oprobit")
      res$tau <- matrix(out$tauO, byrow = TRUE, ncol = ncat - 1)
    if (Ymiss > 0) {
      res$coefR <- matrix(out$coefR, byrow = TRUE, ncol = nfixedR)
      res$PsiR <- matrix(out$sPsiR, byrow = TRUE, ncol = nrandomR*(nrandomR+1)/2)
      colnames(res$coefR) <- colnames(Xr)
    }
    if (model.o == "gaussian")
      res$sig2 <- out$ssig2
  }
  QoI <- matrix(out$QoI, byrow = TRUE, ncol = nqoi)
  if (model.o == "oprobit") {
    res$ITT <- QoI[,1:((ncat-1)*(ngrp+1))]
    res$CACE <- QoI[,((ncat-1)*(ngrp+1)+1):(2*(ngrp+1)*(ncat-1))]
    res$Y1barC <- QoI[,(2*(ngrp+1)*(ncat-1)+1):(3*(ngrp+1)*(ncat-1))]
    res$Y0barC <- QoI[,(3*(ngrp+1)*(ncat-1)+1):(4*(ngrp+1)*(ncat-1))]
    res$YbarN <- QoI[,(4*(ngrp+1)*(ncat-1)+1):(5*(ngrp+1)*(ncat-1))]
    res$pC <- QoI[,(5*(ngrp+1)*(ncat-1)+1):(5*(ngrp+1)*(ncat-1)+ngrp+1)]
    res$pN <- QoI[,(5*(ngrp+1)*(ncat-1)+ngrp+2):(5*(ngrp+1)*(ncat-1)+2*(ngrp+1))]
    if (AT) 
      res$YbarA <- QoI[,(5*(ngrp+1)*(ncat-1)+2*(ngrp+1)+1):(6*(ngrp+1)*(ncat-1)+2*(ngrp+1))]
   } else {
     res$ITT <- QoI[,1:(ngrp+1)]
     res$CACE <- QoI[,((ngrp+1)+1):(2*(ngrp+1))]
     res$pC <- QoI[,(2*(ngrp+1)+1):(3*(ngrp+1))]
     res$pN <- QoI[,(3*(ngrp+1)+1):(4*(ngrp+1))]
     res$Y1barC <- QoI[,(4*(ngrp+1)+1):(5*(ngrp+1))]
     res$Y0barC <- QoI[,(5*(ngrp+1)+1):(6*(ngrp+1))]
     res$YbarN <- QoI[,(6*(ngrp+1)+1):(7*(ngrp+1))]
   }
  if (AT) 
    res$pA <- 1-res$pC-res$pN
  
  return(res)
}
