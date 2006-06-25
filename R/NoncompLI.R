Noncomp.bprobit <- function(formulae, Z, D, data = parent.frame(),
                            n.draws = 5000, param = TRUE,
                            model = "probit", propVar = 1,
                            p.mean.c = 0, p.var.c = 1000, p.mean.o = 0,
                            p.var.o = 1000, p.mean.r = 0, p.var.r = 1000,
                            mda = TRUE, coef.start.c = 0,
                            coef.start.o = 0, coef.start.r = 0, burnin = 0,
                            thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()
  ## outcome model
  mf <- model.frame(formulae[[1]], data=data, na.action='na.pass')
  Xo <- model.matrix(formulae[[1]], data=mf)
  if (sum(is.na(Xo)) > 0)
    stop("missing values not allowed in covariates")
  Y <- as.integer(model.response(mf))
  ## compliance model
  mf <- model.frame(formulae[[2]], data=data, na.action='na.fail')
  Xc <- model.matrix(formulae[[2]], data=mf)
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)
  if (sum(is.na(Z)) > 0)
    stop("missing values not allowed in the encouragement variable")
  if (sum(is.na(D)) > 0)
    stop("missing values not allowed in the treatment variable")
  if (model == "logit") 
    logit <- TRUE
  else
    logit <- FALSE
  
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
    if (logit)
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
  res <- list(call = call, Y = Y, R = R, Xo = Xo, Xc = Xc, A = A, C = C,
              D = D, Z = Z, n.draws = n.draws)
  
  ## Random starting values for missing compliance status
  if (AT) {
    if (logit)
      C[is.na(C) & Z == 1] <- 1 + (runif(sum(is.na(C) & Z == 1)) > 0.5)*1
    else
      C[is.na(C) & Z == 1] <- (runif(sum(is.na(C) & Z == 1)) > 0.5)*1
    C[is.na(C) & Z == 0] <- (runif(sum(is.na(C) & Z == 0)) > 0.5)*1
    if (logit)
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
  if (AT) { # when some always-takers exist
  ## Xo = [c1 c0 a X] where c1 for compliers with encouragement
  ##                        c0 for compliers without encouragement
  ##                        a for always-takers with/without encouragement  
    Xo <- cbind(0, 0, 0, X)
    Xo[A == 1, 3] <- 1
    colnames(Xo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X))
  } else { # when always-takers do not exist
  ## Xo = [c1 c0 X] where c1 for compliers with encouragement
  ##                      c0 for compliers without encouragement
    Xo <- cbind(0, 0, X)
    colnames(Xo) <- c("Complier1", "Complier0", colnames(X))
  }
  Xo[C == 1 & Z == 1, 1] <- 1
  Xo[C == 1 & Z == 0, 2] <- 1
  
  ## dimensions
  ncovC <- ncol(Xc)
  ncovO <- ncol(Xo)
  if (AT)
    nqoi <- 8
  else
    nqoi <- 7
  
  ## starting values and prior
  if (logit & AT) {
    if(length(p.mean.c) != ncovC*2)
      p.mean.c <- rep(p.mean.c, ncovC*2)
    if(length(coef.start.c) != ncovC*2)
      coef.start.c <- rep(coef.start.c, ncovC*2)
  } else {
    if(length(p.mean.c) != ncovC)
      p.mean.c <- rep(p.mean.c, ncovC)
    if(length(coef.start.c) != ncovC)
      coef.start.c <- rep(coef.start.c, ncovC)
  }

  if(length(coef.start.o) != ncovO)
    coef.start.o <- rep(coef.start.o, ncovO)
  if(length(p.mean.o) != ncovO)
    p.mean.o <- rep(p.mean.o, ncovO)

  if(length(coef.start.r) != ncovO)
    coef.start.r <- rep(coef.start.r, ncovO)
  if(length(p.mean.r) != ncovO)
    p.mean.r <- rep(p.mean.r, ncovO)

  if(!is.matrix(p.var.c)) 
    if (logit & AT)
      p.var.c <- diag(p.var.c, ncovC*2)
    else
      p.var.c <- diag(p.var.c, ncovC)
  if(!is.matrix(p.var.o))
    p.var.o <- diag(p.var.o, ncovO)
  if(!is.matrix(p.var.r))
    p.var.r <- diag(p.var.r, ncovO)

  ## proposal variance for logit
  if (AT) {
    if (length(propVar) != ncovC*2 )
      propVar <- rep(propVar, ncovC*2)
  } else {
    if (length(propVar) != ncovC )
      propVar <- rep(propVar, ncovC)
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
  out <- .C("LIbprobit",
            as.integer(Y), as.integer(R), as.integer(Z),
            as.integer(D), as.integer(C), as.integer(A),
            as.integer(Ymiss), as.integer(AT),
            as.double(Xc), as.double(Xo),
            as.double(coef.start.c), as.double(coef.start.c),
            as.double(coef.start.o), as.double(coef.start.r),
            as.integer(N), as.integer(n.draws),
            as.integer(ncovC), as.integer(ncovO),
            as.double(p.mean.c), as.double(p.mean.o),
            as.double(p.mean.r),
            as.double(solve(p.var.c)), as.double(solve(p.var.o)),
            as.double(solve(p.var.r)),
            as.double(propVar), as.integer(logit),
            as.integer(param), as.integer(mda), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
            coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
            coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
            coefR = double(ncovO*(ceiling((n.draws-burnin)/keep))),
            QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="are")

  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = ncovC)
    colnames(res$coefC) <- colnames(Xc)
    if (AT) {
      res$coefA <- matrix(out$coefA, byrow = TRUE, ncol = ncovC)
      colnames(res$coefA) <- colnames(Xc)
    }
    res$coefO <- matrix(out$coefO, byrow = TRUE, ncol = ncovO)
    colnames(res$coefO) <- colnames(Xo)
    if (Ymiss > 0) {
      res$coefR <- matrix(out$coefR, byrow = TRUE, ncol = ncovO)
      colnames(res$coefR) <- colnames(Xo)
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
  print(logit)
  if (AT) 
    res$YbarA <- QoI[,8]

  return(res)
}
