Noncomp.bprobit <- function(formulae, Z, D, data = parent.frame(),
                            n.draws = 5000, param = TRUE,
                            p.mean.c = 0, p.var.c = 1000, p.mean.o = 0,
                            p.var.o = 1000, mda = TRUE, coef.start.c = 0,
                            coef.start.o = 0, burnin = 0,
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
  
  ## Random starting values for missing Y using Bernoulli(0.5)
  R <- is.na(Y)
  Y[R] <- (runif(sum(R)) > 0.5)*1
  
  ## Compliance status: 0 = noncomplier, 1 = complier
  C <- rep(NA, N)
  C[Z == 1 &  D == 0] <- 0 # never-takers

  ## Always-takers: 0 = never-takers/compliers, 1 = always-takers
  if (sum(Z == 0 & D == 1)>0) { # some always-takers
    AT <- TRUE
    A <- rep(NA, N)
    C[Z == 0 & D == 1] <- 0 # always-takers
    A[Z == 0 & D == 1] <- 1  
    A[Z == 1 & D == 0] <- 0 # never-takers
  } else { # no always-takers
    A <- rep(0, N)
    C[Z == 1 & D == 1] <- 1 # compliers
  }
  res <- list(call = call, Y = Y, Xo = Xo, Xc = Xc, A = A, C = C,
              D = D, Z = Z, n.draws = n.draws)
  
  ## Random starting values for missing compliance status
  C[is.na(C)] <- (runif(sum(is.na(C))) > 0.5)*1
  if (AT)
    A[is.na(A)] <- (runif(sum(is.na(A))) > 0.5)*1

  ## Completing the outcome model matrix 
  ## The default category is never-takers
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
    
  ## starting values
  if(length(coef.start.c) != ncovC)
    coef.start.c <- rep(coef.start.c, ncovC)
  if(length(coef.start.o) != ncovO)
    coef.start.o <- rep(coef.start.o, ncovO)
 
  ## prior
  if(length(p.mean.c) != ncovC)
    p.mean.c <- rep(p.mean.c, ncovC)
  if(length(p.mean.o) != ncovO)
    p.mean.o <- rep(p.mean.o, ncovO)
  if(!is.matrix(p.var.c))
    p.var.c <- diag(p.var.c, ncovC)
  if(!is.matrix(p.var.o))
    p.var.o <- diag(p.var.o, ncovO)

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
            as.interger(AT), as.double(Xc), as.double(Xo),
            as.double(coef.start.c), as.double(coef.start.o),
            as.integer(N), as.integer(n.draws),
            as.integer(ncovC), as.integer(ncovO),
            as.double(p.mean.c), as.double(p.mean.o),
            as.double(solve(p.var.c)), as.double(solve(p.var.o)),
            as.integer(param), as.integer(mda), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
            coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
            coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
            QoI = if(AT) double(8*(ceiling((n.draws-burnin)/keep)))
            else double(7*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="are")

  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = ncovC)
    colnames(res$coefC) <- colnames(Xc)
    if (AT) {
      res$coefA <- matrix(out$coefC, byrow = TRUE, ncol = ncovC)
      colnames(res$coefA) <- colnames(Xc)
    }
    res$coefO <- matrix(out$coefC, byrow = TRUE, ncol = ncovO)
    colnames(res$coefO) <- colnames(Xo)
  }
  QoI <- matrix(out$QoI, byrow = TRUE, ncol = if (AT) 8 else 7)
  res$ITT <- QoI[,1]
  res$CACE <- QoI[,2]
  res$pC <- QoI[,3]
  res$pN <- QoI[,4]
  if (AT)
    res$pA <- 1-QoI[,3]-QoI[,4]
  res$Y1barC <- res$QoI[,5]
  res$Y0barC <- res$QoI[,6]
  res$YbarN <- res$QoI[,7]
  if (AT)
    res$YbarA <- res$QoI[,8]

  return(res)
}
