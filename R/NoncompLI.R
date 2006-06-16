Noncomp.bprobit <- function(formulae, Z, D, data = parent.frame(),
                            n.draws = 5000, insample = TRUE, param = TRUE,
                            p.mean.c = 0, p.var.c = 100, p.mean.o = 0,
                            p.var.o = 100, mda = TRUE, coef.start.c = 0,
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
  ncovc <- ncol(Xc)
  ncovo <- ncol(Xo)
    
  ## starting values
  if(length(coef.start.c) != ncovc)
    coef.start.c <- rep(coef.start.c, ncovc)
  if(length(coef.start.o) != ncovo)
    coef.start.o <- rep(coef.start.o, ncovo)
 
  ## prior
  if(length(p.mean.c) != ncovc)
    p.mean.c <- rep(p.mean.c, ncovc)
  if(length(p.mean.o) != ncovo)
    p.mean.o <- rep(p.mean.o, ncovo)
  if(!is.matrix(p.var.c))
    p.var.c <- diag(p.var.c, ncovc)
  if(!is.matrix(p.var.o))
    p.var.o <- diag(p.var.o, ncovo)

  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1
  
  ## calling C function
  out <- .C("probit",
            as.integer(Y), as.integer(R), as.integer(Ymax),
            as.integer(Z), as.integer(D), as.integer(C), 
            as.double(X), as.double(Xo),
            as.double(coef.start.c), as.double(coef.start.o),
            as.integer(N), as.integer(n.draws),
            as.integer(ncov), as.integer(ncovo), as.integer(ncovX),
            as.integer(N11),
            as.double(p.mean.c), as.double(p.mean.o),
            as.double(solve(p.var.c)), as.double(solve(p.var.o)),
            as.integer(insample), as.integer(varT),
            as.integer(param), as.integer(mda), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            pdStore = double(allpar*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="are")

  return(res)
}
