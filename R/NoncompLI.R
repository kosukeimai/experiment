#' Bayesian Analysis of Randomized Experiments with Noncompliance and Missing
#' Outcomes Under the Assumption of Latent Ignorability
#' 
#' This function estimates the average causal effects for randomized
#' experiments with noncompliance and missing outcomes under the assumption of
#' latent ignorability (Frangakis and Rubin, 1999). The models are based on
#' Bayesian generalized linear models and are fitted using the Markov chain
#' Monte Carlo algorithms. Various types of the outcome variables can be
#' analyzed to estimate the Intention-to-Treat effect and Complier Average
#' Causal Effect.
#' 
#' For the details of the model being fitted, see the references. Note that
#' when always-takers exist we fit either two logistic or two probit models by
#' first modeling whether a unit is a complier or a noncomplier, and then
#' modeling whether a unit is an always-taker or a never-taker for those who
#' are classified as non-compliers.
#' 
#' @param formulae A list of formulae where the first formula specifies the
#' (pre-treatment) covariates in the outcome model (the latent compliance
#' covariate will be added automatically), the second formula specifies the
#' compliance model, and the third formula defines the covariate specification
#' for the model for missing-data mechanism (the latent compliance covariate
#' will be added automatically). For the outcome model, the formula should take
#' the two-sided standard R \command{formula} where the outcome variable is
#' specified in the left hand side of the formula which is then separated by
#' \code{~} from the covariate equation in the right hand side, e.g., \code{y ~
#' x1 + x2}. For the compliance and missing-data mechanism models, the
#' one-sided \command{formula} should be used where the left hand side is left
#' unspecified, e.g., \code{~ x1 + x2}.
#' @param Z A randomized encouragement variable, which should be a binary
#' variable in the specified data frame.
#' @param D A treatment variable, which should be a binary variable in the
#' specified data frame.
#' @param data A data frame which contains the variables that appear in the
#' model formulae (\code{formulae}), the encouragement variable (\code{Z}), and
#' the treatment variable (\code{D}).
#' @param n.draws The number of MCMC draws. The default is \code{5000}.
#' @param param A logical variable indicating whether the Monte Carlo draws of
#' the model parameters should be saved in the output object. The default is
#' \code{TRUE}.
#' @param in.sample A logical variable indicating whether or not the sample
#' average causal effect should be calculated using the observed potential
#' outcome for each unit. If it is set to \code{FALSE}, then the population
#' average causal effect will be calculated. The default is \code{FALSE}.
#' @param model.c The model for compliance. Either \code{logit} or
#' \code{probit} model is allowed. The default is \code{probit}.
#' @param model.o The model for outcome. The following five models are allowed:
#' \code{logit}, \code{probit}, \code{oprobit} (ordered probit regression),
#' \code{gaussian} (gaussian regression), \code{negbin} (negative binomial
#' regression), and \code{twopart} (two part model where the first part is the
#' probit regression for \eqn{Pr(Y>0|X)} and the second part models
#' \eqn{p(log(Y)|X, Y>0)} using the gaussian regression). The default is
#' \code{probit}.
#' @param model.r The model for (non)response. Either \code{logit} or
#' \code{probit} model is allowed. The default is \code{probit}.
#' @param tune.c Tuning constants for fitting the compliance model. These
#' positive constants are used to tune the (random-walk) Metropolis-Hastings
#' algorithm to fit the logit model. Use either a scalar or a vector of
#' constants whose length equals that of the coefficient vector. The default is
#' \code{0.01}.
#' @param tune.o Tuning constants for fitting the outcome model. These positive
#' constants are used to tune the (random-walk) Metropolis-Hastings algorithm
#' to fit logit, ordered probit, and negative binomial models. Use either a
#' scalar or a vector of constants whose length equals that of the coefficient
#' vector for logit and negative binomial models. For the ordered probit model,
#' use either a scalar or a vector of constants whose length equals that of
#' cut-point parameters to be estimated. The default is \code{0.01}.
#' @param tune.r Tuning constants for fitting the (non)response model. These
#' positive constants are used to tune the (random-walk) Metropolis-Hastings
#' algorithm to fit the logit model. Use either a scalar or a vector of
#' constants whose length equals that of the coefficient vector. The default is
#' \code{0.01}.
#' @param tune.v A scalar tuning constant for fitting the variance component of
#' the negative binomial (outcome) model. The default is \code{0.01}.
#' 
#' @param p.mean.c Prior mean for the compliance model. It should be either a
#' scalar or a vector of appropriate length. The default is \code{0}.
#' @param p.prec.c Prior precision for the compliance model. It should be
#' either a positive scalar or a positive semi-definite matrix of appropriate
#' size. The default is \code{0.001}.
#' @param p.mean.o Prior mean for the outcome model. It should be either a
#' scalar or a vector of appropriate length. The default is \code{0}.
#' @param p.prec.o Prior precision for the outcome model. It should be either a
#' positive scalar or a positive semi-definite matrix of appropriate size. The
#' default is \code{0.001}.
#' @param p.mean.r Prior mean for the (non)response model. It should be either
#' a scalar or a vector of appropriate length. The default is \code{0}.
#' @param p.prec.r Prior precision for the (non)response model. It should be
#' either a positive scalar or a positive semi-definite matrix of appropriate
#' size. The default is \code{0.001}.
#' @param p.df.o A positive integer. Prior degrees of freedom parameter for the
#' inverse chisquare distribution in the gaussian and twopart (outcome) models.
#' The default is \code{10}.
#' @param p.scale.o A positive scalar. Prior scale parameter for the inverse
#' chisquare distribution (for the variance) in the gaussian and twopart
#' (outcome) models. For the negative binomial (outcome) model, this is used
#' for the scale parameter of the inverse gamma distribution. The default is
#' \code{1}.
#' @param p.shape.o A positive scalar. Prior shape for the inverse chisquare
#' distribution in the negative binomial (outcome) model. The default is
#' \code{1}.
#' @param mda.probit A logical variable indicating whether to use marginal data
#' augmentation for probit models. The default is \code{TRUE}.
#' @param coef.start.c Starting values for coefficients of the compliance
#' model.  It should be either a scalar or a vector of appropriate length. The
#' default is \code{0}.
#' @param coef.start.o Starting values for coefficients of the outcome model.
#' It should be either a scalar or a vector of appropriate length. The default
#' is \code{0}.
#' @param coef.start.r Starting values for coefficients of the (non)response
#' model.  It should be either a scalar or a vector of appropriate length. The
#' default is \code{0}.
#' @param tau.start.o Starting values for thresholds of the ordered probit
#' (outcome) model.  If it is set to \code{NULL}, then the starting values will
#' be a sequence starting from 0 and then incrementing by 0.1. The default is
#' \code{NULL}.
#' @param var.start.o A positive scalar starting value for the variance of the
#' gaussian, negative binomial, and twopart (outcome) models. The default is
#' \code{1}.
#' @param burnin The number of initial burnins for the Markov chain. The
#' default is \code{0}.
#' @param thin The size of thinning interval for the Markov chain. The default
#' is \code{0}.
#' @param verbose A logical variable indicating whether additional progress
#' reports should be prited while running the code. The default is \code{TRUE}.
#' @return An object of class \code{NoncompLI} which contains the following
#' elements as a list: \item{call}{The matched call.} \item{Y}{The outcome
#' variable.} \item{D}{The treatment variable.} \item{Z}{The (randomized)
#' encouragement variable.} \item{R}{The response indicator variable for
#' \code{Y}.} \item{A}{The indicator variable for (known) always-takers, i.e.,
#' the control units who received the treatment.} \item{C}{The indicator
#' variable for (known) compliers, i.e., the encouraged units who received the
#' treatment when there is no always-takers.} \item{Xo}{The matrix of
#' covariates used for the outcome model.} \item{Xc}{The matrix of covariates
#' used for the compliance model.} \item{Xr}{The matrix of covariates used for
#' the (non)response model.} \item{n.draws}{The number of MCMC draws.}
#' \item{QoI}{The Monte carlo draws of quantities of interest from their
#' posterior distributions. Quantities of interest include \code{ITT}
#' (intention-to-treat) effect, \code{CACE} (complier average causal effect),
#' \code{Y1barC} (The mean outcome value under the treatment for compliers),
#' \code{Y0barC} (The mean outcome value under the control for compliers),
#' \code{YbarN} (The mean outcome value for never-takers), \code{YbarA} (The
#' mean outcome value for always-takers), \code{pC} (The proportion of
#' compliers), \code{pN} (The proportion of never-takers), \code{pA} (The
#' proportion of always-takers) } If \code{param} is set to \code{TRUE}, the
#' following elments are also included: \item{coefO}{The Monte carlo draws of
#' coefficients of the outcome model from their posterior distribution.}
#' \item{coefO1}{If \code{model = "twopart"}, this element contains the Monte
#' carlo draws of coefficients of the outcome model for \eqn{p(log(Y)|X, Y >
#' 0)} from their posterior distribution.} \item{coefC}{The Monte carlo draws
#' of coefficients of the compliance model from their posterior distribution.}
#' \item{coefA}{If always-takers exist, then this element contains the Monte
#' carlo draws of coefficients of the compliance model for always-takers from
#' their posterior distribution.} \item{coefR}{The Monte carlo draws of
#' coefficients of the (non)response model from their posterior distribution.}
#' \item{sig2}{The Monte carlo draws of the variance parameter for the
#' gaussian, negative binomial, and twopart (outcome) models.}
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' @references Frangakis, Constantine E. and Donald B. Rubin. (1999).
#' \dQuote{Addressing Complications of Intention-to-Treat Analysis in the
#' Combined Presence of All-or-None Treatment Noncompliance and Subsequent
#' Missing Outcomes.} \emph{Biometrika}, Vol. 86, No. 2, pp. 365-379.
#' 
#' Hirano, Keisuke, Guido W. Imbens, Donald B. Rubin, and Xiao-Hua Zhou.
#' (2000).  \dQuote{Assessing the Effect of an Influenza Vaccine in an
#' Encouragement Design.} \emph{Biostatistics}, Vol. 1, No. 1, pp. 69-88.
#' 
#' Barnard, John, Constantine E. Frangakis, Jennifer L. Hill, and Donald B.
#' Rubin. (2003).  \dQuote{Principal Stratification Approach to Broken
#' Randomized Experiments: A Case Study of School Choice Vouchers in New York
#' (with Discussion)}, \emph{Journal of the American Statistical Association},
#' Vol. 98, No. 462, pp299--311.
#' 
#' Horiuchi, Yusaku, Kosuke Imai, and Naoko Taniguchi (2007). \dQuote{Designing
#' and Analyzing Randomized Experiments: Application to a Japanese Election
#' Survey Experiment.} \emph{American Journal of Political Science}, Vol. 51,
#' No. 3 (July), pp. 669-687.
#' @keywords models
NoncompLI <- function(formulae, Z, D, data = parent.frame(), n.draws = 5000,
                      param = TRUE, in.sample = FALSE, model.c = "probit",
                      model.o = "probit", model.r = "probit", 
                      tune.c = 0.01, tune.o = 0.01, tune.r = 0.01,
                      tune.v = 0.01, p.mean.c = 0, p.mean.o = 0,
                      p.mean.r = 0, p.prec.c = 0.001,
                      p.prec.o = 0.001, p.prec.r = 0.001,
                      p.df.o = 10, p.scale.o = 1, p.shape.o = 1,
                      mda.probit = TRUE, coef.start.c = 0,
                      coef.start.o = 0, tau.start.o = NULL,
                      coef.start.r = 0, var.start.o = 1,
                      burnin = 0, thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()

  ## model types
  if (!(model.c %in% c("logit", "probit"))) 
    stop("no such model is supported for the compliance model.")
  
  if (!(model.o %in% c("logit", "probit", "oprobit", "gaussian",
                       "negbin", "twopart"))) 
    stop("no such model is supported for the outcome model.")
  
  if (!(model.r %in% c("logit", "probit"))) 
    stop("no such model is supported for the response model.")    

  ## outcome model
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
    Y <- as.integer(factor(model.response(mf)))
  else
    Y <- as.integer(model.response(mf))

  ## compliance model
  Xc <- model.matrix(formulae[[2]], data=data)

  ## response model
  if (any(is.na(Y)))
    Xr <- model.matrix(formulae[[3]], data=data)
  else
    Xr <- model.matrix(~ 1, data = data)
  
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)
  if (sum(is.na(Z)) > 0)
    stop("missing values not allowed in the encouragement variable")

  res <- list(call = call, Y = Y, Xo = Xo, Xc = Xc, Xr = Xr,
              D = D, Z = Z, n.draws = n.draws)
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

  ## Random starting values for missing Y using Bernoulli(0.5) for
  ## binary and ordinal
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
    C[Z == 1 & D == 1] <- 1   # compliers
  }
  res$R <- R
  res$A <- A
  res$C <- C
  
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
  X1 <- Xr
  if (AT) { # when some always-takers exist
  ## Xo = [c1 c0 a X] where c1 for compliers with encouragement
  ##                        c0 for compliers without encouragement
  ##                        a for always-takers with/without encouragement  
    Xo <- cbind(0, 0, 0, X)
    Xr <- cbind(0, 0, 0, X1)
    Xo[A == 1, 3] <- 1
    Xr[A == 1, 3] <- 1
    colnames(Xo) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", "AlwaysTaker",
                      colnames(X1))
  } else { # when always-takers do not exist
  ## Xo = [c1 c0 X] where c1 for compliers with encouragement
  ##                      c0 for compliers without encouragement
    Xo <- cbind(0, 0, X)
    Xr <- cbind(0, 0, X1)
    colnames(Xo) <- c("Complier1", "Complier0", colnames(X))
    colnames(Xr) <- c("Complier1", "Complier0", colnames(X1))
  }
  Xo[C == 1 & Z == 1, 1] <- 1
  Xo[C == 1 & Z == 0, 2] <- 1
  Xr[C == 1 & Z == 1, 1] <- 1
  Xr[C == 1 & Z == 0, 2] <- 1
  
  ## dimensions
  ncovC <- ncol(Xc)
  ncovO <- ncol(Xo)
  ncovR <- ncol(Xr)
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

  ## checking starting values and prior
  if ((model.c == "logit") & AT) {
    if(length(p.mean.c) != ncovC*2)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, ncovC*2)
      else
        stop(paste("the length of p.mean.c should be", ncovC*2))    
    if(length(coef.start.c) != ncovC*2)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, ncovC*2)
      else
        stop(paste("the length of coef.start.c should be", ncovC*2))        
  } else {
    if(length(p.mean.c) != ncovC)
      if (length(p.mean.c) == 1)
        p.mean.c <- rep(p.mean.c, ncovC)
      else
        stop(paste("the length of p.mean.c should be", ncovC))        
    if(length(coef.start.c) != ncovC)
      if (length(coef.start.c) == 1)
        coef.start.c <- rep(coef.start.c, ncovC)
      else
        stop(paste("the length of coef.start.c should be", ncovC))    
  }

  if(length(coef.start.o) != ncovO)
    if (length(coef.start.o) == 1)
      coef.start.o <- rep(coef.start.o, ncovO)
    else
      stop(paste("the length of coef.start.o should be", ncovO))      
  if(length(p.mean.o) != ncovO)
    if (length(p.mean.o) == 1)
      p.mean.o <- rep(p.mean.o, ncovO)
    else
      stop(paste("the length of p.mean.o should be", ncovO))    

  if(length(coef.start.r) != ncovR)
    if (length(coef.start.r) == 1)
      coef.start.r <- rep(coef.start.r, ncovR)
    else
      stop(paste("the length of coef.start.r should be", ncovR))    
  if(length(p.mean.r) != ncovR)
    if (length(p.mean.r) == 1)
      p.mean.r <- rep(p.mean.r, ncovR)
    else
      stop(paste("the length of p.mean.r should be", ncovR))    

  if(is.matrix(p.prec.c)) {
    if (dim(p.prec.c) != rep(ncovC*2, 2))
        stop(paste("the dimension of p.prec.c should be",
                   rep(ncovC*2, 2)))    
  } else if (length(p.prec.c) == 1){
    if (model.c == "logit" & AT)
      p.prec.c <- diag(p.prec.c, ncovC*2)
    else
      p.prec.c <- diag(p.prec.c, ncovC)
  } else {
    stop("Incorrect input for p.prec.c")
  }

  if(is.matrix(p.prec.o)) {
    if (dim(p.prec.o) != rep(ncovO, 2))
      stop(paste("the dimension of p.prec.o should be",
                 rep(ncovO, 2)))    
  } else if (length(p.prec.o) == 1){
    p.prec.o <- diag(p.prec.o, ncovO)
  } else {
    stop("Incorrect input for p.prec.o")
  }

  if(is.matrix(p.prec.r)) {
    if (dim(p.prec.r) != rep(ncovR, 2))
      stop(paste("the dimension of p.prec.r should be",
                 rep(ncovR, 2)))    
  } else if (length(p.prec.r) == 1){
    p.prec.r <- diag(p.prec.r, ncovR)
  } else {
    stop("Incorrect input for p.prec.r")
  }
  
  ## proposal variance for logits
  if (model.c == "logit")
    if (AT) {
      if (length(tune.c) != ncovC*2)
        if (length(tune.c) == 1)
          tune.c <- rep(tune.c, ncovC*2)
        else
          stop(paste("the length of tune.c should be", ncovC*2))
    } else {
      if (length(tune.c) != ncovC)
        if (length(tune.c) == 1)
          tune.c <- rep(tune.c, ncovC)
        else
          stop(paste("the length of tune.c should be", ncovC))
    }
  if (model.o == "logit" || model.o == "negbin")
    if (length(tune.o) != ncovO)
      if (length(tune.o) == 1)
        tune.o <- rep(tune.o, ncovO)
      else
        stop(paste("the length of tune.o should be", ncovO))
  if (model.o == "oprobit")
    if (length(tune.o) != ncat-2)
      if (length(tune.o) == 1)
        tune.o <- rep(tune.o, ncat-2)
      else
        stop(paste("the length of tune.o should be", ncat-2))
  if (model.r == "logit")
    if (length(tune.r) != ncovR)
      if (length(tune.r) == 1)
        tune.r <- rep(tune.r, ncovR)
      else
        stop(paste("the length of tune.r should be", ncovR))
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1

  ## calling C function
  if (model.o == "probit" || model.o == "logit")
    out <- .C("LIbinary",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.double(tune.c), as.double(tune.o), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.o == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "oprobit")
    out <- .C("LIordinal",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(tau.start.o),
              as.double(coef.start.r), 
              as.integer(N), as.integer(n.draws), as.integer(ncat),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r),
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r),
              as.double(tune.c), as.double(tune.o), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              tauO = double((ncat-1)*(ceiling((n.draws-burnin)/keep))),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "gaussian")
    out <- .C("LIgaussian",
              as.double(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(var.start.o),
              as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.integer(p.df.o),
              as.double(p.scale.o),
              as.double(tune.c), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "negbin")
    out <- .C("LIcount",
              as.integer(Y), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(var.start.o),
              as.double(coef.start.r),
              as.integer(N), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.double(p.shape.o),
              as.double(p.scale.o), as.double(tune.c),
              as.double(tune.r), as.double(tune.o), as.double(tune.v),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  else if (model.o == "twopart")
    out <- .C("LItwopart",
              as.integer(Y), as.double(Y1), as.integer(R), as.integer(Z),
              as.integer(D), as.integer(RD), as.integer(C), as.integer(A),
              as.integer(Ymiss), as.integer(AT), as.integer(in.sample), 
              as.double(Xc), as.double(Xo), as.double(Xr),
              as.double(coef.start.c), as.double(coef.start.c),
              as.double(coef.start.o), as.double(coef.start.o),
              as.double(var.start.o), as.double(coef.start.r),
              as.integer(c(N, nsamp1)), as.integer(n.draws),
              as.integer(ncovC), as.integer(ncovO), as.integer(ncovR),
              as.double(p.mean.c), as.double(p.mean.o),
              as.double(p.mean.r), 
              as.double(p.prec.c), as.double(p.prec.o),
              as.double(p.prec.r), as.integer(p.df.o),
              as.double(p.scale.o),
              as.double(tune.c), as.double(tune.r),
              as.integer(model.c == "logit"),
              as.integer(model.r == "logit"),
              as.integer(param), as.integer(mda.probit), as.integer(burnin),
              as.integer(keep), as.integer(verbose),
              coefC = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefA = double(ncovC*(ceiling((n.draws-burnin)/keep))),
              coefO = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefO1 = double(ncovO*(ceiling((n.draws-burnin)/keep))),
              coefR = double(ncovR*(ceiling((n.draws-burnin)/keep))),
              var = double(ceiling((n.draws-burnin)/keep)),
              QoI = double(nqoi*(ceiling((n.draws-burnin)/keep))),
              PACKAGE = "experiment")
  
  if (param) {
    res$coefC <- matrix(out$coefC, byrow = TRUE, ncol = ncovC)
    colnames(res$coefC) <- colnames(Xc)
    if (AT) {
      res$coefA <- matrix(out$coefA, byrow = TRUE, ncol = ncovC)
      colnames(res$coefA) <- colnames(Xc)
    }
    res$coefO <- matrix(out$coefO, byrow = TRUE, ncol = ncovO)
    colnames(res$coefO) <- colnames(Xo)
    if (model.o == "twopart") {
      res$coefO1 <- matrix(out$coefO1, byrow = TRUE, ncol = ncovO)
      colnames(res$coefO1) <- colnames(Xo)
    }
    if (model.o == "oprobit")
      res$tau <- matrix(out$tauO, byrow = TRUE, ncol = ncat - 1)
    if (Ymiss > 0) {
      res$coefR <- matrix(out$coefR, byrow = TRUE, ncol = ncovR)
      colnames(res$coefR) <- colnames(Xr)
    }
    if (model.o %in% c("gaussian", "negbin", "twopart"))
      res$sig2 <- out$var
  }

  QoI <- matrix(out$QoI, byrow = TRUE, ncol = nqoi)
  if (model.o == "oprobit") {
    res$ITT <- QoI[,1:(ncat-1)]
    res$CACE <- QoI[,ncat:(2*(ncat-1))]
    res$Y1barC <- QoI[,(2*(ncat-1)+1):(3*(ncat-1))]
    res$Y0barC <- QoI[,(3*(ncat-1)+1):(4*(ncat-1))]
    res$YbarN <- QoI[,(4*(ncat-1)+1):(5*(ncat-1))]
    res$pC <- QoI[,(5*(ncat-1)+1)]
    res$pN <- QoI[,(5*(ncat-1)+2)]
    if (AT) 
      res$YbarA <- QoI[,(5*(ncat-1)+3):(6*(ncat-1)+2)]
  } else {
    res$ITT <- QoI[,1]
    res$CACE <- QoI[,2]
    res$pC <- QoI[,3]
    res$pN <- QoI[,4]
    res$Y1barC <- QoI[,5]
    res$Y0barC <- QoI[,6]
    res$YbarN <- QoI[,7]
    if (AT) 
      res$YbarA <- QoI[,8]
  }
  if (AT) 
    res$pA <- 1-res$pC-res$pN

  class(res) <- "NoncompLI"
  return(res)
}
