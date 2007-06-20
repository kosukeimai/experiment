\name{NoncompLI}

\alias{NoncompLI}

\title{Bayesian Analysis of Randomized Experiments with Noncompliance
  and Missing Outcomes Under the Assumption of Latent Ignorability}

\description{
  This function estimates the average causal effects for randomized
  experiments with noncompliance and missing outcomes under the
  assumption of latent ignorability (Frangakis and Rubin, 1999). The
  models are based on Bayesian generalized linear models and are fitted
  using the Markov chain Monte Carlo algorithms. Various types of the
  outcome variables can be analyzed to estimate the Intention-to-Treat
  effect and Complier Average Causal Effect.
}

\usage{
NoncompLI(formulae, Z, D, data = parent.frame(), n.draws = 5000,
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
          burnin = 0, thin = 0, verbose = TRUE)
}

\arguments{
  \item{formulae}{
    A list of formulae where the first formula specifies the
    (pre-treatment) covariates in the outcome model (the latent
    compliance covariate will be added automatically), the second
    formula specifies the compliance model, and the third formula
    defines the covariate specification for the model for missing-data
    mechanism (the latent compliance covariate will be added
    automatically). For the outcome model, the formula should take the
    two-sided standard R \command{formula} where the outcome variable is
    specified in the left hand side of the formula which is then
    separated by \code{~} from the covariate equation in the right hand
    side, e.g., \code{y ~ x1 + x2}. For the compliance and missing-data
    mechanism models, the one-sided \command{formula} should be used
    where the left hand side is left unspecified, e.g., \code{~ x1 + x2}.
  }
  \item{Z}{
    A randomized encouragement variable, which should be a binary
    variable in the specified data frame.
  }
  \item{D}{
    A treatment variable, which should be a binary variable in the
    specified data frame.
  }
  \item{data}{
    A data frame which contains the variables that appear in
    the model formulae (\code{formulae}), the encouragement variable
    (\code{Z}), and the treatment variable (\code{D}).
  }
  \item{n.draws}{
    The number of MCMC draws. The default is \code{5000}.
  }
  \item{param}{
    A logical variable indicating whether the Monte Carlo draws of the
    model parameters should be saved in the output object. The default
    is \code{TRUE}.
  }
  \item{in.sample}{
    A logical variable indicating whether or not the sample average
    causal effect should be calculated using the observed
    potential outcome for each unit. If it is set to \code{FALSE},
    then the population average causal effect will be calculated. The
    default is \code{FALSE}.
  }
  \item{model.c}{
    The model for compliance. Either \code{logit} or \code{probit} model
    is allowed. The default is \code{probit}.
  }
  \item{model.o}{
    The model for outcome. The following five models are allowed:
    \code{logit}, \code{probit}, \code{oprobit} (ordered probit regression),
    \code{gaussian} (gaussian regression), \code{negbin} (negative
    binomial regression), and \code{twopart} (two part model where the
    first part is the probit regression for \eqn{Pr(Y>0|X)} and the
    second part models \eqn{p(log(Y)|X, Y>0)} using the gaussian
    regression). The default is \code{probit}.
  }
  \item{model.r}{
    The model for (non)response. Either \code{logit} or \code{probit}
    model is allowed. The default is \code{probit}.
  }
  \item{tune.c}{
    Tuning constants for fitting the compliance model. These
    positive constants are used to tune the (random-walk)
    Metropolis-Hastings algorithm to fit the logit model. Use either a
    scalar or a vector of constants whose length equals that of the
    coefficient vector. The default is \code{0.01}.
  }
  \item{tune.o}{
    Tuning constants for fitting the outcome model. These
    positive constants are used to tune the (random-walk)
    Metropolis-Hastings algorithm to fit logit, ordered
    probit, and negative binomial models. Use either a
    scalar or a vector of constants whose length equals that of the
    coefficient vector for logit and negative binomial models. For the
    ordered probit model, use either a scalar or a vector of constants
    whose length equals that of cut-point parameters to be
    estimated. The default is \code{0.01}.  
  }
  \item{tune.r}{
    Tuning constants for fitting the (non)response model. These
    positive constants are used to tune the (random-walk)
    Metropolis-Hastings algorithm to fit the logit model. Use either a
    scalar or a vector of constants whose length equals that of the
    coefficient vector. The default is \code{0.01}.
  }
  \item{tune.v}{
    A scalar tuning constant for fitting the variance component of the
    negative binomial (outcome) model. The default is \code{0.01}.
    
  }
  \item{p.mean.c}{
    Prior mean for the compliance model. It should be either a scalar or
    a vector of appropriate length. The default is \code{0}.
  }
  \item{p.prec.c}{
    Prior precision for the compliance model. It should be either a
    positive scalar or a positive semi-definite matrix of appropriate
    size. The default is \code{0.001}.
  }
  \item{p.mean.o}{
    Prior mean for the outcome model. It should be either a scalar or
    a vector of appropriate length. The default is \code{0}.
  }
  \item{p.prec.o}{
    Prior precision for the outcome model. It should be either a
    positive scalar or a positive semi-definite matrix of appropriate
    size. The default is \code{0.001}.
  }
  \item{p.mean.r}{
    Prior mean for the (non)response model. It should be either a scalar or
    a vector of appropriate length. The default is \code{0}.
  }
  \item{p.prec.r}{
    Prior precision for the (non)response model. It should be either a
    positive scalar or a positive semi-definite matrix of appropriate
    size. The default is \code{0.001}.
  }
  \item{p.df.o}{
    A positive integer. Prior degrees of freedom parameter for the inverse
    chisquare distribution in the gaussian and twopart (outcome) models. The
    default is \code{10}.
  }
  \item{p.scale.o}{
    A positive scalar. Prior scale parameter for the inverse chisquare
    distribution (for the variance) in the gaussian and twopart
    (outcome) models. For the negative binomial (outcome) model, this is
    used for the scale parameter of the inverse gamma distribution. The
    default is \code{1}.
  }
  \item{p.shape.o}{
    A positive scalar. Prior shape for the inverse chisquare
    distribution in the negative binomial (outcome) model. The default
    is \code{1}.
  }
  \item{mda.probit}{
    A logical variable indicating whether to use marginal data
    augmentation for probit models. The default is \code{TRUE}.
  }
  \item{coef.start.c}{
    Starting values for coefficients of the compliance model.     
    It should be either a scalar or a vector of appropriate length. The
    default is \code{0}. 
  }
  \item{coef.start.o}{
    Starting values for coefficients of the outcome model.     
    It should be either a scalar or a vector of appropriate length. The
    default is \code{0}. 
  }
  \item{coef.start.r}{
    Starting values for coefficients of the (non)response model.     
    It should be either a scalar or a vector of appropriate length. The
    default is \code{0}. 
  }
  \item{tau.start.o}{
    Starting values for thresholds of the ordered probit (outcome) model.
    It should be either a scalar or a vector of appropriate length.
    The default is \code{0}. 
  }
  \item{var.start.o}{
    starting values for variances of the outcome model.
  }
  \item{burnin}{
    number of burnins.
  }
  \item{thin}{
    number of thinning.
  }
  \item{verbose}{
    some print outs.
  }
  
}

\examples{
  n <- 5
  Y <- rep(c(0,1), n)
  X <- rnorm(n*2)
  #res <- Noncomp.bayes(list(Y ~ X, ~ X, ~ X), D = 1-Y, Z = Y)
}

\references{
  Frangakis, C. E. and Rubin, D. B. (1999). \dQuote{Addressing
    Complications of Intention-to-Treat Analysis in the Combined
    Presence of All-or-None Treatment Noncompliance and Subsequent
    Missing Outcomes.} \emph{Biometrika}, Vol. 86, No. 2, pp. 365-379.
  
  Horiuchi, Yusaku, Kosuke Imai, and Naoko Taniguchi (2007). \dQuote{Designing
    and Analyzing Randomized Experiments: Application to a Japanese Election
    Survey Experiment.} \emph{American Journal of Political Science}, Vol. 51,
  No. 3 (July), pp. 669-687.
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
}

\keyword{models}





