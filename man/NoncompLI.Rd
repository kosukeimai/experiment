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
    model parameters should be saved in the output object.
  }
  \item{in.sample}{
    in sample qoi?
  }
  \item{model.c}{
    model for compliance
  }
  \item{model.o}{
    model for outcome
  }
  \item{model.r}{
    model for response
  }
  \item{tune.c}{
    tuning constants for compliance model.
  }
  \item{tune.o}{
    tuning constants for outcome model.
  }
  \item{tune.r}{
    tuning constants for response model.
  }
  \item{tune.v}{
    tuning constants for variance.
  }
  \item{p.mean.c}{
    prior mean for the compliance model.
  }
  \item{p.prec.c}{
    prior precision for the compliance model.
  }
  \item{p.mean.o}{
    prior mean for the outcome model.
  }
  \item{p.prec.o}{
    prior precision for the outcome model.
  }
  \item{p.mean.r}{
    prior mean for the response model.
  }
  \item{p.prec.r}{
    prior precision for the response model.
  }
  \item{p.df.o}{
    prior degrees of freedom for the outcome model.
  }
  \item{p.scale.o}{
    prior scale for the outcome model.
  }
  \item{p.shape.o}{
    prior shape for the outcome model.
  }
  \item{mda.probit}{
    Use marginal data augmentation for probit models?
  }
  \item{coef.start.c}{
    starting values for coefficients of the compliance model.
  }
  \item{coef.start.o}{
    starting values for coefficients of the outcome model.
  }
  \item{coef.start.r}{
    starting values for coefficients of the response model.
  }
  \item{tau.start.o}{
    starting values for thresholds of the outcome model.
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

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
}

\keyword{models}





