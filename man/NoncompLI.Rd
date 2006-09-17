\name{Noncomp.bayes}

\alias{Noncomp.bayes}

\title{Analysis of Randomized Experiments with Noncompliance and
  Nonignorable Missing Outcomes}

\description{
  This function estimates the average causal effects for randomized
  experiments with noncompliance and nonignorable missing outcomes.
}

\usage{
  Noncomp.bayes(formulae, Z, D, data = parent.frame(),
                          n.draws = 5000, param = TRUE,
                          in.sample = FALSE, model.c = "probit", 
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
    a list of formulae.
  }
  \item{Z}{
    randomized encouragement.
  }
  \item{D}{
    Binary treatment variable.
  }
  \item{data}{
    a data frame.
  }
  \item{n.draws}{
    number of MCMC draws.
  }
  \item{param}{
    store parameters?
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

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
}

\keyword{models}





