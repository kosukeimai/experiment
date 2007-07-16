\name{ATEbounds}

\alias{ATEbounds}

\title{Bounding the Average Treatment Effect when some of the Outcome
  Data are Missing
}

\description{
  This function computes the sharp bounds on the average treatment
  effect when some of the outcome data are missing. The confidence
  intervals for the bounds are also computed. 
}

\usage{
  ATEbounds(formula, data = parent.frame(), maxY = NULL,
            minY = NULL, alpha = 0.05, n.reps = 0,
            strata = NULL, ratio = NULL, survey = NULL, ...)
}

\arguments{
  \item{formula}{A formula of the form \code{Y ~ X} where \code{Y} is
    the name of the outcome variable and \code{X} is the
    name of the (randomized) treatment variable. \code{X} should be a
    factor variable but its value can take more than two levels. The missing
    values for \code{Y} should be coded as \code{NA}.
  }
  \item{data}{
    A data frame containing the relevant variables. 
  }
  \item{maxY}{
    A scalar. The maximum value of the outcome variable. The default is the
    maximum sample value.
  }
  \item{minY}{
    A scalar. The minimum value of the outcome variable. The default is the
    minimum sample value.
  }
  \item{alpha}{
    A positive scalar that is less than or equal to 0.5. This will
    determine the (1-\code{alpha}) level of confidence intervals. The
    default is \code{0.05}.
  }
  \item{strata}{
    The variable name indicating strata. If this is specified, the
    quantities of interest will be first calculated within each strata
    and then aggregated. The default is \code{NULL}.
  }
  \item{ratio}{
    A \eqn{J \times M} matrix of probabilities where \eqn{J} is the
    number of strata and \eqn{M} is the number of treatment and control
    groups. Each element of the matrix specifies the probability of a
    unit falling into that category. The default is \code{NULL} in which
    case the sample estimates of these probabilities are used for
    computation. 
  }
  \item{survey}{
    The variable name for survey weights. The default is \code{NULL}.
  }
  \item{n.reps}{
    A positive integer. The number of bootstrap replicates used for
    the construction of confidence intervals via B-method of Berran
    (1988). If it equals zero, the confidence intervals will not be
    constructed. 
  }
  \item{...}{
    The arguments passed to other functions.
  }
}

\details{
  For the details of the method implemented by this function, see the
  references. 
}

\value{
  A list of class \code{ATEbounds} which contains the following items:
  \item{call}{
    The matched call.
  }
  \item{Y}{
    The outcome variable.
  }
  \item{D}{
    The treatment variable.
  }
  \item{bounds}{
    The point estimates of the sharp bounds on the average treatment
    effect.
  }
  \item{bounds.Y}{
    The point estimates of the sharp bounds on the outcome variable
    within each treatment/control group.
  }
  \item{bmethod.ci}{
    The B-method confidence interval of the bounds on the average
    treatment effect.
  }
  \item{bonf.ci}{
    The Bonferroni confidence interval of the bounds on the average
    treatment effect.
  }
  \item{bonf.ci.Y}{
    The Bonferroni confidence interval of the bounds on the outcome
    variable within each treatment/control group.
  }
  \item{bmethod.ci.Y}{
    The B-method confidence interval of the bounds on the outcome
    variable within each treatment/control group.
  }
  \item{maxY}{
    The maximum value of the outcome variable used in the computation.
  }
  \item{minY}{
    The minimum value of the outcome variable used in the computation.
  }
  \item{nobs}{
    The number of observations.
  }
  \item{nobs.Y}{
    The number of observations within each treatment/control group.
  }
  \item{ratio}{
    The probability of treatment assignment (within each strata if
    \code{strata} is specified) used in the computation.
  }
    
}

\references{
  Horowitz, Joel L. and Charles F. Manski. (1998).
  \dQuote{Censoring of Outcomes and Regressors due to Survey
    Nonresponse: Identification and Estimation Using Weights and
    Imputations.} \emph{Journal of Econometrics}, Vol. 84, pp.37-58.
  
  Horowitz, Joel L. and Charles F. Manski. (2000).
  \dQuote{Nonparametric Analysis of Randomized Experiments
    With Missing Covariate and Outcome Data.} \emph{Journal of the
    Americal Statistical Association}, Vol. 95, No. 449, pp.77-84.

  Harris-Lacewell, Melissa, Kosuke Imai, and Teppei Yamamoto. (2007).
    \dQuote{Racial Gaps in the Responses to Hurricane Katrina:
      An Experimental Study},
   \emph{Technical Report}. Department of Politics, Princeton University.
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
}

\keyword{design}





