\name{ATEbounds}

\alias{ATEbounds}

\title{Bounding the Average Treatment Effect when the Outcome Data are
       Missing}

\description{
  This function computes the sharp bounds on the average treatment
  effect when the outcome data are missing.
}

\usage{
  ATEbounds(formula, data = parent.frame(), maxY = NULL,
            minY = NULL, alpha = 0.05, n.reps = 0,
            strata = NULL, ratio = NULL, survey = NULL, ...)
}

\arguments{
  \item{data}{
    A data frame containing the observations to which the treatments
    are randomly assigned.
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
    determine the (1-alpha) level of confidence intervals.
  }
  \item{n.reps}{
    A positive integer. The number of bootstrap replicates used for
    the construction of confidence intervals via B-method of Berran
    (1988). If it equals zero, the confidence intervals will not be
    constructed. 
  }
}

\details{
  Randomized-block designs refer to the complete randomization of the
  treatment within the pre-specified blocks which contain multiple
  observations. Matched-pair designs refer to the randomization of the
  binary treatment variable within the pre-specified pair of
  observations.
}

\value{
  A list of class \code{randomize} which contains the following items:
  \item{call}{
    the matched call.
  }
  \item{treatment}{
    The vector of randomized treatments.
  }
  \item{data}{
    The data frame that was used to conduct the randomization.
  }
  \item{block}{
    The blocking variable that was used to implement randomized-block designs.
  }
  \item{match}{
    The matching variable that was used to implement matched-pair
    designs.
  }
  \item{block.id}{
    The variable indicating which observations belong to which blocks in
    randomized-block designs.
  }
  \item{match.id}{
    The variable indicating which observations belong to which pairs in
    matched-pair designs.
  }
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://imai.princeton.edu};
}

\keyword{design}





