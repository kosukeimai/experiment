# experiment: R package for designing and analyzing randomized experiments [![Build Status](https://travis-ci.org/kosukeimai/experiment.svg?branch=master)](https://travis-ci.org/kosukeimai/experiment) [![CRAN Version](http://www.r-pkg.org/badges/version/experiment)](https://CRAN.R-project.org/package=experiment)
R package experiment Provides various statistical methods for designing and analyzing randomized experiments. One functionality of the package is the implementation of randomized-block and matched-pair designs based on possibly multivariate pre-treatment covariates. The package also provides the tools to analyze various randomized experiments including cluster randomized experiments, two-stage randomized experiments, randomized experiments with noncompliance, and randomized experiments with missing data.

## Examples

Complete Randomization is supported via the `completeRandomization` function, which takes `data` as an argument and assigns `n` units to treatment:

```
## assigns half of the units to treatment, half to control
completeRandomization(data) 

## assigns 10 units to treatment, remainder to control
completeRandomization(data, n = 10) 

## 10 units to treatment, 10 to control
completeRandomization(data, n = c(10, 10))

## assigns 10 units to first group, 5 to other two groups.
completeRandomization(data, n = c(10, 5, 5))
```