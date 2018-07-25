# experiment: R package for designing and analyzing randomized experiments [![Build Status](https://travis-ci.org/kosukeimai/experiment.svg?branch=master)](https://travis-ci.org/kosukeimai/experiment) [![CRAN Version](http://www.r-pkg.org/badges/version/experiment)](https://CRAN.R-project.org/package=experiment)
R package experiment Provides various statistical methods for designing and analyzing randomized experiments. One functionality of the package is the implementation of randomized-block and matched-pair designs based on possibly multivariate pre-treatment covariates. The package also provides the tools to analyze various randomized experiments including cluster randomized experiments, two-stage randomized experiments, randomized experiments with noncompliance, and randomized experiments with missing data.

## Randomization

Complete Randomization is supported via the `completeRandomization` function, which takes `data` as an argument and assigns `n` units to treatment. Multiple treatment arms are supported by providing a vector for `n`. For a sample dataset, `data`:

```
## example data
data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

## assigns half of the units to treatment, half to control
completeRandomization(data) 

## assigns 10 units to treatment, remainder to control
completeRandomization(data, n = 10) 

## 10 units to treatment, 10 to control 
## (equivalent to above line if nrow(data) is 20)
completeRandomization(data, n = c(10, 10))

## assigns 10 units to first group, 5 to other two groups.
completeRandomization(data, n = c(10, 5, 5))
```

## Analysis

The analysis of complete randomized experiments is provided via the `estimateComplete` function, which takes a formula and provides several different options for standard errors:

```{r}
## PlantGrowth data provided with R
data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

## HC3 standard errors by default
estimateComplete(weight ~ group, PlantGrowth, std.error = "HC2")

## specify homoskedasticity
estimateComplete(weight ~ group, PlantGrowth, std.error = "const")
```