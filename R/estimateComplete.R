#' Estimator for Completely Randomized Experiments
#' 
#' This function conducts an analysis of a completely randomized experiment
#' and returns the effect estimate and standard errors.
#' 
#' Full provided list of \code{std.error} as provided in \code{sandwich} (Zeileis 2004) is: 
#' const (the default R \code{lm} implementation), HC, HC1, HC2, HC3 (MacKinnon and White 1985), 
#' HC4 (Cribari-Neto 2004), HC4m (Cribari-Neto, Souza & Vasconcellos 2007), 
#' HC5 (HC5: Cribari-Neto & Da Silva 2011). Further details are provided in the 
#' provided references and \code{sandwich} package documentation.
#'  
#' @aliases estimateComplete estimateCompleteRandomization
#' @param data A data frame containing treatment and covariates.
#' @param formula An object of class formula.
#' @param std.error Defaults to HC3 (MacKinnon and White 1985), as the homoskedasticity 
#' assumption is not often defensible. User can supply options from \code{sandwich} package.
#' @return A list of class \code{estimateComplete} which contains the following items:
#' \item{call}{ the matched call. } \item{estimates}{The estimated coefficients and
#' standard errors.}\item{data}{The data frame that was used in the analysis.}
#' @author Kosuke Imai, Professor of Government and Statistics, Harvard University
#' \email{imai@@harvard.edu};
#' @author Tyler Simko, PhD Student, Harvard University
#' \email{tsimko@@g.harvard.edu};
#' @importFrom sandwich vcovHC
#' @references 
#' Achim Zeileis (2004). Econometric Computing with HC and HAC 
#' Covariance Matrix Estimators. Journal of Statistical Software 11(10), 1-17. 
#' URL \url{http://www.jstatsoft.org/v11/i10/.}
#' 
#' Cribari-Neto F. (2004), Asymptotic Inference under Heteroskedasticity of Unknown Form. 
#' Computational Statistics & Data Analysis 45, 215-233. 
#' \url{http://doi.org/10.1016/S0167-9473(02)00366-3}.
#' 
#' Cribari-Neto F., Da Silva W.B. (2011), A New Heteroskedasticity-Consistent Covariance Matrix 
#' Estimator for the Linear Regression Model. Advances in Statistical Analysis, 95(2), 129-146.
#' \url{https://doi.org/10.1007/s10182-010-0141-2}.
#' 
#' Cribari-Neto F., Souza T.C., Vasconcellos, K.L.P. (2007), Inference under Heteroskedasticity 
#' and Leveraged Data. Communications in Statistics - Theory and Methods, 36, 1877-1888. Errata: 
#' 37, 3329-3330, 2008. \url{https://doi.org/10.1080/03610920601126589}.
#' 
#' MacKinnon, James, and Halbert White. (1985), "Some Heteroskedasticity-Consistent 
#' Covariance Matrix Estimators with Improved Finite Sample Properties." Journal 
#' of Econometrics 29 (3): 305-25. \url{https://doi.org/10.1016/0304-4076(85)90158-7}.
#' 
#' White H. (1980), A Heteroskedasticity-Consistent Covariance Matrix and a Direct Test 
#' for Heteroskedasticity. Econometrica 48, 817-838. \url{http://doi.org/10.2307/1912934}.
#' 
#' @examples
#' ## PlantGrowth data provided with R, HC3 standard errors
#' estimateComplete(weight ~ group, PlantGrowth)
#' 
#' ## specify standard errors
#' estimateComplete(weight ~ group, PlantGrowth, std.error = "HC2")
#' 
#' ## specify homoskedasticity
#' estimateComplete(weight ~ group, PlantGrowth, std.error = "const")
#' @export estimateComplete

estimateComplete <- function(formula, data, std.error = "HC3") {
  call <- match.call()
  
  ## check if valid option for standard errors
  se.options <- c("const", "HC", "HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")
  if (!(std.error %in% se.options)) 
    stop("Invalid option for standard errors. Run ?estimateComplete for list of options.")
  
  ## fit model
  lm.fit <- lm(formula, data = data)
  sum <- summary(lm.fit)
  
  ## calculate appropriate standard errors
  se <- vcovHC(lm.fit, type = std.error)
  
  ## degrees of freedom
  df <- df.residual(lm.fit)
  
  ## t-values
  t.val <- coef(lm.fit) / sqrt(diag(se))
  
  ## p-values
  p.val <- 2 * pt(-abs(t.val), df)
  
  ## collect results
  est <- data.frame(Estimates = lm.fit$coefficients, 
                    SE = sqrt(diag(se)),
                    t.val = t.val,
                    p.val = p.val)
  
  ## put standard error name in data frame for clarity
  se.name <- paste("SE (", std.error, ")", sep = "")
  attr(est, "names")[2] <- se.name
  
  res <- list(call = call, estimates = est, 
              data = data, df = df,
              r.squared = sum$r.squared,
              adj.r.squared = sum$adj.r.squared)
  
  class(res) <- "estimateComplete"
  
  return(res)
}