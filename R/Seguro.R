#' Data from the Mexican universal health insurance program, Seguro Popular.
#'
#' This data set contains the outcome, missing indicator and the treatment for the application 
#' in  Imai and Jiang (2018)
#'
#' @format A data frame with 14,902 rows and 6 variables:
#' \describe{
#'   \item{Ya}{Satisfaction for the first unit in the matched pairs}
#'   \item{Yb}{Satisfaction for the second unit in the matched pairs}
#'   \item{Ra}{Missing indicator for the first unit in the matched pairs}
#'   \item{Rb}{Missing indicator for the second unit in the matched pairs}
#'   \item{Ta}{Treatment assignment for the first unit in the matched pairs}
#'   \item{Tb}{Treatment assignment for the second unit in the matched pairs}
#'   #' }
#'
#' @docType data
#' @keywords datasets
#' @name seguro
#' @examples
#' data(seguro)
#' attach(seguro)
#' ATOPnoassumption(Ya,Yb,Ra,Rb,Ta,Tb,l=0,u=1,alpha=0.05,rep=1000)
#' ATOPsens(Ya,Yb,Ra,Rb,Ta,Tb,gamma=0.95,l=0,u=1,alpha=0.05,rep=1000)
#' ATOPobs(Ya,Yb,Ra,Rb,Ta,Tb,gamma=0.95,kappa1=1,kappa0=1,l=0,u=1,alpha=0.05,rep=1000)
"seguro"