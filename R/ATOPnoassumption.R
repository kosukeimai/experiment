###
### Calculates the no-assumption bounds for the  ATOP in the presence of missing
### response  under matched-pairs design
###



#' Bounding the  ATOP when some of the Outcome Data are
#' Missing Under the Matched-Pairs Design
#' 
#' This function computes the no assumption bounds on the average treatment effect among always-observed pairs (ATOP) when
#' some of the outcome data are missing. The confidence intervals for the
#' ATOP are also computed.
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param Ya  A vector of the outcomes of the first unit in the matched pairs. The missing values for \code{Ya} should be coded
#' as \code{NA}.
#' @param Yb  A vector of the outcomes of the second unit in the matched pairs. The missing values for \code{Yb} should be coded
#' as \code{NA}.
#' @param Ra A vector of the missing data indicators of the first unit in the matched pairs. 
#' @param Rb A vector of the missing data indicators of the second unit in the matched pairs.
#' @param Ta A vector of the treatment conditions of the first unit in the matched pairs. 
#' @param Tb A vector of the treatment conditions of the second unit in the matched pairs.
#' @param l The lower limit of the outcome.
#' @param u The upper limit of the outcome.
#' @param alpha A positive scalar that is less than or equal to 0.5. This will
#' determine the (1-\code{alpha}) level of confidence intervals. The default is
#' \code{0.05}.
#' @param rep The number of repetitions for bootstraping.
#' @return A list of class \code{ATOPnoassumption} which contains the following items:
#' \item{LB}{ The lower bound for the ATOP.  } \item{UB}{ The upper bound for the ATOP.   }
#' \item{LB.CI}{ The lower limit of the confidence interval for the ATOP.   }
#'\item{UB.CI}{ The upper limit of the confidence interval for the ATOP.    } 
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu}; 
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai and Zhichao Jiang (2018).
#' \dQuote{A Sensitivity Analysis for Missing Outcomes Due to 
#' Truncation-by-Death under the Matched-Pairs Design}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords matched-pairs design
#' @examples 
#' data(seguro)
#' attach(seguro)
#' ATOPnoassumption(Ya,Yb,Ra,Rb,Ta,Tb,l=0,u=1,alpha=0.05,rep=100)
#' @export ATOPnoassumption

ATOPnoassumption = function(Ya,Yb,Ra,Rb,Ta,Tb,l,u,alpha,rep){
  if(!(is.vector(Ya)&is.vector(Yb)&is.vector(Ra)&is.vector(Rb)&is.vector(Ta)&is.vector(Tb)))
    stop('Data should be input as vectors')
  if (length(unique(apply(cbind(Ya,Yb,Ra,Rb,Ta,Tb),2,length)))!=1)
    stop('Vectors should have the same length')
  
  omega1 <- mean(c(Ya[Ta==1&Ra==1&Rb==1],Yb[Tb==1&Rb==1&Ra==1]))   	                      
  omega0 <- mean(c(Ya[Ta==0&Ra==1&Rb==1],Yb[Tb==0&Rb==1&Ra==1])) 
  
  pi <- mean(Ra==1&Rb==1)
  N <- length(Ya)
  
  #####  
  
  indL1 <- 1
  indL2 <- 1
  indU1 <- 1
  indU2 <- 1
  if (pi<1/2 | f1(omega1,2-1/pi,u)<l){indL1=2}
  if (pi<1/2 | f2(omega0,2-1/pi,l)>u){indL2=2}
  if (pi<1/2| f2(omega1,2-1/pi,l)>u){indU1=2}
  if (pi<1/2 | f1(omega0,2-1/pi,u)<l){indU2=2}
  
  
  LB <- ifelse(indL1==1,f1(omega1,2-1/pi,u),l)-ifelse(indL2==1,f2(omega0,2-1/pi,l),u)
  UB <- ifelse(indU1==1,f2(omega1,2-1/pi,l),u)-ifelse(indU2==1,f1(omega0,2-1/pi,u),l)
  
  ####  bootstrap
  LB.boots <- rep(0,rep)
  UB.boots <- rep(0,rep)
  for (i in 1 :rep){
    index <- sample(1:N,N,replace=TRUE)
    Ya.ind <- Ya[index]
    Yb.ind <- Yb[index]
    Ra.ind <- Ra[index]
    Rb.ind <- Rb[index]
    Ta.ind <- Ta[index]
    Tb.ind <- Tb[index]
    omega1.boots <- mean(c(Ya.ind[Ta.ind==1&Ra.ind==1&Rb.ind==1],Yb.ind[Tb.ind==1&Rb.ind==1&Ra.ind==1]))
    omega0.boots <- mean(c(Ya.ind[Ta.ind==0&Ra.ind==1&Rb.ind==1],Yb.ind[Tb.ind==0&Rb.ind==1&Ra.ind==1]))
    pi.boots <- mean(Ra.ind==1&Rb.ind==1)
    LB.boots[i] <- ifelse(indL1==1,f1(omega1.boots,2-1/pi.boots,u),l)-ifelse(indL2==1,f2(omega0.boots,2-1/pi.boots,l),u)
    UB.boots[i] <- ifelse(indU1==1,f2(omega1.boots, 2-1/pi.boots,l),u)-ifelse(indU2==1,f1(omega0.boots,2-1/pi.boots,u),l)
  }
  #### Calculate CI
  sd.LB <- sd(LB.boots)
  sd.UB <- sd(UB.boots)
  C <- uniroot(CalC,hat.LB=LB,hat.UB=UB,sigma.LB=sd.LB,sigma.UB=sd.UB,alpha=alpha,lower=0,upper=100)$root
  LB.CI <- max(LB-C*sd.LB,l-u)
  UB.CI <- min(UB+C*sd.UB,u-l)
  return(list(LB=LB,UB=UB,LB.CI=LB.CI,UB.CI=UB.CI))
}