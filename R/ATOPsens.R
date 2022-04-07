###
### Sensitivity analysis for the  ATOP in the presence of missing
### response  under matched-pairs design
###



#' Sensitivity analysis for the ATOP when some of the Outcome Data are
#' Missing Under the Matched-Pairs Design
#' 
#' This function computes the bounds on the average treatment effect among always-observed pairs (ATOP)
#'  with pre-specified sensivity parameters when
#' some of the outcome data are missing. The sensivity parameter characterizes the degree of the within-pair similarity.
#' The confidence intervals for the
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
#' @param gamma The sensitivity parameter which charaterizes the degree of the within-pair similarity.
#' @param l The lower limit of the outcome.
#' @param u The upper limit of the outcome.
#' @param alpha A positive scalar that is less than or equal to 0.5. This will
#' determine the (1-\code{alpha}) level of confidence intervals. The default is
#' \code{0.05}.
#' @param rep The number of repetitions for bootstraping.
#' @return A list of class \code{ATOPsens} which contains the following items:
#' \item{LB}{ The lower bound for the ATOP.  } \item{UB}{ The upper bound for the ATOP.   }
#' \item{LB.CI}{ The lower limit of the confidence interval for the ATOP.   }
#'\item{UB.CI}{ The upper limit of the confidence interval for the ATOP.    } 
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai and Zhichao Jiang (2018).
#' \dQuote{A Sensitivity Analysis for Missing Outcomes Due to 
#' Truncation-by-Death under the Matched-Pairs Design}, \emph{Statistics in Medicine}.
#' @keywords matched-pairs design
#' @examples 
#' data(seguro)
#' attach(seguro)
#' ATOPobs(Ya,Yb,Ra,Rb,Ta,Tb,gamma=0.95,kappa1=1,kappa0=1,l=0,u=1,alpha=0.05,rep=100)
#' @export ATOPsens


ATOPsens <- function(Ya,Yb,Ra,Rb,Ta,Tb,gamma,l,u,alpha,rep){
  if(!(is.vector(Ya)&is.vector(Yb)&is.vector(Ra)&is.vector(Rb)&is.vector(Ta)&is.vector(Tb)))
    stop('Data should be input as vectors')
  if (length(unique(apply(cbind(Ya,Yb,Ra,Rb,Ta,Tb),2,length)))!=1)
    stop('Vectors should have the same length')
  
  omega1 <- mean(c(Ya[Ta==1&Ra==1&Rb==1],Yb[Tb==1&Rb==1&Ra==1]))   	                      
  omega0 <- mean(c(Ya[Ta==0&Ra==1&Rb==1],Yb[Tb==0&Rb==1&Ra==1])) 
  
  alpha1 <- (sum(Ta==1&Ra==1)+sum(Tb==1&Rb==1))/(sum(Ta==1)+sum(Tb==1))
  alpha0 <- (sum(Ta==0&Ra==1)+sum(Tb==0&Rb==1))/(sum(Ta==0)+sum(Tb==0))
  
  pi <- mean(Ra==1&Rb==1)
  N<-length(Ya)

  Delta1 <- max(2*pi-1+gamma*(1-alpha1),2*pi-1+gamma*(1-alpha0),pi-(1-gamma)*(alpha1+alpha0), 2*pi-(2-gamma)*alpha1,  2*pi-(2-gamma)*alpha0,pi-(1-gamma)*(2-alpha1-alpha0),pi-(1-gamma)*(1-abs(alpha1-alpha0)))/pi
  
  ind <- which.max(c(2*pi-1+gamma*(1-alpha1),2*pi-1+gamma*(1-alpha0),pi-(1-gamma)*(alpha1+alpha0),2*pi-(2-gamma)*alpha1,  2*pi-(2-gamma)*alpha0,pi-(1-gamma)*(2-alpha1-alpha0),pi-(1-gamma)*(1-abs(alpha1-alpha0)) ))
  
  
  #####  
  
  indL1 <- ind
  indL2 <- ind
  indU1 <- ind
  indU2 <- ind
  if (Delta1<0 | f1(omega1,Delta1,u)<l){indL1=0}
  if (Delta1<0 | f2(omega0,Delta1,l)>u){indL2=0}
  if (Delta1<0| f2(omega1,Delta1,l)>u){indU1=0}
  if (Delta1<0 | f1(omega0,Delta1,u)<l){indU2=0}
  
  LB <- max(ifelse(indL1!=0,f1(omega1,Delta1,u),l)-ifelse(indL2!=0,f2(omega0,Delta1,l),u),l-u)
  UB <- min(ifelse(indU1!=0,f2(omega1,Delta1,l),u)-ifelse(indU2!=0,f1(omega0,Delta1,u),l),u-l)
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
    alpha1.boots <- (sum(Ta.ind==1&Ra.ind==1)+sum(Tb.ind==1&Rb.ind==1))/(sum(Ta.ind==1)+sum(Tb.ind==1))
    alpha0.boots <- (sum(Ta.ind==0&Ra.ind==1)+sum(Tb.ind==0&Rb.ind==1))/(sum(Ta.ind==0)+sum(Tb.ind==0))
    pi.boots <- mean(Ra.ind==1&Rb.ind==1)
    Delta1.boots <- c(2*pi.boots-1+gamma*(1-alpha1.boots),2*pi.boots-1+gamma*(1-alpha0.boots),pi.boots-(1-gamma)*(alpha1.boots+alpha0.boots),2*pi.boots-(2-gamma)*alpha1.boots,  2*pi.boots-(2-gamma)*alpha0.boots,pi.boots-(1-gamma)*(2-alpha1.boots-alpha0.boots),pi.boots-(1-gamma)*(1-abs(alpha1.boots-alpha0.boots)) )[ind]/pi.boots
    
    
    
    LB.boots[i] <- max(ifelse(indL1!=0&Delta1.boots>0,f1(omega1.boots,Delta1.boots,u),l)-ifelse(indL2!=0&Delta1.boots>0,f2(omega0.boots,Delta1.boots,l),u),l-u)
    UB.boots[i] <- min(ifelse(indU1!=0&Delta1.boots>0,f2(omega1.boots, Delta1.boots,l),u)-ifelse(indU2!=0&Delta1.boots>0,f1(omega0.boots,Delta1.boots,u),l),u-l)
  }
  #### Calculate CI
  sd.LB <- sd(LB.boots)
  sd.UB <- sd(UB.boots)
  C <- uniroot(CalC,hat.LB=LB,hat.UB=UB,sigma.LB=sd.LB,sigma.UB=sd.UB,alpha=alpha,lower=0,upper=100)$root
  LB.CI <- max(LB-C*sd.LB,l-u)
  UB.CI <- min(UB+C*sd.UB,u-l)
  return(list(LB=LB,UB=UB,LB.CI=LB.CI,UB.CI=UB.CI))
}
