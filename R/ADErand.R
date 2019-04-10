###
### Randomization-based method for the intention-to-treat (ITT)  average direct and spillover effects
### 
###



#' Randomization-based method for the intention-to-treat (ITT)  average direct and spillover effects
#'
#' 
#' This function computes the point estimates and variance estimates of the  average direct effect (ADE)  and the average spillover effect (ASE).
#' The estimators calculated using this function are either individual weighted or cluster-weighted. The point estimates and variances of ITT effects are also included. 
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @param individual  A binary variable with TRUE for  individual-weighted estimators and FALSE for cluster-weighted estimators.
#' @return A list of class \code{ADErand} which contains the following items:
#'\item{ADE1}{ The point estimate of ADE(1).  } \item{ADE0}{ The point estimate of ADE(0).  }  
#'\item{var.ADE1}{ The  variance estimate of ADE(1).   } 
#'\item{var.ADE0}{ The  variance estimate of ADE(0).   } 
#'\item{ASE1}{ The point estimate of ASE(1).  } \item{ASE0}{ The point estimate of ASE(0).  } 
#'\item{var.ASE1}{ The  variance estimate of ASE(1).   } 
#'\item{var.ASE0}{ The  variance estimate of ASE(0).   } 
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@Princeton.Edu}, \url{https://imai.princeton.edu};
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai, Zhichao Jiang and Anup Malani (2018).
#' \dQuote{Causal Inference with Interference and Noncompliance in the Two-Stage Randomized Experiments}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords two-stage randomized experiments
#' @export ADErand


ADErand=function(data,individual=1){
	## transform the data into list 	
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id=unique(data$id)	
  n.cluster=length(cluster.id)	
  Z=vector("list", n.cluster) 	
  Y=vector("list", n.cluster) 
  A=rep(0,n.cluster)
  for (i in 1:n.cluster){
  	Z[[i]]=as.numeric(data$Z[data$id==cluster.id[i]])
  	Y[[i]]=data$Y[data$id==cluster.id[i]]
  	if (length(unique(data$A[data$id==cluster.id[i]]))!=1){
  		stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
  	}
  	A[i]=data$A[data$id==cluster.id[i]][1]
  }

	
	
  n=sapply(Z,length)
  N=sum(n)
  J=length(n)
  if (individual==1){
    for ( i in 1:n.cluster){
      Y[[i]]=Y[[i]]*n[i]*J/N
    }
  }
  
  
  est.Yj00=rep(0,J)
  est.Yj00= sapply(Difflist(Y,Productlist(Y,Z)),sum)/(n-sapply(Z,sum))*(1-A)
  est.Yj01=rep(0,J)
  est.Yj01= sapply(Difflist(Y,Productlist(Y,Z)),sum)/(n-sapply(Z,sum))*(A)
  est.Yj10=rep(0,J)
  est.Yj10= sapply(Productlist(Y,Z),sum)/(sapply(Z,sum))*(1-A)
  est.Yj11=rep(0,J)
  est.Yj11= sapply(Productlist(Y,Z),sum)/(sapply(Z,sum))*(A)
  est.Y00= sum(est.Yj00*(1-A))/sum(1-A)
  est.Y10= sum(est.Yj10*(1-A))/sum(1-A)
  est.Y01= sum(est.Yj01*(A))/sum(A)
  est.Y11= sum(est.Yj11*(A))/sum(A)
  est.DEYj0=est.Yj10-est.Yj00
  est.DEYj1=est.Yj11-est.Yj01
  est.DEY1=est.Y11-est.Y01
  est.DEY0=est.Y10-est.Y00
  est.SEY1=est.Y11-est.Y10
  est.SEY0=est.Y01-est.Y00
  ### variance
  est.sigmaDE0=sum((est.DEYj0-est.DEY0)^2*(1-A))/(sum(1-A)-1)
  est.sigmaDE1=sum((est.DEYj1-est.DEY1)^2*(A))/(sum(A)-1)
  est.sigmaj01=rep(0,J)
  est.sigmaj00=rep(0,J)
  est.sigmaj11=rep(0,J)
  est.sigmaj10=rep(0,J)
  est.sigmab00=sum((est.Yj00-est.Y00)^2*(1-A))/(sum(1-A)-1)
  est.sigmab10=sum((est.Yj10-est.Y10)^2*(1-A))/(sum(1-A)-1)
  est.sigmab01=sum((est.Yj01-est.Y01)^2*(A))/(sum(A)-1)
  est.sigmab11=sum((est.Yj11-est.Y11)^2*(A))/(sum(A)-1)
  for (j in 1:J){
    est.sigmaj01[j]=sum((Y[[j]]-est.Yj01[j])^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.sigmaj00[j]=sum((Y[[j]]-est.Yj00[j])^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.sigmaj11[j]=sum((Y[[j]]-est.Yj11[j])^2*(Z[[j]]))/(sum(Z[[j]])-1)
    est.sigmaj10[j]=sum((Y[[j]]-est.Yj10[j])^2*(Z[[j]]))/(sum(Z[[j]])-1)
  }
  var.DEY0=est.sigmaDE0*(1/sum(1-A)-1/J)+ sum((est.sigmaj00/(n-sapply(Z,sum))+est.sigmaj10/sapply(Z,sum))*(1-A))/J/sum(1-A)
  var.DEY1=est.sigmaDE1*(1/sum(A)-1/J)+ sum((est.sigmaj01/(n-sapply(Z,sum))+est.sigmaj11/sapply(Z,sum))*(A))/J/sum(A)
  var.SEY1=est.sigmab11/sum(A)+est.sigmab10/sum(1-A)
  var.SEY0=est.sigmab01/sum(A)+est.sigmab00/sum(1-A)
  
 
  return(list( ADE1=est.DEY1,ADE0=est.DEY0,
  var.ADE1=var.DEY1,var.ADE0=var.DEY0,
  ASE1=est.SEY1,ASE0=est.SEY0,
  var.ASE1=var.SEY1,var.ASE0=var.SEY0
  ))
}
