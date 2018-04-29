###
### Randomization-based method for the complier average direct effect and the complier average spillover effect
### 
###



#' Randomization-based method for the complier average direct effect and the complier average spillover effect
#'
#' 
#' This function computes the point estimates and variance estimates of the complier average direct effect (CADE)  and the complier average spillover effect (CASE).
#' The estimators calculated using this function are either individual weighted or cluster-weighted. The point estimates and variances of ITT effects are also included. 
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @param individual  A binary variable with TRUE for  individual-weighted estimators and FALSE for cluster-weighted estimators.
#' @return A list of class \code{CADErand} which contains the following items:
#' \item{CADE1}{ The point estimate of CADE(1).  } \item{CADE0}{ The point estimate of CADE(0).  } 
#'\item{CADE1}{ The point estimate of CASE(1).  } \item{CASE0}{ The point estimate of CASE(0).  } 
#'\item{var.CADE1}{ The  variance estimate of CADE(1).   } 
#'\item{var.CADE0}{ The  variance estimate of CADE(0).   } 
#'\item{var.CASE1}{ The  variance estimate of CASE(1).   } 
#'\item{var.CASE0}{ The  variance estimate of CASE(0).   } 
#'\item{DEY1}{ The point estimate of DEY(1).  } \item{DEY0}{ The point estimate of DEY(0).  } 
#'\item{DED1}{ The point estimate of DED(1).  } \item{DED0}{ The point estimate of DED(0).  } 
#'\item{var.DEY1}{ The  variance estimate of DEY(1).   } 
#'\item{var.DEY0}{ The  variance estimate of DEY(0).   } 
#'\item{var.DED1}{ The  variance estimate of DED(1).   } 
#'\item{var.DED0}{ The  variance estimate of DED(0).   } 
#'\item{SEY1}{ The point estimate of SEY(1).  } \item{SEY0}{ The point estimate of SEY(0).  } 
#'\item{SED1}{ The point estimate of SED(1).  } \item{SED0}{ The point estimate of SED(0).  } 
#'\item{var.SEY1}{ The  variance estimate of SEY(1).   } 
#'\item{var.SEY0}{ The  variance estimate of SEY(0).   } 
#'\item{var.SED1}{ The  variance estimate of SED(1).   } 
#'\item{var.SED0}{ The  variance estimate of SED(0).   } 
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@Princeton.Edu}, \url{https://imai.princeton.edu};
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai, Zhichao Jiang and Anup Malani (2018).
#' \dQuote{Causal Inference with Interference and Noncompliance in the Two-Stage Randomized Experiments}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords two-stage randomized experiments
#' @export CADErand


CADErand=function(data,individual=1){
	## transform the data into list 	
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id=unique(data$id)	
  n.cluster=length(cluster.id)	
  Z=vector("list", n.cluster) 	
  D=vector("list", n.cluster) 
  Y=vector("list", n.cluster) 
  A=rep(0,n.cluster)
  for (i in 1:n.cluster){
  	Z[[i]]=as.numeric(data$Z[data$id==cluster.id[i]])
  	D[[i]]=as.numeric(data$D[data$id==cluster.id[i]])
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
      D[[i]]=D[[i]]*n[i]*J/N
    }
  }
  est.Dj00=rep(0,J)
  est.Dj00= sapply(Difflist(D,Productlist(D,Z)),sum)/(n-sapply(Z,sum))*(1-A)
  est.Dj01=rep(0,J)
  est.Dj01= sapply(Difflist(D,Productlist(D,Z)),sum)/(n-sapply(Z,sum))*(A)
  est.Dj10=rep(0,J)
  est.Dj10=sapply(Productlist(D,Z),sum)/(sapply(Z,sum))*(1-A)
  est.Dj11=rep(0,J)
  est.Dj11= sapply(Productlist(D,Z),sum)/(sapply(Z,sum))*(A)
  est.D00= sum(est.Dj00*(1-A))/sum(1-A)
  est.D10= sum(est.Dj10*(1-A))/sum(1-A)
  est.D01= sum(est.Dj01*(A))/sum(A)
  est.D11= sum(est.Dj11*(A))/sum(A)
  est.DEDj0=est.Dj10-est.Dj00
  est.DEDj1=est.Dj11-est.Dj01
  est.DED1=est.D11-est.D01
  est.DED0=est.D10-est.D00
  est.SED1=est.D11-est.D10
  est.SED0=est.D01-est.D00
  ### variance
  est.xiDE0=sum((est.DEDj0-est.DED0)^2*(1-A))/(sum(1-A)-1)
  est.xiDE1=sum((est.DEDj1-est.DED1)^2*(A))/(sum(A)-1)
  est.xij01=rep(0,J)
  est.xij00=rep(0,J)
  est.xij11=rep(0,J)
  est.xij10=rep(0,J)
  est.xib00=sum((est.Dj00-est.D00)^2*(1-A))/(sum(1-A)-1)
  est.xib10=sum((est.Dj10-est.D10)^2*(1-A))/(sum(1-A)-1)
  est.xib01=sum((est.Dj01-est.D01)^2*(A))/(sum(A)-1)
  est.xib11=sum((est.Dj11-est.D11)^2*(A))/(sum(A)-1)
  for (j in 1:J){
    est.xij01[j]=sum((D[[j]]-est.Dj01[j])^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.xij00[j]=sum((D[[j]]-est.Dj00[j])^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.xij11[j]=sum((D[[j]]-est.Dj11[j])^2*(Z[[j]]))/(sum(Z[[j]])-1)
    est.xij10[j]=sum((D[[j]]-est.Dj10[j])^2*(Z[[j]]))/(sum(Z[[j]])-1)
  }
  var.DED0=est.xiDE0*(1/sum(1-A)-1/J)+ sum((est.xij00/(n-sapply(Z,sum))+est.xij10/sapply(Z,sum))*(1-A))/J/sum(1-A)
  var.DED1=est.xiDE1*(1/sum(A)-1/J)+ sum((est.xij01/(n-sapply(Z,sum))+est.xij11/sapply(Z,sum))*(A))/J/sum(A)
  var.SED1=est.xib11/sum(A)+est.xib10/sum(1-A)
  var.SED0=est.xib01/sum(A)+est.xib00/sum(1-A)
  
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
  
  ### covariance
  est.zetaDE0=sum((est.DEYj0-est.DEY0)*(est.DEDj0-est.DED0)*(1-A))/(sum(1-A)-1)
  est.zetaDE1=sum((est.DEYj1-est.DEY1)*(est.DEDj1-est.DED1)*(A))/(sum(A)-1)
  est.zetaj01=rep(0,J)
  est.zetaj00=rep(0,J)
  est.zetaj11=rep(0,J)
  est.zetaj10=rep(0,J)
  est.zetab00=sum((est.Yj00-est.Y00)*(est.Dj00-est.D00)*(1-A))/(sum(1-A)-1)
  est.zetab10=sum((est.Yj10-est.Y10)*(est.Dj10-est.D10)*(1-A))/(sum(1-A)-1)
  est.zetab01=sum((est.Yj01-est.Y01)*(est.Dj01-est.D01)*(A))/(sum(A)-1)
  est.zetab11=sum((est.Yj11-est.Y11)*(est.Dj11-est.D11)*(A))/(sum(A)-1)
  for (j in 1:J){
    est.zetaj01[j]=sum((Y[[j]]-est.Yj01[j])*(D[[j]]-est.Dj01[j])*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.zetaj00[j]=sum((Y[[j]]-est.Yj00[j])*(D[[j]]-est.Dj00[j])*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    est.zetaj11[j]=sum((Y[[j]]-est.Yj11[j])*(D[[j]]-est.Dj11[j])*(Z[[j]]))/(sum(Z[[j]])-1)
    est.zetaj10[j]=sum((Y[[j]]-est.Yj10[j])*(D[[j]]-est.Dj10[j])*(Z[[j]]))/(sum(Z[[j]])-1)
  }
  est.zeta0=est.zetaDE0*(1/sum(1-A)-1/J)+ sum((est.zetaj00/(n-sapply(Z,sum))+est.zetaj10/sapply(Z,sum))*(1-A))/J/sum(1-A)
  est.zeta1=est.zetaDE1*(1/sum(A)-1/J)+ sum((est.zetaj01/(n-sapply(Z,sum))+est.zetaj11/sapply(Z,sum))*(A))/J/sum(A)
  est.zetab1=est.zetab11/sum(A)+est.zetab10/sum(1-A)
  est.zetab0=est.zetab01/sum(A)+est.zetab00/sum(1-A)
  
  
  #### CADE
  est.CADE1=est.DEY1/est.DED1
  est.CADE0=est.DEY0/est.DED0
  est.CASE1=est.SEY1/est.SED1
  est.CASE0=est.SEY0/est.SED0
  est.varCADE1=   (var.DEY1-2*est.CADE1*est.zeta1+est.CADE1^2*var.DED1)/est.DED1^2
  est.varCADE0=   (var.DEY0-2*est.CADE0*est.zeta0+est.CADE0^2*var.DED0)/est.DED0^2
  est.varCASE1=   (var.SEY1-2*est.CASE1*est.zetab1+est.CASE1^2*var.SED1)/est.SED1^2
  est.varCASE0=   (var.SEY0-2*est.CASE0*est.zetab0+est.CASE0^2*var.SED0)/est.SED0^2
  return(list(CADE1=est.CADE1,CADE0=est.CADE0,CASE1=est.CASE1,CASE0=est.CASE0, var.CADE1=est.varCADE1,var.CADE0=est.varCADE0,var.CASE1=est.varCASE1,var.CASE0=est.varCASE0,
  DEY1=est.DEY1,DEY0=est.DEY0,DED1=est.DED1,DED0=est.DED0,
  var.DEY1=var.DEY1,var.DEY0=var.DEY0,var.DED1=var.DED1,var.DED0=var.DED0,
  SEY1=est.SEY1,SEY0=est.SEY0,SED1=est.SED1,SED0=est.SED0,
  var.SEY1=var.SEY1,var.SEY0=var.SEY0,var.SED1=var.SED1,var.SED0=var.SED0
  ))
}
