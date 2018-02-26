###
### Regression-based method for the complier average direct effect
### 
###



#' Regression-based method for the complier average direct effect
#'
#' 
#' This function computes the point estimates of the complier average direct effect (CADE) and four
#'  different variance estimates: the HC2 variance, the cluster-robust variance, the cluster-robust HC2
#'  variance and the variance proposed in the reference. The estimators calculated using this function
#'  are cluster-weighted, i.e., the weights are equal for each cluster. To obtain the indivudal-weighted
#'  estimators, please multiply the recieved treatment and the outcome by \code{n_jJ/N}, where
#'  \code{n_j} is the number of individuals in cluster \code{j}, \code{J} is the number of clusters and 
#'  \code{N} is the total number of individuals. 
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @return A list of class \code{CADEreg} which contains the following items:
#' \item{CADE1}{ The point estimate of CADE(1).  } \item{CADE0}{ The point estimate of CADE(0).  } 
#' \item{var1.clu}{ The cluster-robust variance of CADE(1).   } \item{var0.clu}{ The cluster-robust variance of CADE(0).  }
#'\item{var1.clu.hc2}{ The cluster-robust HC2 variance of CADE(1).   } 
#'\item{var0.clu.hc2}{ The cluster-robust HC2 variance of CADE(0).    } 
#'\item{var1.hc2}{ The  HC2 variance of CADE(1).    } 
#'\item{var0.hc2}{ The  HC2 variance of CADE(0).    } 
#'\item{var1.ind}{ The  individual-robust variance of CADE(1).    } 
#'\item{var0.ind}{ The  individual-robust variance of CADE(0).    } 
#'\item{var1.reg}{ The  proposed variance of CADE(1).    } 
#'\item{var0.reg}{ The  proposed variance of CADE(0).    } 
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu};
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai, Zhichao Jiang and Anup Malani (2018).
#' \dQuote{Causal Inference with Interference and Noncompliance in the Two-Stage Randomized Experiments}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords two-stage randomized experiments
#' @export CADEreg


CADEreg=function(data){
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
  J=length(n)
  # weights 
  W=sum(n)
  J1=sum(A)
  J0=J-J1
  n1=sapply(Z,sum)
  n0=n-n1
  
  index.l=rep(1,J)
  index.r=rep(1,J)   
  for(j in 2:J){
    index.l[j]=1+sum(n[1:(j-1)])
  } 
  for(j in 1:J){
    index.r[j]=sum(n[1:j])
  }     
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    W[index]=ifelse(A[j]==1, 1/J1,1/(J-J1)) *   ifelse(Z[[j]]==1, 1/n1[j],1/(n0[j]))
  }
  ## Design matrix in the fist stage
  A.reg= rep(0,sum(n))
  Z.reg=rep(0,sum(n))
  D.reg=rep(0,sum(n))
  Y.reg=rep(0,sum(n))
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    A.reg[index]=A[j]
    Z.reg[index]=Z[[j]]
    D.reg[index]=D[[j]]
    Y.reg[index]=Y[[j]]
  }
  
  X= cbind(A.reg, 1-A.reg,  Z.reg*A.reg, Z.reg*(1-A.reg)  )
  
  reg1s=lm(D.reg~0+X,weights=W)
  D.hat=X%*%reg1s$coefficients
  
  M= cbind(A.reg, 1-A.reg,  D.hat*A.reg, D.hat*(1-A.reg)  )
  reg2s=lm(Y.reg~0+M,weights=as.vector(W))
  res= Y.reg-cbind(A.reg, 1-A.reg,  D.reg*A.reg, D.reg*(1-A.reg)  )%*%reg2s$coefficients
  ###  variance 
  
  ## cluster robust variance
  MM=t(M)%*%diag(W)%*%M
  
  var.cluster.med=array(0,dim=c(4,4))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]	
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index],0)
      
      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0,W[index],   0, W[index]*D.hat[index])
      
      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
      
    }
  }
  var.cluster=solve(MM)%*%var.cluster.med%*%solve(MM)
  
  
  ## cluster robust hc2 variance
  MM=t(M)%*%diag(W)%*%M
  
  var.cluster.med=array(0,dim=c(4,4))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]	
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index],0)*sqrt(J1/(J1-1)) 
      
      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0,W[index],   0, W[index]*D.hat[index])*sqrt((J-J1)/(J-J1-1))
      
      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
      
    }
  }
  
  var.cluster.hc2=solve(MM)%*%var.cluster.med%*%solve(MM)
  ### individual robust hc2
  res.ind=rep(0,sum(n))
  var.ind.med=array(0,dim=c(4,4))
  for (j in 1:J){
    index= index.l[j]:index.r[j]
    adj1=   sum(res[index]*Z[[j]]/sum(Z[[j]]))
    adj0=   sum(res[index]*(1-Z[[j]])/sum(1-Z[[j]]))
    
    res.ind[index]=res[index] - ifelse(Z[[j]]==1,adj1,adj0)
  }
  
  for (j in 1:J){
    for (i in 1:n[j]){
      index=index.l[j]-1+i
      var.ind.med=var.ind.med+(M[index,])%*% t( M[index,])  *W[index]^2 * ifelse(Z.reg[index]==1, n1[j]/(n1[j]-1),(n0[j])/(n0[j]-1))*res.ind[index]^2
    }
  }
  var.ind=solve(MM)%*%var.ind.med%*%solve(MM)
  
  
  ### traditional hc2 variance
  var.hc2.med=array(0,dim=c(4,4))
  for (j in 1:J){
    for (i in 1:n[j]){
      index=index.l[j]-1+i
      if (A[j]==1){
        constant=ifelse(Z.reg[index]==1, J1*n1[j]/(J1*n1[j]-1),J1*n0[j]/(J1*n0[j]-1))
      }else{
        constant=ifelse(Z.reg[index]==1, J0*n1[j]/(J1*n1[j]-1),J0*n0[j]/(J1*n0[j]-1))
      }
      
      var.hc2.med=var.hc2.med+(M[index,])%*% t( M[index,])  *W[index]^2 * constant*res[index]^2
    }
  }
  var.hc2=solve(MM)%*%var.hc2.med%*%solve(MM)
  
  ## results
  est.CADE1=reg2s$coefficients[3]
  est.CADE0=reg2s$coefficients[4]
  var1.cluster=var.cluster[3,3]
  var0.cluster=var.cluster[4,4]
  var1.cluster.hc2=var.cluster.hc2[3,3]
  var0.cluster.hc2=var.cluster.hc2[4,4]
  var1.ind=var.ind[3,3]
  var0.ind=var.ind[4,4]
  var1.reg=(1-J1/J)*var.cluster.hc2[3,3]+(J1/J)*var.ind[3,3]
  var0.reg=(J1/J)*var.cluster.hc2[4,4]+(1-J1/J)*var.ind[4,4]
  var1.hc2=var.hc2[3,3]
  var0.hc2=var.hc2[4,4]
  return(list(CADE1=est.CADE1,CADE0=est.CADE0, var1.clu=var1.cluster,var0.clu=var0.cluster,var1.clu.hc2=var1.cluster.hc2,var0.clu.hc2=var0.cluster.hc2,   var1.ind=var1.ind,var0.ind=var0.ind,var1.reg=var1.reg,var0.reg=var0.reg,var1.hc2=var1.hc2,var0.hc2=var0.hc2))
}