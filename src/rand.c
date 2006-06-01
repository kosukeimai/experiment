/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai, Jordan R. Vance, and 
  David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"


/* Sample from a univariate truncated Normal distribution 
   (truncated both from above and below): choose either inverse cdf
   method or rejection sampling method. For rejection sampling, 
   if the range is too far from mu, it uses standard rejection
   sampling algorithm with exponential envelope function. */ 
double TruncNorm(
		 double lb,  /* lower bound */ 
		 double ub,  /* upper bound */
		 double mu,  /* mean */
		 double var, /* variance */
		 int invcdf  /* use inverse cdf method? */
		 ) {
  
  double z;
  double sigma = sqrt(var);
  double stlb = (lb-mu)/sigma;  /* standardized lower bound */
  double stub = (ub-mu)/sigma;  /* standardized upper bound */
  if(stlb >= stub)
    error("TruncNorm: lower bound is greater than upper bound\n");
  if (invcdf) {  /* inverse cdf method */
    z = qnorm(runif(pnorm(stlb, 0, 1, 1, 0), pnorm(stub, 0, 1, 1, 0)),
	      0, 1, 1, 0); 
  }
  else { /* rejection sampling method */
    double tol=2.0;
    double temp, M, u, exp_par;
    int flag=0;  /* 1 if stlb, stub <-tol */
    if(stub<=-tol){
      flag=1;
      temp=stub;
      stub=-stlb;
      stlb=-temp;
    }
    if(stlb>=tol){
      exp_par=stlb;
      while(pexp(stub,1/exp_par,1,0) - pexp(stlb,1/exp_par,1,0) < 0.000001) 
	exp_par/=2.0;
      if(dnorm(stlb,0,1,1) - dexp(stlb,1/exp_par,1) >=
	 dnorm(stub,0,1,1) - dexp(stub,1/exp_par,1)) 
	M=exp(dnorm(stlb,0,1,1) - dexp(stlb,1/exp_par,1));
      else
	M=exp(dnorm(stub,0,1,1) - dexp(stub,1/exp_par,1));
      do{ 
	u=unif_rand();
	z=-log(1-u*(pexp(stub,1/exp_par,1,0)-pexp(stlb,1/exp_par,1,0))
	       -pexp(stlb,1/exp_par,1,0))/exp_par;
      }while(unif_rand() > exp(dnorm(z,0,1,1)-dexp(z,1/exp_par,1))/M );  
      if(flag==1) z=-z;
    } 
    else{ 
      do z=norm_rand();
      while( z<stlb || z>stub ); 
    }
  }
  return(z*sigma + mu); 
}


/* Sample from the MVN dist */
void rMVN(                      
	  double *Sample,         /* Vector for the sample */
	  double *mean,           /* The vector of means */
	  double **Var,           /* The matrix Variance */
	  int size)               /* The dimension */
{
  int j,k;
  double **Model = doubleMatrix(size+1, size+1);
  double cond_mean;
    
  /* draw from mult. normal using SWP */
  for(j=1;j<=size;j++){       
    for(k=1;k<=size;k++)
      Model[j][k]=Var[j-1][k-1];
    Model[0][j]=mean[j-1];
    Model[j][0]=mean[j-1];
  }
  Model[0][0]=-1;
  Sample[0]=(double)norm_rand()*sqrt(Model[1][1])+Model[0][1];
  for(j=2;j<=size;j++){
    SWP(Model,j-1,size+1);
    cond_mean=Model[j][0];
    for(k=1;k<j;k++) cond_mean+=Sample[k-1]*Model[j][k];
    Sample[j-1]=(double)norm_rand()*sqrt(Model[j][j])+cond_mean;
  }
  
  FreeMatrix(Model,size+1);
}


/* Sample from a wish dist */
/* Odell, P. L. and Feiveson, A. H. ``A Numerical Procedure to Generate
   a Sample Covariance Matrix'' Journal of the American Statistical
   Association, Vol. 61, No. 313. (Mar., 1966), pp. 199-203. */

void rWish(                  
	   double **Sample,        /* The matrix with to hold the sample */
	   double **S,             /* The parameter */
	   int df,                 /* the degrees of freedom */
	   int size)               /* The dimension */
{
  int i,j,k;
  double *V = doubleArray(size);
  double **B = doubleMatrix(size, size);
  double **C = doubleMatrix(size, size);
  double **N = doubleMatrix(size, size);
  double **mtemp = doubleMatrix(size, size);
  
  for(i=0;i<size;i++) {
    V[i]=rchisq((double) df-i-1);
    B[i][i]=V[i];
    for(j=(i+1);j<size;j++)
      N[i][j]=norm_rand();
  }

  for(i=0;i<size;i++) {
    for(j=i;j<size;j++) {
      Sample[i][j]=0;
      Sample[j][i]=0;
      mtemp[i][j]=0;
      mtemp[j][i]=0;
      if(i==j) {
	if(i>0)
	  for(k=0;k<j;k++)
	    B[j][j]+=N[k][j]*N[k][j];
      }
      else { 
	B[i][j]=N[i][j]*sqrt(V[i]);
	if(i>0)
	  for(k=0;k<i;k++)
	    B[i][j]+=N[k][i]*N[k][j];
      }
      B[j][i]=B[i][j];
    }
  }
  
  dcholdc(S, size, C);
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
	mtemp[i][j]+=C[i][k]*B[k][j];
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
	Sample[i][j]+=mtemp[i][k]*C[j][k];

  free(V);
  FreeMatrix(B, size);
  FreeMatrix(C, size);
  FreeMatrix(N, size);
  FreeMatrix(mtemp, size);
}

