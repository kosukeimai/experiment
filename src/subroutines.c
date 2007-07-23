
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include "vector.h"
#include "rand.h"

/*  The Sweep operator */
void SWP(
	 double **X,             /* The Matrix to work on */
	 int k,                  /* The row to sweep */
	 int size)               /* The dim. of X */
{
  int i,j;

  if (X[k][k] < 10e-20) 
    error("SWP: singular matrix.\n");
  else
    X[k][k]=-1/X[k][k];
  for(i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      if(i!=k && j!=k)
	X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];
  
}


/* inverting a matrix */
void dinv(double **X,
	  int	size,
	  double **X_inv)
{
  int i,j, k, errorM;
  double *pdInv = doubleArray(size*size);

  for (i = 0, j = 0; j < size; j++) 
    for (k = 0; k <= j; k++) 
      pdInv[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM);
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM);
    if (errorM) {
      Rprintf("LAPACK dpptri failed, %d\n", errorM);
      error("Exiting from dinv().\n");
    }
  }
  else {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dinv().\n");
  }
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      X_inv[j][k] = pdInv[i];
      X_inv[k][j] = pdInv[i++];
    }
  }

  free(pdInv);
}


/* Cholesky decomposition */
/* returns lower triangular matrix */
void dcholdc(double **X, int size, double **L)
{
  int i, j, k, errorM;
  double *pdTemp = doubleArray(size*size);

  for (j = 0, i = 0; j < size; j++) 
    for (k = 0; k <= j; k++) 
      pdTemp[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM);
  if (errorM) {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dcholdc().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
	L[j][k] = 0.0;
      else
	L[j][k] = pdTemp[i++];
    }
  }

  free(pdTemp);
} 

/* calculate the determinant of the positive definite symmetric matrix
   using the Cholesky decomposition  */
double ddet(double **X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);
  
  dcholdc(X, size, pdTemp);
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}
