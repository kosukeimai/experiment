
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void bprobit(int *Y,         /* binary outcome variable */ 
	     int *R,         /* recording indicator for Y */
	     int *D,         /* treatment status */ 
	     double *dXo,    /* covariates */
	     double *dXr,    /* covariates */
	     double *beta,   /* coefficients */
	     double *delta,  /* coefficients */
	     int *insamp,    /* # of obs */ 
	     int *incov,     /* # of covariates */
	     double *beta0,  /* prior mean */
	     double *delta0, /* prior mean */
	     double *dAo,    /* prior precision */
	     double *dAr,    /* prior precision */
	     int *param,     /* store parameters? */ 
	     int *mda,       /* marginal data augmentation? */ 
	     int *ndraws,    /* # of gibbs draws */
	     int *iBurnin,   /* # of burnin */
	     int *iKeep,     /* every ?th draws to keep */
	     int *verbose,  
	     double *coefo,  /* storage for coefficients */ 
	     double *coefr,  /* storage for coefficients */ 
	     double *ATE     /* storage for ATE */
	     ) {
  
  /*** counters ***/
  int n_samp = *insamp;      /* sample size */
  int n_gen = *ndraws;       /* number of gibbs draws */
  int n_cov = *incov;        /* number of covariates */
  
  /*** data ***/
  /* covariates for the response model */
  double **Xr = doubleMatrix(n_samp+n_cov, n_cov+1);
  /* covariates for the outcome model */     
  double **Xo = doubleMatrix(n_samp+n_cov, n_cov+1);
  double *W = doubleArray(n_samp); /* latent variable */

  /*** model parameters ***/
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta and delta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta and delta */
  double **Ao = doubleMatrix(n_cov, n_cov);
  double **Ar = doubleMatrix(n_cov, n_cov);
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, main_loop;  
  int itemp, itemp1, itemp2, itempP = ftrunc((double) n_gen/10);
  double dtemp, dtemp0, dtemp1, p0, p1, r0, r1;

  /*** marginal data augmentation ***/
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;

  /*** get random seed **/
  GetRNGstate();

  /*** read the data ***/
  itemp = 0;
  for (j = 0; j < n_cov; j++)
    for (i = 0; i < n_samp; i++) 
      Xo[i][j] = dXo[itemp++];
  itemp = 0;
  for (j = 0; j < n_cov; j++)
    for (i = 0; i < n_samp; i++) 
      Xr[i][j] = dXr[itemp++];
  
  /*** read the prior and it as additional data points ***/ 
  itemp = 0;
  for (k = 0; k < n_cov; k++)
    for (j = 0; j < n_cov; j++)
      Ao[j][k] = dAo[itemp++];

  itemp = 0;
  for (k = 0; k < n_cov; k++)
    for (j = 0; j < n_cov; j++)
      Ar[j][k] = dAr[itemp++];

  dcholdc(Ao, n_cov, mtemp);
  for(i = 0; i < n_cov; i++) {
    Xo[n_samp+i][n_cov] = 0;
    for(j = 0; j < n_cov; j++) {
      Xo[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
      Xo[n_samp+i][j] = mtemp[i][j];
    }
  }

  dcholdc(Ar, n_cov, mtemp);
  for(i = 0; i < n_cov; i++) {
    Xr[n_samp+i][n_cov] = 0;
    for(j = 0; j < n_cov; j++) {
      Xr[n_samp+i][n_cov] += mtemp[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtemp[i][j];
    }
  }

  /*** Gibbs Sampler! ***/
  itemp = 0; itemp1 = 0; itemp2 = 0;     
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /** Response Model **/    
    if (*mda) sig2 = s0/rchisq((double)nu0);

    for(i = 0; i < n_samp; i++){
      dtemp = 0;
      for(j = 0; j < n_cov; j++) 
	dtemp += Xr[i][j]*delta[j];
      if(R[i]==0) 
	W[i] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+1000,dtemp,1,0);
      Xr[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0; i < n_samp+n_cov; i++)
      for(j = 0; j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += Xr[i][j]*Xr[i][k];
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);
    /* draw delta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (*mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) 
	V[j][k] = -SS[j][k]*sig2;
    rMVN(delta, mean, V, n_cov);
    /* rescale the parameters */
    if (*mda) 
      for (i = 0; i < n_cov; i++) delta[i] /= sqrt(sig2);

    /** Imputing the missing data **/
    for (i = 0; i < n_samp; i++) {
      if (R[i] == 0) {
	p0 = beta[0];
	p1 = beta[1];
	r1 = delta[1];
	r0 = delta[0];
	for (j = 2; j < n_cov; j++) {
	  p0 += Xo[i][j]*beta[j];
	  p1 += Xo[i][j]*beta[j];
	  r0 += Xr[i][j]*delta[j];
	  r1 += Xr[i][j]*delta[j];
	}
	p0 = pnorm(p0, 0, 1, 1, 0);
	p1 = pnorm(p1, 0, 1, 1, 0);
	r0 = pnorm(r0, 0, 1, 1, 0);
	r1 = pnorm(r1, 0, 1, 1, 0);
	if (D[i] == 0) 
	  if (unif_rand() < (1-r1)*p0/((1-r1)*p0+(1-r0)*(1-p0)))
	    Y[i] = 1;
	  else
	    Y[i] = 0;
	else
	  if (unif_rand() < (1-r1)*p1/((1-r1)*p1+(1-r0)*(1-p1)))
	    Y[i] = 1;
	  else
	    Y[i] = 0;
	if (Y[i] == 0) {
	  Xr[i][0] = 1;
	  Xr[i][1] = 0;
	}
	else {
	  Xr[i][0] = 0;
	  Xr[i][1] = 1;
	}
      }
    }
      
    /** Outcome Model **/
    if (*mda) sig2 = s0/rchisq((double)nu0);
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_cov; j++) 
	dtemp += Xo[i][j]*beta[j]; 
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+1000,dtemp,1,0);
      Xo[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }
    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += Xo[i][j]*Xo[i][k];
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);

    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (*mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mean, V, n_cov); 
    /* rescaling the parameters */
    if(*mda) 
      for (i = 0; i < n_cov; i++) beta[i] /= sqrt(sig2);

    /** Compute quantities of interest **/
    p0 = 0; p1 = 0;
    for (i = 0; i < n_samp; i++) {
      dtemp0 = beta[0];
      dtemp1 = beta[1];
      for (j = 2; j < n_cov; j++) {
	dtemp0 += Xo[i][j]*beta[j];
	dtemp1 += Xo[i][j]*beta[j];
      }
      p0 += pnorm(dtemp0, 0, 1, 1, 0);
      p1 += pnorm(dtemp1, 0, 1, 1, 0);
    }
    p0 /= n_samp; p1 /= n_samp;
    
    /** Storing the results **/
    if (main_loop > *iBurnin) {
      if (keep == *iKeep) {
	ATE[itemp++] = p0;
	ATE[itemp++] = p1;
	ATE[itemp++] = p1-p0;
	if (*param) {
	  for (i = 0; i < n_cov; i++) {
	    coefo[itemp1++] = beta[i];
	    coefr[itemp2++] = delta[i];
	  }
	}
	keep = 1;
      }
      else
	keep++;
    }

    if(*verbose) {
      if(main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	itempP += ftrunc((double) n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  FreeMatrix(Xr, n_samp+n_cov);
  FreeMatrix(Xo, n_samp+n_cov);
  free(W);
  FreeMatrix(SS, n_cov+1);
  free(mean);
  FreeMatrix(V, n_cov);
  FreeMatrix(Ao, n_cov);
  FreeMatrix(Ar, n_cov);
  FreeMatrix(mtemp, n_cov);

} /* main */



