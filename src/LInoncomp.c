/****

     This file contains models for randomized experiments with
     noncompliance under the assumptions of latent ignorability of
     Frangakis and Rubin (1999, Biometrika).

****/

#include <string.h>
#include <stdio.h>      
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/*
  Read the data etc.
*/

void Prep(double *dXc, double **Xc, double *dXo, double **Xo,
	  double *dXr, double **Xr, double **Xobs, int *R,
	  int n_samp, int n_obs, int n_covC, int n_covO,
	  int n_covR, int logitC, int AT, double *dA0C,
	  double **A0C, double *dA0O, double **A0O, int priorO,
	  double *dA0R, double **A0R, double *beta0, double *delta0,
	  double *gamma0, double *pC, double *pN, double *pA,
	  double *prC, double *prN, double *prA, int *acceptC,
	  int *acceptR, int n_miss){
  int i, j, k;
  int itemp;
  double **mtempC = doubleMatrix(n_covC, n_covC); 
  double **mtempO = doubleMatrix(n_covO, n_covO); 
  double **mtempR = doubleMatrix(n_covR, n_covR); 

  /*** read the data ***/
  itemp = 0;
  for (j = 0; j < n_covC; j++)
    for (i = 0; i < n_samp; i++)
      Xc[i][j] = dXc[itemp++];

  itemp = 0;
  for (j = 0; j < n_covO; j++)
    for (i = 0; i < n_samp; i++)
      Xo[i][j] = dXo[itemp++];

  itemp = 0;
  for (j = 0; j < n_covR; j++)
    for (i = 0; i < n_samp; i++)
      Xr[i][j] = dXr[itemp++];

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) {
      for (j = 0; j < n_covO; j++)
	Xobs[itemp][j] = Xo[i][j];
      itemp++;
    }

  /*** read the prior as additional data points ***/ 
  itemp = 0; 
  if ((logitC == 1) && (AT == 1))
    for (k = 0; k < n_covC*2; k++)
      for (j = 0; j < n_covC*2; j++)
	A0C[j][k] = dA0C[itemp++];
  else
    for (k = 0; k < n_covC; k++)
      for (j = 0; j < n_covC; j++)
	A0C[j][k] = dA0C[itemp++];

  itemp = 0;
  for (k = 0; k < n_covO; k++)
    for (j = 0; j < n_covO; j++)
      A0O[j][k] = dA0O[itemp++];

  itemp = 0;
  for (k = 0; k < n_covR; k++)
    for (j = 0; j < n_covR; j++)
      A0R[j][k] = dA0R[itemp++];

  if (logitC != 1) {
    dcholdc(A0C, n_covC, mtempC);
    for (i = 0; i < n_covC; i++) {
      Xc[n_samp+i][n_covC]=0;
      for (j = 0; j < n_covC; j++) {
	Xc[n_samp+i][n_covC] += mtempC[i][j]*beta0[j];
	Xc[n_samp+i][j] = mtempC[i][j];
      }
    }
  }

  if (priorO) {
    dcholdc(A0O, n_covO, mtempO);
    for (i = 0; i < n_covO; i++) {
      Xobs[n_obs+i][n_covO]=0;
      for (j = 0; j < n_covO; j++) {
	Xobs[n_obs+i][n_covO] += mtempO[i][j]*gamma0[j];
	Xobs[n_obs+i][j] = mtempO[i][j];
      }
    }
  }
  
  dcholdc(A0R, n_covR, mtempR);
  for (i = 0; i < n_covR; i++) {
    Xr[n_samp+i][n_covR]=0;
    for (j = 0; j < n_covR; j++) {
      Xr[n_samp+i][n_covR] += mtempR[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtempR[i][j];
    }
  }

  /*** starting values for probabilities ***/
  for (i = 0; i < n_samp; i++) {
    pC[i] = 1;
    pN[i] = 1;
    pA[i] = 1;
    prC[i] = 1; 
    prN[i] = 1; 
    prA[i] = 1;
  }

  for (j = 0; j < n_covC*2; j++)
    acceptC[j] = 0;
  for (j = 0; j < n_covR; j++)
    acceptR[j] = 0;

  FreeMatrix(mtempC, n_covC);
  FreeMatrix(mtempO, n_covO);
  FreeMatrix(mtempR, n_covR);
 
} /* end of Prep */


/* 
   Response function 
*/

void Response(int logitR, int *R, double **Xr, double *delta,
	      int n_samp, int n_covR, double *delta0, double **A0R,
	      double *VarR, int *acceptR, int mda, int AT,
	      int *Z, int *D, double *prC, double *prN, double *prA){
  double dtemp;
  int i, j;

  if (logitR)
    logitMetro(R, Xr, delta, n_samp, 1, n_covR, delta0, A0R, VarR,
	       1, acceptR);
  else
    bprobitGibbs(R, Xr, delta, n_samp, n_covR, 0, delta0, A0R, mda, 1);
  
  /* Compute probabilities of R = Robs */ 
  for (i = 0; i < n_samp; i++) {
    dtemp = 0;
    if (AT) { /* always-takers */
      for (j = 3; j < n_covR; j++)
	dtemp += Xr[i][j]*delta[j];
      if ((Z[i] == 0) && (D[i] == 0)) {
	if (logitR) {
	  prC[i] = R[i]/(1+exp(-dtemp-delta[1])) + 
	    (1-R[i])/(1+exp(dtemp+delta[1]));
	  prN[i] = R[i]/(1+exp(-dtemp)) + 
	    (1-R[i])/(1+exp(dtemp));
	} else {
	  prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	  prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	} 
      }
      if ((Z[i] == 1) && (D[i] == 1)) {
	if (logitR) {
	  prC[i] = R[i]/(1+exp(-dtemp-delta[0])) + 
	    (1-R[i])/(1+exp(dtemp+delta[0]));
	  prA[i] = R[i]/(1+exp(-dtemp-delta[2])) + 
	    (1-R[i])/(1+exp(dtemp+delta[2]));
	} else {
	  prC[i] = R[i]*pnorm(dtemp+delta[0], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[0], 0, 1, 0, 0);
	  prA[i] = R[i]*pnorm(dtemp+delta[2], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[2], 0, 1, 0, 0);
	}
      }
    } else { /* no always-takers */
      for (j = 2; j < n_covR; j++)
	dtemp += Xr[i][j]*delta[j];
      if (Z[i] == 0) {
	if (logitR) {
	  prC[i] = R[i]/(1+exp(-dtemp-delta[1])) +
	    (1-R[i])/(1+exp(dtemp+delta[1]));
	  prN[i] = R[i]/(1+exp(-dtemp)) +
	    (1-R[i])/(1+exp(dtemp));
	} else { 
	  prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) + 
	    (1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	  prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	}
      }
    } 
  }
} /* end of Response */


/* 
   Compliance prediction
*/

void Compliance(int logitC, int AT, int *C, double **Xc,
		double *betaC, int n_samp, int n_covC, double *beta0,
		double **A0C, double *betaA, double *VarC, int *acceptC,
		int mda, int *A){
  int i, j; 
  int itemp;
  /* subset of the data */
  int *Atemp = intArray(n_samp);
  double **Xtemp = doubleMatrix(n_samp+n_covC, n_covC+1);
  
  if (logitC) 
    if (AT) 
      logitMetro(C, Xc, betaC, n_samp, 2, n_covC, beta0, A0C, VarC, 1,
		 acceptC); 
    else 
      logitMetro(C, Xc, betaC, n_samp, 1, n_covC, beta0, A0C, VarC, 1,
		 acceptC);  
  else {
    /* complier vs. noncomplier */
    bprobitGibbs(C, Xc, betaC, n_samp, n_covC, 0, beta0, A0C,
		 mda, 1);
    if (AT){
      /* never-taker vs. always-taker */
      /* subset the data */
      itemp = 0;
      for (i = 0; i < n_samp; i++)
	if (C[i] == 0) {
	  Atemp[itemp] = A[i];
	  for (j = 0; j < n_covC; j++)
	    Xtemp[itemp][j] = Xc[i][j];
	  itemp++;
	}
      for (i = n_samp; i < n_samp + n_covC; i++) {
	for (j = 0; j <= n_covC; j++)
	  Xtemp[itemp][j] = Xc[i][j];
	itemp++;
      }
      bprobitGibbs(Atemp, Xtemp, betaA, itemp-n_covC, n_covC, 0,
		   beta0, A0C, mda, 1); 
    }      
  }
  free(Atemp);
  FreeMatrix(Xtemp, n_samp+n_covC);
}

/* 
   Sampling Compliance Status 
*/

void SampleComp(int n_samp, int n_covC,	int AT,	double **Xc,
		double **Xo, double **Xr, double **Xobs, double *betaC,
		double *betaA, int logitC, double *qC, double *qN,
		int *Z,	int *D,	int *R,	int *RD, int *C, int *A,
		double *pC, double *pN,	double *pA, double *prA,
		double *prN, double *prC){

  int i, j, itemp;
  double dtemp, dtemp1 , dtemp2;
  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);
  double *meana = doubleArray(n_samp);

  itemp = 0;
  for (i = 0; i < n_samp; i++) {
    meanc[i] = 0;
    for (j = 0; j < n_covC; j++) 
      meanc[i] += Xc[i][j]*betaC[j];
    if (AT) { /* some always-takers */
      meana[i] = 0;
      if (logitC) { /* if logistic regression is used */
	for (j = 0; j < n_covC; j++) 
	  meana[i] += Xc[i][j]*betaC[j+n_covC];
	qC[i] = exp(meanc[i])/(1 + exp(meanc[i]) + exp(meana[i]));
	qN[i] = 1/(1 + exp(meanc[i]) + exp(meana[i]));
      } else { /* double probit regressions */
	for (j = 0; j < n_covC; j++) 
	  meana[i] += Xc[i][j]*betaA[j];
	qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	qN[i] = (1-qC[i])*pnorm(meana[i], 0, 1, 0, 0);
      }
      if (RD[i] == 0) {
	if (R[i] == 1) {
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	  dtemp1 = (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]) / 
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	} else { 
	  dtemp = qC[i]*prC[i] / 
	    (qC[i]*prC[i]+qN[i]*prN[i]+(1-qC[i]-qN[i])*prA[i]);
	  dtemp1 = (qC[i]*prC[i]+qN[i]*prN[i]) / 
	    (qC[i]*prC[i]+qN[i]*prN[i]+(1-qC[i]-qN[i])*prA[i]);
	}
	dtemp2 = unif_rand();
	if (dtemp2 < dtemp) {
	  C[i] = 1; A[i] = 0; D[i] = Z[i];
	  Xo[i][1-Z[i]] = 1; Xr[i][1-Z[i]] = 1;
	  Xo[i][Z[i]] = 0; Xr[i][Z[i]] = 0;
	  Xo[i][2] = 0; Xr[i][2] = 0;
	} else if (dtemp2 < dtemp1) {
	  C[i] = 0; A[i] = 0; D[i] = 0; 
	  Xo[i][0] = 0; Xr[i][0] = 0;
	  Xo[i][1] = 0; Xr[i][1] = 0;
	  Xo[i][2] = 0; Xr[i][2] = 0;
	}  else {
	  if (logitC)
	    C[i] = 2; 
	  else
	    C[i] = 0;
	  A[i] = 1; D[i] = 1; 
	  Xo[i][0] = 0; Xr[i][0] = 0;
	  Xo[i][1] = 0; Xr[i][1] = 0;
	  Xo[i][2] = 1; Xr[i][2] = 1;
	}
      } else if ((Z[i] == 0) && (D[i] == 0)){
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]);
	else 
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+qN[i]*prN[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	} else {
	  C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	}  
      } else if ((Z[i] == 1) && (D[i] == 1)){
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	else
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i]-qN[i])*prA[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; Xo[i][0] = 1; Xr[i][0] = 1;
	  A[i] = 0; Xo[i][2] = 0; Xr[i][2] = 0;
	} else {
	  if (logitC)
	    C[i] = 2;
	  else
	    C[i] = 0; 
	  A[i] = 1; Xo[i][0] = 0; Xr[i][0] = 0; Xo[i][2] = 1; Xr[i][2] = 1;
	}  
      }
      if (R[i] == 1) {
	Xobs[itemp][0] = Xo[i][0];
	Xobs[itemp][1] = Xo[i][1];
	Xobs[itemp][2] = Xo[i][2];
      }
    } else { /* no always-takers */
      if ((Z[i] == 0) || (RD[i] == 0)) {
	if (logitC)
	  qC[i] = 1/(1+exp(-meanc[i]));
	else
	  qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+(1-qC[i])*pN[i]*prN[i]);
	else
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i])*prN[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; D[i] = Z[i];
	  Xo[i][1-Z[i]] = 1; Xr[i][1-Z[i]] = 1; 
	  Xo[i][Z[i]] = 0; Xr[i][Z[i]] = 0; 
	} else {
	  C[i] = 0; D[i] = 0;
	  Xo[i][0] = 0; Xr[i][0] = 0;
	  Xo[i][1] = 0; Xr[i][1] = 0;
	}
      }
      if (R[i] == 1) {
	Xobs[itemp][0] = Xo[i][0]; 
	Xobs[itemp][1] = Xo[i][1]; 
      }
    }
    if (R[i] == 1) itemp++;
  }
  
  free(meanc);  
  free(meana);
}


/* 
   Binary outcomes (logit and probit)
*/

void LIbinary(int *Y,         /* binary outcome variable */ 
	      int *R,         /* recording indicator for Y */
	      int *Z,         /* treatment assignment */
	      int *D,         /* treatment status */ 
	      int *RD,        /* recording indicator for D */
	      int *C,         /* compliance status; 
				 for probit, complier = 1,
				 noncomplier = 0
				 for logit, never-taker = 0,
				 complier = 1, always-taker = 2
			      */
	      int *A,         /* always-takers; always-taker = 1, others
				 = 0 */
	      int *Ymiss,     /* number of missing obs in Y */
	      int *AT,        /* Are there always-takers? */
	      int *Insample,  /* Insample (=1) or population QoI? */
	      double *dXc,    /* model matrix for compliance model */
	      double *dXo,    /* model matrix for outcome model */
	      double *dXr,    /* model matrix for response model */
	      double *betaC,  /* coefficients for compliance model */
	      double *betaA,  /* coefficients for always-takers model */
	      double *gamma,  /* coefficients for outcome model */
	      double *delta,  /* coefficients for response model */
	      int *in_samp,   /* # of observations */
	      int *n_gen,     /* # of Gibbs draws */
	      int *in_covC,   /* # of covariates for compliance model */ 
	      int *in_covO,   /* # of covariates for outcome model */
	      int *in_covR,   /* # of covariates for response model */
	      double *beta0,  /* prior mean for betaC and betaA */ 
	      double *gamma0, /* prior mean for gamma */
	      double *delta0, /* prior mean for delta */
	      double *dA0C,   /* prior precision for betaC and betaA */ 
	      double *dA0O,   /* prior precision for gamma */
	      double *dA0R,   /* prior precision for delta */
	      double *VarC,    /* proposal variance for compliance
				  model */
	      double *VarO,   /* proposal variance for outcome model */
	      double *VarR,   /* proposal variance for response model */
	      int *logitC,    /* Use logistic regression for the
				 compliance model? */
	      int *logitO,    /* Use logistic regression for the
				 outcome model? */
	      int *logitR,    /* Use logistic regression for the
				 response model? */
	      int *param,     /* Want to keep paramters? */
	      int *mda,       /* Want to use marginal data
				 augmentation for probit regressions? */
	      int *burnin,   /* number of burnin */
	      int *iKeep,     /* keep ?th draws */
	      int *verbose,   /* print out messages */
	      double *coefC,  /* Storage for coefficients of the
				 compliance model */
	      double *coefA,  /* Storage for coefficients of the
				 always-takers model */
	      double *coefO,  /* Storage for coefficients of the
				 outcome model */
	      double *coefR,  /* Storage for coefficients of the
				 response model */	      
	      double *QoI     /* Storage of quantities of interest */
	      ) {
  /** counters **/
  int n_samp = *in_samp;
  int n_covC = *in_covC;
  int n_covO = *in_covO;
  int n_covR = *in_covR;
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  /*** observed Y ***/
  int *Yobs = intArray(n_obs);
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */
  double **Xo = doubleMatrix(n_samp+n_covO, n_covO+1);    
  /* covariates for the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_covO, n_covO+1);    
  /* covariates for the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_covR, n_covR+1);    
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);
  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  /* prior precision matrices */
  double **A0C = doubleMatrix(n_covC*2, n_covC*2);
  double **A0O = doubleMatrix(n_covO, n_covO);
  double **A0R = doubleMatrix(n_covR, n_covR);
  
  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int *acceptC = intArray(n_covC*2);      /* number of acceptance */
  int *acceptO = intArray(n_covO);      /* number of acceptance */
  int *acceptR = intArray(n_covR);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  double dtemp, dtemp1;

  /*** get random seed ***/
  GetRNGstate();

  /*** Preparing ***/
  Prep(dXc, Xc, dXo, Xo, dXr, Xr, Xobs, R, n_samp, n_obs, n_covC,
       n_covO, n_covR, *logitC, *AT, dA0C, A0C, dA0O, A0O, 1, dA0R, 
       A0R, beta0, delta0, gamma0, pC, pN, pA, prC, prN, prA, acceptC,
       acceptR, n_miss); 

  /*** observed Y ***/
  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];
  
  for (j = 0; j < n_covO; j++)
    acceptO[j] = 0;

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      Response(*logitR, R, Xr, delta, n_samp, n_covR, delta0, A0R, VarR,
	       acceptR, *mda, *AT, Z, D, prC, prN, prA);

    /** Step 2: COMPLIANCE MODEL **/
    Compliance(*logitC, *AT, C, Xc, betaC, n_samp, n_covC, beta0, A0C, 
	       betaA, VarC, acceptC, *mda, A);

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
    SampleComp(n_samp, n_covC, *AT, Xc, Xo, Xr, Xobs, betaC, betaA,
	       *logitC, qC, qN, Z, D, R, RD, C, A, pC, pN, pA, prA, prN, prC);

    /** Step 4: OUTCOME MODEL **/
    if (*logitO)
      logitMetro(Yobs, Xobs, gamma, n_obs, 1, n_covO, gamma0, A0O,
		 VarO, 1, acceptO);
    else
      bprobitGibbs(Yobs, Xobs, gamma, n_obs, n_covO, 0, gamma0, A0O, *mda, 1);

    /** Compute probabilities of Y = Yobs **/
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (*logitO) {
	      pC[i] = Y[i]/(1+exp(-meano[i]-gamma[1-Z[i]])) + (1-Y[i])/(1+exp(meano[i]+gamma[1-Z[i]]));
	      pA[i] = Y[i]/(1+exp(-meano[i]-gamma[2])) + (1-Y[i])/(1+exp(meano[i]+gamma[2]));
	      pN[i] = Y[i]/(1+exp(-meano[i])) + (1-Y[i])/(1+exp(meano[i]));
	    } else {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pA[i] = Y[i]*pnorm(meano[i]+gamma[2], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	      pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	    }
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (*logitO) {
	      pC[i] = Y[i]/(1+exp(-meano[i]-gamma[1-Z[i]])) + (1-Y[i])/(1+exp(meano[i]+gamma[1-Z[i]]));
	      pN[i] = Y[i]/(1+exp(-meano[i])) + (1-Y[i])/(1+exp(meano[i]));
	    } else {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0) + 
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	    }
	  }
      } 
    }
    
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	n_comp[0] = 0; n_comp[1] = 0; 
	n_never[0] = 0; n_never[1] = 0;
	n_always[0] = 0; n_always[1] = 0;
	p_comp = 0; p_never = 0; ITT = 0;
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i]; p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	  } else if (A[i] == 1) {
	    if (Z[i] == 1)
	      n_always[1]++;
	    else
	      n_always[0]++;
	  } else {
	    if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	  }
	  if (*Insample) { /* insample QoI */
	    if (C[i] == 1) {
	      if (*logitO) {
		dtemp = (double)(1/(1+exp(-meano[i]-gamma[0])) > unif_rand());
		dtemp1 = (double)(1/(1+exp(-meano[i]-gamma[1])) > unif_rand());
	      } else {
		dtemp = (double)((meano[i]+gamma[0]+norm_rand()) > 0);
		dtemp1 = (double)((meano[i]+gamma[1]+norm_rand()) > 0);
	      }
	      if (R[i] == 1) {
		if (Z[i] == 1) {
		  dtemp = (double)Y[i];
		} else { 
		  dtemp1 = (double)Y[i];
		}
	      }
	      Y1barC += dtemp; Y0barC += dtemp1;
	    } else if (A[i] == 1) {
	      if (R[i] == 1)
		YbarA += (double)Y[i];
	      else {
		if (*logitO)
		  YbarA += (double)(1/(1+exp(-meano[i]-gamma[2])) > unif_rand());
		else
		  YbarA += (double)((meano[i]+gamma[2]+norm_rand()) > 0);
	      }
	    } else {
	      if (R[i] == 1)
		YbarN += (double)Y[i];
	      else {
		if (*logitO)
		  YbarN += (double)(1/(1+exp(-meano[i])) > unif_rand());
		else
		  YbarN += (double)((meano[i]+norm_rand()) > 0);
	      }
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    if (*logitO) {
	      dtemp = 1/(1+exp(-meano[i]-gamma[0]));
	      dtemp1 = 1/(1+exp(-meano[i]-gamma[1])); 
	    } else {
	      dtemp = pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
	      dtemp1 = pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	    }
	    Y1barC += dtemp; Y0barC += dtemp1;
	    ITT += ((dtemp-dtemp1)*qC[i]);
	    /* always-takers */
	    if (*logitO)
	      YbarA += 1/(1+exp(-meano[i]-gamma[2]));
	    else
	      YbarA += pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	    /* never-takers */
	    if (*logitO)
	      YbarN += 1/(1+exp(-meano[i]));
	    else
	      YbarN += pnorm(meano[i], 0, 1, 1, 0);
	  }
	}

	p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
				      a complier */ 
	p_never /= (double)n_samp; /* Prob. of being a never-taker */
	if (*Insample) {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  YbarN /= (double)(n_never[0]+n_never[1]);
	  if (*AT)
	    YbarA /= (double)(n_always[0]+n_always[1]);
	} else {
	  ITT /= (double)n_samp;     /* ITT effect */
	  Y1barC /= (double)n_samp;
	  Y0barC /= (double)n_samp;
	  YbarN /= (double)n_samp;
	  if (*AT)
	    YbarA /= (double)n_samp;
	}
	CACE = Y1barC-Y0barC;    /* CACE */

	QoI[itempQ++] = ITT;   
	QoI[itempQ++] = CACE;   
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	QoI[itempQ++] = Y1barC;
	QoI[itempQ++] = Y0barC;
	QoI[itempQ++] = YbarN;
	if (*AT)
	  QoI[itempQ++] = YbarA;

	if (*param) {
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaC[j+n_covC];
	    } else {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_covR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_covC*2; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_covC; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	}
	if (*logitO) {
	  Rprintf("\n  Current Acceptance Ratio for the outcome model:");
	  for (j = 0; j < n_covO; j++)
	    Rprintf("%10g", (double)acceptO[j]/(double)main_loop);
	  Rprintf("\n");
	}
	if (*logitR) {
	  Rprintf("\n  Current Acceptance Ratio for the response model:");
	  for (j = 0; j < n_covR; j++)
	    Rprintf("%10g", (double)acceptR[j]/(double)main_loop);
	  Rprintf("\n");
	}
	itempP += ftrunc((double) *n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp+n_covO);
  FreeMatrix(Xobs, n_obs+n_covO);
  FreeMatrix(Xr, n_samp+n_covR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_covC*2);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(A0R, n_covR);
  free(acceptC);
  free(acceptO);
  free(acceptR);
  free(n_comp);
  free(n_never);
  free(n_always);

} /* end of LIbinary */


/*
  Gausssian outcome
*/

void LIgaussian(double *Y,      /* gaussian outcome variable */ 
		int *R,         /* recording indicator for Y */
		int *Z,         /* treatment assignment */
		int *D,         /* treatment status */ 
		int *RD,        /* recording indicator for D */
		int *C,         /* compliance status; 
				   for probit, complier = 1,
				 noncomplier = 0
				 for logit, never-taker = 0,
				 complier = 1, always-taker = 2
				*/
		int *A,         /* always-takers; always-taker = 1, others
				   = 0 */
		int *Ymiss,     /* number of missing obs in Y */
		int *AT,        /* Are there always-takers? */
		int *Insample,  /* Insample (=1) or population QoI? */
		double *dXc,    /* model matrix for compliance model */
		double *dXo,    /* model matrix for outcome model */
		double *dXr,    /* model matrix for response model */
		double *betaC,  /* coefficients for compliance model */
		double *betaA,  /* coefficients for always-takers model */
		double *gamma,  /* coefficients for outcome model */
		double *sig2,   /* variance for outcome model */
		double *delta,  /* coefficients for response model */
		int *in_samp,   /* # of observations */
		int *n_gen,     /* # of Gibbs draws */
		int *in_covC,   /* # of covariates for compliance model */ 
		int *in_covO,   /* # of covariates for outcome model */
		int *in_covR,   /* # of covariates for response model */
		double *beta0,  /* prior mean for betaC and betaA */ 
		double *gamma0, /* prior mean for gamma */
		double *delta0, /* prior mean for delta */
		double *dA0C,   /* prior precision for betaC and betaA */ 
		double *dA0O,   /* prior precision for gamma */
		double *dA0R,   /* prior precision for delta */
		int *nu0,       /* prior df for sig2 */
		double *s0,     /* prior scale for sig2 */
		double *VarC,    /* proposal variance for compliance
				    model */
		double *VarR,   /* proposal variance for response model */
		int *logitC,    /* Use logistic regression for the
				   compliance model? */
		int *logitR,    /* Use logistic regression for the
				   response model? */
		int *param,     /* Want to keep paramters? */
		int *mda,       /* Want to use marginal data
				   augmentation for probit regressions? */
		int *burnin,   /* number of burnin */
		int *iKeep,     /* keep ?th draws */
		int *verbose,   /* print out messages */
		double *coefC,  /* Storage for coefficients of the
				   compliance model */
		double *coefA,  /* Storage for coefficients of the
				   always-takers model */
		double *coefO,  /* Storage for coefficients of the
				   outcome model */
		double *coefR,  /* Storage for coefficients of the
				   response model */	      
		double *var,    /* Storage for sig2 */
		double *QoI     /* Storage of quantities of interest */
	      ) {
  /** counters **/
  int n_samp = *in_samp;
  int n_covC = *in_covC;
  int n_covO = *in_covO;
  int n_covR = *in_covR;
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */
  double **Xo = doubleMatrix(n_samp+n_covO, n_covO+1);    
  /* covariates for the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_covO, n_covO+1);    
  /* covariates for the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_covR, n_covR+1);    
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);
  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  /* prior precision matrices */
  double **A0C = doubleMatrix(n_covC*2, n_covC*2);
  double **A0O = doubleMatrix(n_covO, n_covO);
  double **A0R = doubleMatrix(n_covR, n_covR);
  
  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int *acceptC = intArray(n_covC*2);      /* number of acceptance */
  int *acceptR = intArray(n_covR);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR, itempS;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Preparing ***/
  Prep(dXc, Xc, dXo, Xo, dXr, Xr, Xobs, R, n_samp, n_obs, n_covC,
       n_covO, n_covR, *logitC, *AT, dA0C, A0C, dA0O, A0O, 1, dA0R, 
       A0R, beta0, delta0, gamma0, pC, pN, pA, prC, prN, prA, acceptC,
       acceptR, n_miss); 

  /*** observed Y ***/
  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Xobs[itemp++][n_covO] = Y[i];

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0; itempS = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      Response(*logitR, R, Xr, delta, n_samp, n_covR, delta0, A0R, VarR,
	       acceptR, *mda, *AT, Z, D, prC, prN, prA);

    /** Step 2: COMPLIANCE MODEL **/
    Compliance(*logitC, *AT, C, Xc, betaC, n_samp, n_covC, beta0, A0C, 
	       betaA, VarC, acceptC, *mda, A);

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
    SampleComp(n_samp, n_covC, *AT, Xc, Xo, Xr, Xobs, betaC, betaA,
	       *logitC, qC, qN, Z, D, R, RD, C, A, pC, pN, pA, prA, prN, prC);

    /** Step 4: OUTCOME MODEL **/
    bNormalReg(Xobs, gamma, sig2, n_obs, n_covO, 0, 1, gamma0, A0O, 1,
	       *nu0, *s0, 0);

    /** Compute probabilities of Y = Yobs **/
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]], sqrt(*sig2), 0);
	    pN[i] = dnorm(Y[i], meano[i], sqrt(*sig2), 0);
	    pA[i] = dnorm(Y[i], meano[i]+gamma[2], sqrt(*sig2), 0);
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)){
	    pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]], sqrt(*sig2), 0);
	    pN[i] = dnorm(Y[i], meano[i], sqrt(*sig2), 0);
	  }
      } 
    }
    
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	n_comp[0] = 0; n_comp[1] = 0; 
	n_never[0] = 0; n_never[1] = 0;
	n_always[0] = 0; n_always[1] = 0;
	p_comp = 0; p_never = 0; ITT = 0;
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i]; p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	  } else if (A[i] == 1) {
	    if (Z[i] == 1)
	      n_always[1]++;
	    else
	      n_always[0]++;
	  } else {
	    if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	  }
	  if (*Insample) { /* insample QoI */
	    if (C[i] == 1) {
	      dtemp = rnorm(meano[i]+gamma[0], sqrt(*sig2));
	      dtemp1 = rnorm(meano[i]+gamma[1], sqrt(*sig2));
	      if (R[i] == 1) {
		if (Z[i] == 1) {
		  dtemp = Y[i];
		} else { 
		  dtemp1 = Y[i];
		}
	      }
	      Y1barC += dtemp; Y0barC += dtemp1;
	    } else if (A[i] == 1) {
	      if (R[i] == 1)
		YbarA += Y[i];
	      else 
		YbarA += rnorm(meano[i]+gamma[2], sqrt(*sig2));
	    } else {
	      if (R[i] == 1)
		YbarN += Y[i];
	      else 
		YbarN += rnorm(meano[i], sqrt(*sig2));
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    dtemp = (meano[i]+gamma[0]);
	    dtemp1 = (meano[i]+gamma[1]); 
	    Y1barC += dtemp; Y0barC += dtemp1;
	    ITT += ((dtemp-dtemp1)*qC[i]);
	    /* always-takers */
	    YbarA += (meano[i]+gamma[2]);
	    /* never-takers */
	    YbarN += meano[i];
	  }
	}

	p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
				      a complier */ 
	p_never /= (double)n_samp; /* Prob. of being a never-taker */
	if (*Insample) {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  YbarN /= (double)(n_never[0]+n_never[1]);
	  if (*AT)
	    YbarA /= (double)(n_always[0]+n_always[1]);
	} else {
	  ITT /= (double)n_samp;     /* ITT effect */
	  Y1barC /= (double)n_samp;
	  Y0barC /= (double)n_samp;
	  YbarN /= (double)n_samp;
	  if (*AT)
	    YbarA /= (double)n_samp;
	}
	CACE = Y1barC-Y0barC;    /* CACE */

	QoI[itempQ++] = ITT;   
	QoI[itempQ++] = CACE;   
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	QoI[itempQ++] = Y1barC;
	QoI[itempQ++] = Y0barC;
	QoI[itempQ++] = YbarN;
	if (*AT)
	  QoI[itempQ++] = YbarA;

	if (*param) {
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  var[itempS++] = sig2[0];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaC[j+n_covC];
	    } else {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_covR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_covC*2; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_covC; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	}
	if (*logitR) {
	  Rprintf("\n  Current Acceptance Ratio for the response model:");
	  for (j = 0; j < n_covR; j++)
	    Rprintf("%10g", (double)acceptR[j]/(double)main_loop);
	  Rprintf("\n");
	}
	itempP += ftrunc((double) *n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp+n_covO);
  FreeMatrix(Xobs, n_obs+n_covO);
  FreeMatrix(Xr, n_samp+n_covR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_covC*2);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(A0R, n_covR);
  free(acceptC);
  free(acceptR);
  free(n_comp);
  free(n_never);
  free(n_always);

} /* end of LIgaussian */



/* 
   Ordinal outcomes (probit)
*/

void LIordinal(int *Y,         /* binary outcome variable */ 
	       int *R,         /* recording indicator for Y */
	       int *Z,         /* treatment assignment */
	       int *D,         /* treatment status */ 
	       int *RD,        /* recording indicator for D */
	       int *C,         /* compliance status; 
				  for probit, complier = 1,
				  noncomplier = 0
				  for logit, never-taker = 0,
				  complier = 1, always-taker = 2
			       */
	       int *A,         /* always-takers; always-taker = 1, others
				  = 0 */
	       int *Ymiss,     /* number of missing obs in Y */
	       int *AT,        /* Are there always-takers? */
	       int *Insample,  /* Insample (=1) or population QoI? */
	       double *dXc,    /* model matrix for compliance model */
	       double *dXo,    /* model matrix for outcome model */
	       double *dXr,    /* model matrix for response model */
	       double *betaC,  /* coefficients for compliance model */
	       double *betaA,  /* coefficients for always-takers model */
	       double *gamma,  /* coefficients for outcome model */
	       double *tau,    /* J cutpoints where the first
				  cutpoint is set to 0 and last set to
				  tau_{J-1} + 1000 */
	       double *delta,  /* coefficients for response model */
	       int *in_samp,   /* # of observations */
	       int *n_gen,     /* # of Gibbs draws */
	       int *n_cat,     /* # of outcome categories: J */
	       int *in_covC,   /* # of covariates for compliance model */ 
	       int *in_covO,   /* # of covariates for outcome model */
	       int *in_covR,   /* # of covariates for response model */
	       double *beta0,  /* prior mean for betaC and betaA */ 
	       double *gamma0, /* prior mean for gamma */
	       double *delta0, /* prior mean for delta */
	       double *dA0C,   /* prior precision for betaC and betaA */ 
	       double *dA0O,   /* prior precision for gamma */
	       double *dA0R,   /* prior precision for delta */
	       double *VarC,   /* proposal variance for compliance
				  model */
	       double *VarO,   /* proposal variance for taus */
	       double *VarR,   /* proposal variance for response model */
	       int *logitC,    /* Use logistic regression for the
				  compliance model? */
	       int *logitR,    /* Use logistic regression for the
				  response model? */
	       int *param,     /* Want to keep paramters? */
	       int *mda,       /* Want to use marginal data
				  augmentation for probit regressions? */
	       int *burnin,   /* number of burnin */
	       int *iKeep,     /* keep ?th draws */
	       int *verbose,   /* print out messages */
	       double *coefC,  /* Storage for coefficients of the
				  compliance model */
	       double *coefA,  /* Storage for coefficients of the
				  always-takers model */
	       double *coefO,  /* Storage for coefficients of the
				  outcome model */
	       double *coefR,  /* Storage for coefficients of the
				  response model */	      
	       double *tauO,   /* Storage for taus */
	       double *QoI     /* Storage of quantities of interest */
	       ) {
  /** counters **/
  int n_samp = *in_samp;
  int n_covC = *in_covC;
  int n_covO = *in_covO;
  int n_covR = *in_covR;
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  /*** observed Y ***/
  int *Yobs = intArray(n_obs);
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */
  double **Xo = doubleMatrix(n_samp+n_covO, n_covO+1);    
  /* covariates for the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_covO, n_covO+1);    
  /* covariates for the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_covR, n_covR+1);    
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);
  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  /* prior precision matrices */
  double **A0C = doubleMatrix(n_covC*2, n_covC*2);
  double **A0O = doubleMatrix(n_covO, n_covO);
  double **A0R = doubleMatrix(n_covR, n_covR);
  
  /* quantities of interest: ITT, CACE  */
  double *ITT = doubleArray(*n_cat-1);
  double *CACE = doubleArray(*n_cat-1);
  double *Y1barC = doubleArray(*n_cat-1); 
  double *Y0barC = doubleArray(*n_cat-1);
  double *YbarN = doubleArray(*n_cat-1);
  double *YbarA = doubleArray(*n_cat-1);
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress; progress = 1;
  int keep; keep = 1;
  int *acceptC = intArray(n_covC*2);    /* number of acceptance */
  int *acceptR = intArray(n_covR);      /* number of acceptance */
  int *acceptO = intArray(1); acceptO[0] = 0;
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR, itempT;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Preparing ***/
  Prep(dXc, Xc, dXo, Xo, dXr, Xr, Xobs, R, n_samp, n_obs, n_covC,
       n_covO, n_covR, *logitC, *AT, dA0C, A0C, dA0O, A0O, 1, dA0R, 
       A0R, beta0, delta0, gamma0, pC, pN, pA, prC, prN, prA, acceptC,
       acceptR, n_miss); 

  /*** observed Y ***/
  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];
  
  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0; itempT = 0;  
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      Response(*logitR, R, Xr, delta, n_samp, n_covR, delta0, A0R, VarR,
	       acceptR, *mda, *AT, Z, D, prC, prN, prA);

    /** Step 2: COMPLIANCE MODEL **/    
    Compliance(*logitC, *AT, C, Xc, betaC, n_samp, n_covC, beta0, A0C, 
	       betaA, VarC, acceptC, *mda, A);

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
    SampleComp(n_samp, n_covC, *AT, Xc, Xo, Xr, Xobs, betaC, betaA,
	       *logitC, qC, qN, Z, D, R, RD, C, A, pC, pN, pA, prA, prN, prC);

    /** Step 4: OUTCOME MODEL **/
    boprobitMCMC(Yobs, Xobs, gamma, tau, n_obs, n_covO, *n_cat,
		 0, gamma0, A0O, *mda, 1, VarO, acceptO, 1);

    /** Compute probabilities of Y = 1 **/
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (Y[i] == 0) {
	      pC[i] = pnorm(tau[0],meano[i]+gamma[1-Z[i]], 1, 1, 0);
	      pA[i] = pnorm(tau[0],meano[i]+gamma[2], 1, 1, 0);
	      pN[i] = pnorm(tau[0],meano[i], 1, 1, 0);
	    } else if (Y[i] == (*n_cat-1)) {
	      pC[i] = pnorm(tau[*n_cat-2],meano[i]+gamma[1-Z[i]], 1, 0, 0);
	      pA[i] = pnorm(tau[*n_cat-2],meano[i]+gamma[2], 1, 0, 0);
	      pN[i] = pnorm(tau[*n_cat-2],meano[i], 1, 0, 0);
	    } else {
	      pC[i] = pnorm(tau[Y[i]],meano[i]+gamma[1-Z[i]], 1, 1, 0) - 
		pnorm(tau[Y[i]-1],meano[i]+gamma[1-Z[i]], 1, 1, 0);
	      pA[i] = pnorm(tau[Y[i]],meano[i]+gamma[2], 1, 1, 0) - 
		pnorm(tau[Y[i]-1],meano[i]+gamma[2], 1, 1, 0);
	      pN[i] = pnorm(tau[Y[i]],meano[i], 1, 1, 0) - 
		pnorm(tau[Y[i]-1],meano[i], 1, 1, 0);
	    }
	  }
	}
      } else { /* no always-takers */
	for (j = 2; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1) 
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (Y[i] == 0) {
	      pC[i] = pnorm(tau[0],meano[i]+gamma[1-Z[i]], 1, 1, 0);
	      pN[i] = pnorm(tau[0],meano[i], 1, 1, 0);
	    } else if (Y[i] == (*n_cat-1)) {
	      pC[i] = pnorm(tau[*n_cat-2],meano[i]+gamma[1-Z[i]], 1, 0, 0);
	      pN[i] = pnorm(tau[*n_cat-2],meano[i], 1, 0, 0);
	    } else {
	      pC[i] = pnorm(tau[Y[i]],meano[i]+gamma[1-Z[i]], 1, 1, 0) - 
		pnorm(tau[Y[i]-1],meano[i]+gamma[1-Z[i]], 1, 1, 0);
	      pN[i] = pnorm(tau[Y[i]],meano[i]+gamma[1], 1, 1, 0) - 
		pnorm(tau[Y[i]-1],meano[i], 1, 1, 0);
	    }
	  }
      }
    }

    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	n_comp[0] = 0; n_comp[1] = 0; 
	n_never[0] = 0; n_never[1] = 0;
	n_always[0] = 0; n_always[1] = 0;
	p_comp = 0; p_never = 0; 
	for (j = 0; j < *n_cat-1; j++) {
	  ITT[j] = 0;
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i]; p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	  } else if (A[i] == 1) {
	    if (Z[i] == 1)
	      n_always[1]++;
	    else
	      n_always[0]++;
	  } else {
	    if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	  }
	  if (*Insample) { /* insample QoI */
	    for (j = 1; j < *n_cat; j++) {
	      if (C[i] == 1) {
		if (j == (*n_cat-1)) { 
		  /* Y1barC[j-1] and Y0barC[j-1] */
		  dtemp = (double)(unif_rand() <
				   pnorm(tau[*n_cat-2], meano[i]+gamma[0], 1, 0, 0));
		  dtemp1 = (double)(unif_rand() <
				    pnorm(tau[*n_cat-2], meano[i]+gamma[1], 1, 0, 0));
		} else {
		  dtemp = (double)(unif_rand() < 
				   (pnorm(tau[j], meano[i]+gamma[0], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[0], 1, 1, 0)));
		  dtemp1 = (double)(unif_rand() < 
				    (pnorm(tau[j], meano[i]+gamma[1], 1, 1, 0) 
				     -pnorm(tau[j-1], meano[i]+gamma[1], 1, 1, 0)));		
		}
		if (R[i] == 1) {
		  if (Z[i] == 1) {
		    dtemp = (double)(Y[i] == j);
		  } else { 
		    dtemp1 = (double)(Y[i] == j);
		  }
		}
		Y1barC[j-1] += dtemp; Y0barC[j-1] += dtemp1;
	      } else if (A[i] == 1) {
		if (R[i] == 1)
		  YbarA[j-1] += (double)(Y[i] == j);
		else
		  if (j == (*n_cat-1))
		    YbarA[j-1] += (double)(unif_rand() <
					   pnorm(tau[*n_cat-2], meano[i]+gamma[2], 1, 0, 0));
		  else
		    YbarA[j-1] += (double)(unif_rand() < 
					   (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
					    -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0)));		
	      } else {
		if (R[i] == 1)
		  YbarN[j-1] += (double)(Y[i] == j);
		else 
		  if (j == (*n_cat-1))
		    YbarN[j-1] += (double)(unif_rand() <
					   pnorm(tau[*n_cat-2], meano[i], 1, 0, 0));
		  else
		    YbarN[j-1] += (double)(unif_rand() < 
					   (pnorm(tau[j], meano[i], 1, 1, 0) 
					    -pnorm(tau[j-1], meano[i], 1, 1, 0)));		
	      }
	    }
	  } else { /* population QoI */
	    for (j = 1; j < *n_cat; j++) {
	      /* compliers */
	      if (j == (*n_cat-1)) {
		dtemp = pnorm(tau[*n_cat-2], meano[i]+gamma[0], 1, 0, 0);
		dtemp1 = pnorm(tau[*n_cat-2], meano[i]+gamma[1], 1, 0, 0);
	      } else {
		dtemp = (pnorm(tau[j], meano[i]+gamma[0], 1, 1, 0) 
			 -pnorm(tau[j-1], meano[i]+gamma[0], 1, 1, 0));
		dtemp1 = (pnorm(tau[j], meano[i]+gamma[1], 1, 1, 0) 
			  -pnorm(tau[j-1], meano[i]+gamma[1], 1, 1, 0));
	      }
	      Y1barC[j-1] += dtemp; Y0barC[j-1] += dtemp1;
	      ITT[j-1] += ((dtemp-dtemp1)*qC[i]);
	      /* always-takers */
	      if (j == (*n_cat-1))
		YbarA[j-1] += pnorm(tau[*n_cat-2], meano[i]+gamma[2], 1, 0, 0);
	      else
		YbarA[j-1] += (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
			       -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0));
	      /* never-takers */
	      if (j == (*n_cat-1))
		YbarN[j-1] += pnorm(tau[*n_cat-2], meano[i], 1, 0, 0);
	      else
		YbarN[j-1] += (pnorm(tau[j], meano[i], 1, 1, 0) 
			       -pnorm(tau[j-1], meano[i], 1, 1, 0));
	    }
	  }
	}

	p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
				      a complier */ 
	p_never /= (double)n_samp; /* Prob. of being a never-taker */
	for (j = 0; j < (*n_cat-1); j++) {
	  if (*Insample) {
	    ITT[j] = (Y1barC[j]-Y0barC[j])/(double)n_samp;     /* ITT effect */
	    Y1barC[j] /= (double)(n_comp[0]+n_comp[1]);  
	    Y0barC[j] /= (double)(n_comp[0]+n_comp[1]); 
	    YbarN[j] /= (double)(n_never[0]+n_never[1]);
	    if (*AT)
	      YbarA[j] /= (double)(n_always[0]+n_always[1]);
	  } else {
	    ITT[j] /= (double)n_samp;     /* ITT effect */
	    Y1barC[j] /= (double)n_samp;
	    Y0barC[j] /= (double)n_samp;
	    YbarN[j] /= (double)n_samp;
	    if (*AT)
	      YbarA[j] /= (double)n_samp;
	  }
	  CACE[j] = Y1barC[j]-Y0barC[j];    /* CACE */
	}

	for (j = 0; j < (*n_cat-1); j++) 
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (*n_cat-1); j++) 
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (*n_cat-1); j++) 
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (*n_cat-1); j++) 
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (*n_cat-1); j++) 
	  QoI[itempQ++] = YbarN[j];
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	if (*AT)
	  for (j = 0; j < (*n_cat-1); j++) 
	    QoI[itempQ++] = YbarA[j];

	if (*param) {
	  for (j = 0; j < (*n_cat-1); j++)
	    tauO[itempT++] = tau[j];
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaC[j+n_covC];
	    } else {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_covR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_covC*2; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_covC; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	}
	if (*logitR) {
	  Rprintf("\n  Current Acceptance Ratio for the response model:");
	  for (j = 0; j < n_covR; j++)
	    Rprintf("%10g", (double)acceptR[j]/(double)main_loop);
	  Rprintf("\n");
	}
	Rprintf("\n  Current Acceptance Ratio for the outcome model:");
	Rprintf("%10g\n", (double)acceptO[0]/(double)main_loop);
	itempP += ftrunc((double) *n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp+n_covO);
  FreeMatrix(Xobs, n_obs+n_covO);
  FreeMatrix(Xr, n_samp+n_covR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_covC*2);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(A0R, n_covR);
  free(ITT);
  free(CACE);
  free(Y1barC);
  free(Y0barC);
  free(YbarN);
  free(YbarA);
  free(n_comp);
  free(n_never);
  free(n_always);
  free(acceptC);
  free(acceptR);
  free(acceptO);

} /* end of LIordinal */



/*
  Count outcome (negative binomial)
*/

void LIcount(int *Y,         /*count outcome variable */ 
	     int *R,         /* recording indicator for Y */
	     int *Z,         /* treatment assignment */
	     int *D,         /* treatment status */ 
	     int *RD,        /* recording indicator for D */
	     int *C,         /* compliance status; 
				for probit, complier = 1,
				noncomplier = 0
				for logit, never-taker = 0,
				complier = 1, always-taker = 2
			     */
	     int *A,         /* always-takers; always-taker = 1, others
				= 0 */
	     int *Ymiss,     /* number of missing obs in Y */
	     int *AT,        /* Are there always-takers? */
	     int *Insample,  /* Insample (=1) or population QoI? */
	     double *dXc,    /* model matrix for compliance model */
	     double *dXo,    /* model matrix for outcome model */
	     double *dXr,    /* model matrix for response model */
	     double *betaC,  /* coefficients for compliance model */
	     double *betaA,  /* coefficients for always-takers model */
	     double *gamma,  /* coefficients for outcome model */
	     double *sig2,   /* dispersion for outcome model */
	     double *delta,  /* coefficients for response model */
	     int *in_samp,   /* # of observations */
	     int *n_gen,     /* # of Gibbs draws */
	     int *in_covC,   /* # of covariates for compliance model */ 
	     int *in_covO,   /* # of covariates for outcome model */
	     int *in_covR,   /* # of covariates for response model */
	     double *beta0,  /* prior mean for betaC and betaA */ 
	     double *gamma0, /* prior mean for gamma */
	     double *delta0, /* prior mean for delta */
	     double *dA0C,   /* prior precision for betaC and betaA */ 
	     double *dA0O,   /* prior precision for gamma */
	     double *dA0R,   /* prior precision for delta */
	     double *a0,     /* prior shape for sig2 */
	     double *b0,     /* prior scale for sig2 */
	     double *VarC,   /* proposal variance for compliance
				 model */
	     double *VarR,   /* proposal variance for response model */
	     double *VarO,   /* proposal variance for outcome model */
	     double *VarS,   /* proposal variance for dispersion */
	     int *logitC,    /* Use logistic regression for the
				compliance model? */
	     int *logitR,    /* Use logistic regression for the
				response model? */
	     int *param,     /* Want to keep paramters? */
	     int *mda,       /* Want to use marginal data
				augmentation for probit regressions? */
	     int *burnin,   /* number of burnin */
	     int *iKeep,     /* keep ?th draws */
	     int *verbose,   /* print out messages */
	     double *coefC,  /* Storage for coefficients of the
				compliance model */
	     double *coefA,  /* Storage for coefficients of the
				always-takers model */
	     double *coefO,  /* Storage for coefficients of the
				outcome model */
	     double *coefR,  /* Storage for coefficients of the
				response model */	      
	     double *var,    /* Storage for sig2 */
	     double *QoI     /* Storage of quantities of interest */
	     ) {
  /** counters **/
  int n_samp = *in_samp;
  int n_covC = *in_covC;
  int n_covO = *in_covO;
  int n_covR = *in_covR;
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  /*** observed Y ***/
  int *Yobs = intArray(n_obs);
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */
  double **Xo = doubleMatrix(n_samp, n_covO);    
  /* covariates for the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs, n_covO);    
  /* covariates for the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_covR, n_covR+1);    
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /*** model parameters ***/
  double *cont = doubleArray(n_obs);
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);
  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  /* prior precision matrices */
  double **A0C = doubleMatrix(n_covC*2, n_covC*2);
  double **A0O = doubleMatrix(n_covO, n_covO);
  double **A0R = doubleMatrix(n_covR, n_covR);
  
  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int *acceptC = intArray(n_covC*2);      /* number of acceptance */
  int *acceptO = intArray(2);
  int *acceptR = intArray(n_covR);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR, itempS;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Preparing ***/
  Prep(dXc, Xc, dXo, Xo, dXr, Xr, Xobs, R, n_samp, n_obs, n_covC,
       n_covO, n_covR, *logitC, *AT, dA0C, A0C, dA0O, A0O, 0, dA0R, 
       A0R, beta0, delta0, gamma0, pC, pN, pA, prC, prN, prA, acceptC,
       acceptR, n_miss); 

  /*** observed Y ***/
  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];

  for (i = 0; i < n_obs; i++)
    cont[i] = 0;
  acceptO[0] = 0; acceptO[1] = 0;

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0; itempS = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      Response(*logitR, R, Xr, delta, n_samp, n_covR, delta0, A0R, VarR,
	       acceptR, *mda, *AT, Z, D, prC, prN, prA);

    /** Step 2: COMPLIANCE MODEL **/
    Compliance(*logitC, *AT, C, Xc, betaC, n_samp, n_covC, beta0, A0C, 
	       betaA, VarC, acceptC, *mda, A);

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
    SampleComp(n_samp, n_covC, *AT, Xc, Xo, Xr, Xobs, betaC, betaA,
	       *logitC, qC, qN, Z, D, R, RD, C, A, pC, pN, pA, prA, prN, prC);

    /** Step 4: OUTCOME MODEL **/
    negbinMetro(Yobs, Xobs, gamma, sig2, n_obs, n_covO, gamma0, A0O,
		*a0, *b0, VarO, *VarS, cont, 1, acceptO, 0);

    /** Compute probabilities of Y = 1 **/
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]), *sig2, 0);
	    pA[i] = dnegbin(Y[i], exp(meano[i]+gamma[2]), *sig2, 0);
	    pN[i] = dnegbin(Y[i], exp(meano[i]), *sig2, 0);
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]), *sig2, 0);
	    pN[i] = dnegbin(Y[i], exp(meano[i]), *sig2, 0);
	  }
      } 
    }

    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	n_comp[0] = 0; n_comp[1] = 0; 
	n_never[0] = 0; n_never[1] = 0;
	n_always[0] = 0; n_always[1] = 0;
	p_comp = 0; p_never = 0; ITT = 0;
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i]; p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	  } else if (A[i] == 1) {
	    if (Z[i] == 1)
	      n_always[1]++;
	    else
	      n_always[0]++;
	  } else {
	    if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	  }
	  if (*Insample) { /* insample QoI */
	    if (C[i] == 1) {
	      dtemp = rnegbin(exp(meano[i]+gamma[0]), *sig2);
	      dtemp1 = rnegbin(exp(meano[i]+gamma[1]), *sig2);
	      if (R[i] == 1) {
		if (Z[i] == 1) {
		  dtemp = (double)Y[i];
		} else { 
		  dtemp1 = (double) Y[i];
		}
	      }
	      Y1barC += dtemp; Y0barC += dtemp1;
	    } else if (A[i] == 1) {
	      if (R[i] == 1)
		YbarA += (double)Y[i];
	      else 
		YbarA += rnegbin(exp(meano[i]+gamma[2]), *sig2);
	    } else {
	      if (R[i] == 1)
		YbarN += (double)Y[i];
	      else 
		YbarN += rnegbin(exp(meano[i]), *sig2);
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    dtemp = exp(meano[i]+gamma[0]);
	    dtemp1 = exp(meano[i]+gamma[1]); 
	    Y1barC += dtemp; Y0barC += dtemp1;
	    ITT += ((dtemp-dtemp1)*qC[i]);
	    /* always-takers */
	    YbarA += exp(meano[i]+gamma[2]);
	    /* never-takers */
	    YbarN += exp(meano[i]);
	  }
	}

	p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
				      a complier */ 
	p_never /= (double)n_samp; /* Prob. of being a never-taker */
	if (*Insample) {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  YbarN /= (double)(n_never[0]+n_never[1]);
	  if (*AT)
	    YbarA /= (double)(n_always[0]+n_always[1]);
	} else {
	  ITT /= (double)n_samp;     /* ITT effect */
	  Y1barC /= (double)n_samp;
	  Y0barC /= (double)n_samp;
	  YbarN /= (double)n_samp;
	  if (*AT)
	    YbarA /= (double)n_samp;
	}
	CACE = Y1barC-Y0barC;    /* CACE */
    
	QoI[itempQ++] = ITT;   
	QoI[itempQ++] = CACE;   
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	QoI[itempQ++] = Y1barC;
	QoI[itempQ++] = Y0barC;
	QoI[itempQ++] = YbarN;
	if (*AT)
	  QoI[itempQ++] = YbarA;

	if (*param) {
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  var[itempS++] = sig2[0];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaC[j+n_covC];
	    } else {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_covR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_covC*2; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_covC; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	}
	if (*logitR) {
	  Rprintf("\n  Current Acceptance Ratio for the response model:");
	  for (j = 0; j < n_covR; j++)
	    Rprintf("%10g", (double)acceptR[j]/(double)main_loop);
	  Rprintf("\n");
	}
	Rprintf("\n  Current Acceptance Ratio for the outcome model:");
	Rprintf("%10g%10g\n", (double)acceptO[0]/(double)main_loop, 
		(double)acceptO[1]/(double)main_loop);
       
	itempP += ftrunc((double) *n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp);
  FreeMatrix(Xobs, n_obs);
  FreeMatrix(Xr, n_samp+n_covR);
  free(meano);
  free(cont);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_covC*2);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(A0R, n_covR);
  free(acceptC);
  free(acceptO);
  free(acceptR);
  free(n_comp);
  free(n_never);
  free(n_always);

} /* end of LIcount */


/*
  Gausssian outcome
*/

void LItwopart(int *Y,      /* indicator variable; Y > 0 */
	       double *Y1,     /* gaussian outcome variable for Y > 0 */ 
	       int *R,         /* recording indicator for Y */
	       int *Z,         /* treatment assignment */
	       int *D,         /* treatment status */ 
	       int *RD,        /* recording indicator for D */
	       int *C,         /* compliance status; 
				  for probit, complier = 1,
				  noncomplier = 0
				  for logit, never-taker = 0,
				  complier = 1, always-taker = 2
			       */
	       int *A,         /* always-takers; always-taker = 1, others
				  = 0 */
	       int *Ymiss,     /* number of missing obs in Y */
	       int *AT,        /* Are there always-takers? */
	       int *Insample,  /* Insample (=1) or population QoI? */
	       double *dXc,    /* model matrix for compliance model */
	       double *dXo,    /* model matrix for outcome model */
	       double *dXr,    /* model matrix for response model */
	       double *betaC,  /* coefficients for compliance model */
	       double *betaA,  /* coefficients for always-takers model */
	       double *gamma,  /* coefficients for outcome model 1 */
	       double *gamma1, /* coeffs for outcoem model 2 */
	       double *sig2,   /* variance for outcome model */
	       double *delta,  /* coefficients for response model */
	       int *in_samp,   /* # of observations (a vector where
				  the first element is total number of obs and the second
				  element is the number of obs with Y > 0 */
	       int *n_gen,     /* # of Gibbs draws */
	       int *in_covC,   /* # of covariates for compliance model */ 
	       int *in_covO,   /* # of covariates for outcome model */
	       int *in_covR,   /* # of covariates for response model */
	       double *beta0,  /* prior mean for betaC and betaA */ 
	       double *gamma0, /* prior mean for gamma and gamma1 */
	       double *delta0, /* prior mean for delta */
	       double *dA0C,   /* prior precision for betaC and betaA */ 
	       double *dA0O,   /* prior precision for gamma and gamma1 */
	       double *dA0R,   /* prior precision for delta */
	       int *nu0,       /* prior df for sig2 */
	       double *s0,     /* prior scale for sig2 */
	       double *VarC,    /* proposal variance for compliance
				   model */
	       double *VarR,   /* proposal variance for response model */
	       int *logitC,    /* Use logistic regression for the
				  compliance model? */
	       int *logitR,    /* Use logistic regression for the
				  response model? */
	       int *param,     /* Want to keep paramters? */
	       int *mda,       /* Want to use marginal data
				  augmentation for probit regressions? */
	       int *burnin,   /* number of burnin */
	       int *iKeep,     /* keep ?th draws */
	       int *verbose,   /* print out messages */
	       double *coefC,  /* Storage for coefficients of the
				  compliance model */
	       double *coefA,  /* Storage for coefficients of the
				  always-takers model */
	       double *coefO,  /* Storage for coefficients of the
				  outcome model; gamma and gamma1 */
	       double *coefO1,  /* Storage for coefficients of the
				  outcome model; gamma and gamma1 */
	       double *coefR,  /* Storage for coefficients of the
				  response model */	      
	       double *var,    /* Storage for sig2 */
	       double *QoI     /* Storage of quantities of interest */
	       ) {
  /** counters **/
  int n_samp = in_samp[0];
  int n_samp1 = in_samp[1];
  int n_covC = *in_covC;
  int n_covO = *in_covO;
  int n_covR = *in_covR;
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */
  double **Xo = doubleMatrix(n_samp+n_covO, n_covO+1);    
  /* covariates for the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_covO, n_covO+1);    
  double **Xobs1 = doubleMatrix(n_samp1+n_covO, n_covO+1);    
  int *Yobs = intArray(n_obs); 
  /* covariates for the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_covR, n_covR+1);    
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);
  double *meano1 = doubleArray(n_samp);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);
  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  /* prior precision matrices */
  double **A0C = doubleMatrix(n_covC*2, n_covC*2);
  double **A0O = doubleMatrix(n_covO, n_covO);
  double **A0R = doubleMatrix(n_covR, n_covR);
  
  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int *acceptC = intArray(n_covC*2);      /* number of acceptance */
  int *acceptR = intArray(n_covR);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempO1, itempQ, itempR, itempS;
  double dtemp, dtemp1;
  double **mtemp = doubleMatrix(n_covO, n_covO);

  /*** get random seed **/
  GetRNGstate();

  /*** Preparing ***/
  Prep(dXc, Xc, dXo, Xo, dXr, Xr, Xobs, R, n_samp, n_obs, n_covC,
       n_covO, n_covR, *logitC, *AT, dA0C, A0C, dA0O, A0O, 1, dA0R, 
       A0R, beta0, delta0, gamma0, pC, pN, pA, prC, prN, prA, acceptC,
       acceptR, n_miss); 

  /*** observed Y ***/
  itemp = 0;
  for (i = 0; i < n_samp; i++) {
    if ((R[i] == 1) && (Y[i] == 1)) {
      for (j = 0; j < n_covO; j++)
	Xobs1[itemp][j] = Xo[i][j];
      Xobs1[itemp++][n_covO] = log(Y1[i]);
    }
  }
  dcholdc(A0O, n_covO, mtemp);
  for (i = 0; i < n_covO; i++) {
    Xobs1[itemp+i][n_covO] = 0; 
    for (j = 0; j < n_covO; j++) {
      Xobs1[itemp+i][n_covO] += mtemp[i][j]*gamma0[j];
      Xobs1[itemp+i][j] = mtemp[i][j];
    }
  }
      
  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempO1 = 0; 
  itempQ = 0; itempR = 0; itempS = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      Response(*logitR, R, Xr, delta, n_samp, n_covR, delta0, A0R, VarR,
	       acceptR, *mda, *AT, Z, D, prC, prN, prA);

    /** Step 2: COMPLIANCE MODEL **/
    Compliance(*logitC, *AT, C, Xc, betaC, n_samp, n_covC, beta0, A0C, 
	       betaA, VarC, acceptC, *mda, A);

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
    SampleComp(n_samp, n_covC, *AT, Xc, Xo, Xr, Xobs, betaC, betaA,
	       *logitC, qC, qN, Z, D, R, RD, C, A, pC, pN, pA, prA, prN, prC);

    /** Step 4: OUTCOME MODEL **/
    bprobitGibbs(Yobs, Xobs, gamma, n_obs, n_covO, 0, gamma0, A0O,
		 *mda, 1);
    bNormalReg(Xobs1, gamma1, sig2, n_samp1, n_covO, 0, 1, gamma0, A0O, 1,
	       *nu0, *s0, 0);

    /** Compute probabilities of Y = Yobs **/
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0; meano1[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_covO; j++) {
	  meano[i] += Xo[i][j]*gamma[j];
	  meano1[i] += Xo[i][j]*gamma1[j];
	}
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (Y[i] == 1) {
	      pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]], sqrt(*sig2), 0) * 
		pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0);
	      pN[i] = dlnorm(Y1[i], meano1[i], sqrt(*sig2), 0) * 
		pnorm(meano[i], 0, 1, 1, 0);
	      pA[i] = dlnorm(Y1[i], meano1[i]+gamma1[2], sqrt(*sig2), 0) * 
		pnorm(meano[i]+gamma[2], 0, 1, 1, 0) ;
	    } else {
	      pC[i] = pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pN[i] = pnorm(meano[i], 0, 1, 0, 0);
	      pA[i] = pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	    }
	  }
	}
      } else { /* no always-takers */
	for (j = 2; j < n_covO; j++) {
	  meano[i] += Xo[i][j]*gamma[j];
	  meano1[i] += Xo[i][j]*gamma1[j];
	}
	if (R[i] == 1) {
	  if ((Z[i] == 0) || (RD[i] == 0)){
	    if (Y[i] == 1) {
	      pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]], sqrt(*sig2), 0) * 
		pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0);
	      pN[i] = dlnorm(Y1[i], meano1[i], sqrt(*sig2), 0) *
		pnorm(meano[i], 0, 1, 1, 0);
	    } else {
	      pC[i] = pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pN[i] = pnorm(meano[i], 0, 1, 0, 0);
	    }
	  }
	} 
      }
    }
    
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	n_comp[0] = 0; n_comp[1] = 0; 
	n_never[0] = 0; n_never[1] = 0;
	n_always[0] = 0; n_always[1] = 0;
	p_comp = 0; p_never = 0; ITT = 0;
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i]; p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	  } else if (A[i] == 1) {
	    if (Z[i] == 1)
	      n_always[1]++;
	    else
	      n_always[0]++;
	  } else {
	    if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	  }
	  if (*Insample) { /* insample QoI */
	    if (C[i] == 1) {
	      dtemp = rlnorm(meano1[i]+gamma1[0], sqrt(*sig2)) *
		((meano[i]+gamma[0]+norm_rand()) > 0);
	      dtemp1 = rlnorm(meano1[i]+gamma1[1], sqrt(*sig2)) *
		((meano[i]+gamma[1]+norm_rand()) > 0);;
	      if (R[i] == 1) {
		if (Z[i] == 1) {
		  dtemp = Y1[i];
		} else { 
		  dtemp1 = Y1[i];
		}
	      }
	      Y1barC += dtemp; Y0barC += dtemp1;
	    } else if (A[i] == 1) {
	      if (R[i] == 1)
		YbarA += Y1[i];
	      else 
		YbarA += rlnorm(meano1[i]+gamma1[2], sqrt(*sig2)) * 
		  ((meano[i]+gamma[2]+norm_rand()) > 0);
	    } else {
	      if (R[i] == 1)
		YbarN += Y1[i];
	      else 
		YbarN += rlnorm(meano1[i], sqrt(*sig2)) *
		  ((meano[i]+gamma[2]+norm_rand()) > 0);
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    dtemp = exp(meano1[i]+gamma1[0]+0.5*sig2[0]) * 
	      pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
	    dtemp1 = exp(meano1[i]+gamma1[1]+0.5*sig2[0]) * 
	      pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	    Y1barC += dtemp; Y0barC += dtemp1;
	    ITT += ((dtemp-dtemp1)*qC[i]);
	    /* always-takers */
	    YbarA += exp(meano1[i]+gamma1[2]+0.5*sig2[0]) * 
	      pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	    /* never-takers */
	    YbarN += exp(meano1[i]+0.5*sig2[0]) * pnorm(meano[i], 0, 1, 1, 0);
	  }
	}

	p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
				      a complier */ 
	p_never /= (double)n_samp; /* Prob. of being a never-taker */
	if (*Insample) {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  YbarN /= (double)(n_never[0]+n_never[1]);
	  if (*AT)
	    YbarA /= (double)(n_always[0]+n_always[1]);
	} else {
	  ITT /= (double)n_samp;     /* ITT effect */
	  Y1barC /= (double)n_samp;
	  Y0barC /= (double)n_samp;
	  YbarN /= (double)n_samp;
	  if (*AT)
	    YbarA /= (double)n_samp;
	}
	CACE = Y1barC-Y0barC;    /* CACE */

	QoI[itempQ++] = ITT;   
	QoI[itempQ++] = CACE;   
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	QoI[itempQ++] = Y1barC;
	QoI[itempQ++] = Y0barC;
	QoI[itempQ++] = YbarN;
	if (*AT)
	  QoI[itempQ++] = YbarA;

	if (*param) {
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  var[itempS++] = sig2[0];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaC[j+n_covC];
	    } else {
	      for (j = 0; j < n_covC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	  for (j = 0; j < n_covO; j++)
	    coefO1[itempO1++] = gamma1[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_covR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_covC*2; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_covC; j++)
	      Rprintf("%10g", (double)acceptC[j]/(double)main_loop);
	}
	if (*logitR) {
	  Rprintf("\n  Current Acceptance Ratio for the response model:");
	  for (j = 0; j < n_covR; j++)
	    Rprintf("%10g", (double)acceptR[j]/(double)main_loop);
	  Rprintf("\n");
	}
	itempP += ftrunc((double) *n_gen/10); 
	progress++;
	R_FlushConsole(); 
      }
    }
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp+n_covO);
  FreeMatrix(Xobs, n_obs+n_covO);
  FreeMatrix(Xobs1, n_samp1+n_covO);
  FreeMatrix(Xr, n_samp+n_covR);
  FreeMatrix(mtemp, n_covO);
  free(Yobs);
  free(meano);
  free(meano1);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_covC*2);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(A0R, n_covR);
  free(acceptC);
  free(acceptR);
  free(n_comp);
  free(n_never);
  free(n_always);

} /* end of LItwopart */


