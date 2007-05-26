/****

     This file contains models for clustered randomized experiments
     (Frangakis, Rubin, and Zhou Biostatistics) with noncompliance
     under the assumptions of latent ignorability of Frangakis and
     Rubin (1999, Biometrika) .

****/

#include <stdio.h>      
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/*
  Preparing the data
*/

void PrepMixed(double *dXc, double *dZc, double *dXo, double *dZo,
	       double *dXr, double *dZr, double **Xc, double **Xo,
	       double **Xr, double **Xobs, double ***Zc,
	       double ***Zo, double ***Zr, int n_samp, int n_grp,
	       int n_obs, int n_miss, int n_fixedC, int n_fixedO,
	       int n_fixedR, int n_randomC, int n_randomO,
	       int n_randomR, double ***Zobs, int *R,
	       int *grp, int *grp_obs, double ***xiC, double **xiO,
	       double **xiR, double *dPsiC, double *dPsiA,
	       double *dPsiO, double *dPsiR, double ***Psi,
	       double **PsiO, double **PsiR, double *dT0C,
	       double **T0C, double *dT0A, double **T0A,
	       double *dT0O, double **T0O, double *dT0R,
	       double **T0R, double *dA0C, double **A0C,
	       double *dA0O, double **A0O, double *dA0R,
	       double **A0R, int logitC, double *pC, double *pN,
	       double *pA, double *prC, double *prN, double *prA,
	       int AT, double *beta0, double *gamma0, double *delta0,
	       int prior) {
  int i, j, k;
  int itemp = 0;
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);
  double **mtempC = doubleMatrix(n_fixedC, n_fixedC); 
  double **mtempO = doubleMatrix(n_fixedO, n_fixedO); 
  double **mtempR = doubleMatrix(n_fixedR, n_fixedR); 

  /*** pack the fixed effects covariates ***/
  itemp = 0;
  for (j = 0; j < n_fixedC; j++)
    for (i = 0; i < n_samp; i++)
      Xc[i][j] = dXc[itemp++];

  itemp = 0;
  for (j = 0; j < n_fixedO; j++)
    for (i = 0; i < n_samp; i++)
      Xo[i][j] = dXo[itemp++];

  itemp = 0;
  for (j = 0; j < n_fixedR; j++)
    for (i = 0; i < n_samp; i++)
      Xr[i][j] = dXr[itemp++];

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) {
      grp_obs[itemp] = grp[i];
      for (j = 0; j < n_fixedO; j++)
	Xobs[itemp][j] = Xo[i][j];
      itemp++;
    }

  /** pack random effects **/
  itemp = 0;
  for (k = 0; k < n_randomC; k++)
    for (j = 0; j < n_grp; j++)
      for (i = 0; i < 2; i++)
	xiC[i][j][k] = norm_rand();

  itemp = 0;
  for (k = 0; k < n_randomO; k++)
    for (j = 0; j < n_grp; j++)
      xiO[j][k] = norm_rand();

  itemp = 0;
  for (k = 0; k < n_randomR; k++)
    for (j = 0; j < n_grp; j++)
      xiR[j][k] = norm_rand();

  /** pack random effects covariates **/
  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    for (j = 0; j < n_randomC; j++)
      Zc[grp[i]][vitemp[grp[i]]][j] = dZc[itemp++];
    vitemp[grp[i]]++;
  }

  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    for (j = 0; j < n_randomO; j++)
      Zo[grp[i]][vitemp[grp[i]]][j] = dZo[itemp++];
    vitemp[grp[i]]++;
  }

  itemp = 0;
  for (j = 0; j < n_grp; j++) {
    vitemp[j] = 0; 
    vitemp1[j] = 0;
  }
  for (i = 0; i < n_samp; i++) { 
    if (R[i] == 1) {
      for (j = 0; j < n_randomO; j++)
	Zobs[grp[i]][vitemp1[grp[i]]][j] = Zo[grp[i]][vitemp[grp[i]]][j];
      vitemp1[grp[i]]++;
    }
    vitemp[grp[i]]++;
  }

  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    for (j = 0; j < n_randomR; j++)
      Zr[grp[i]][vitemp[grp[i]]][j] = dZr[itemp++];
    vitemp[grp[i]]++;
  }

  /* covariance parameters and its prior */
  itemp = 0;
  for (k = 0; k < n_randomC; k++)
    for (j = 0; j < n_randomC; j++) 
      Psi[0][j][k] = dPsiC[itemp++];
 
  itemp = 0;
  for (k = 0; k < n_randomC; k++)
    for (j = 0; j < n_randomC; j++)
      T0C[j][k] = dT0C[itemp++];

  itemp = 0;
  for (k = 0; k < n_randomC; k++)
    for (j = 0; j < n_randomC; j++) 
      Psi[1][j][k] = dPsiA[itemp++];

  itemp = 0;
  for (k = 0; k < n_randomC; k++)
    for (j = 0; j < n_randomC; j++)
      T0A[j][k] = dT0A[itemp++];

  itemp = 0;
  for (k = 0; k < n_randomO; k++)
    for (j = 0; j < n_randomO; j++)
      PsiO[j][k] = dPsiO[itemp++];
  itemp = 0;
  for (k = 0; k < n_randomO; k++)
    for (j = 0; j < n_randomO; j++)
      T0O[j][k] = dT0O[itemp++];

  itemp = 0;
  for (k = 0; k < n_randomR; k++)
    for (j = 0; j < n_randomR; j++)
      PsiR[j][k] = dPsiR[itemp++];
  itemp = 0;
  for (k = 0; k < n_randomR; k++)
    for (j = 0; j < n_randomR; j++)
      T0R[j][k] = dT0R[itemp++];

  itemp = 0; 
  if ((logitC == 1) && (AT == 1))
    for (k = 0; k < n_fixedC*2; k++)
      for (j = 0; j < n_fixedC*2; j++)
	A0C[j][k] = dA0C[itemp++];
  else
    for (k = 0; k < n_fixedC; k++)
      for (j = 0; j < n_fixedC; j++)
	A0C[j][k] = dA0C[itemp++];
  
  itemp = 0;
  for (k = 0; k < n_fixedO; k++)
    for (j = 0; j < n_fixedO; j++)
      A0O[j][k] = dA0O[itemp++];

  itemp = 0;
  for (k = 0; k < n_fixedR; k++)
    for (j = 0; j < n_fixedR; j++)
      A0R[j][k] = dA0R[itemp++];

  /* prior for fixed effects as additional data points */ 
  if (logitC != 1) {
    dcholdc(A0C, n_fixedC, mtempC);
    for (i = 0; i < n_fixedC; i++) {
      Xc[n_samp+i][n_fixedC] = 0;
      for (j = 0; j < n_fixedC; j++) {
	Xc[n_samp+i][n_fixedC] += mtempC[i][j]*beta0[j];
	Xc[n_samp+i][j] = mtempC[i][j];
      }
    }
  }

  dcholdc(A0R, n_fixedR, mtempR);
  for (i = 0; i < n_fixedR; i++) {
    Xr[n_samp+i][n_fixedR] = 0;
    for (j = 0; j < n_fixedR; j++) {
      Xr[n_samp+i][n_fixedR] += mtempR[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtempR[i][j];
    }
  }
  
  if (prior) {
    dcholdc(A0O, n_fixedO, mtempO);
    for (i = 0; i < n_fixedO; i++) {
      Xobs[n_obs+i][n_fixedO] = 0;
      for (j = 0; j < n_fixedO; j++) {
	Xobs[n_obs+i][n_fixedO] += mtempO[i][j]*gamma0[j];
	Xobs[n_obs+i][j] = mtempO[i][j];
      }
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

  free(vitemp);
  free(vitemp1);
  FreeMatrix(mtempC, n_fixedC);
  FreeMatrix(mtempO, n_fixedO);
  FreeMatrix(mtempR, n_fixedR);
}

/* 
   Response probabilities 
*/

void ResponseMixed(int *R, double **Xr, double ***Zr, 
		   int *grp, double *delta, double **xiR, 
		   double **PsiR, int n_samp, int n_fixedR, 
		   int n_randomR, int n_grp, double *delta0, 
		   double **A0R, int *tau0s, double **T0R,
		   int AT, int random, int *Z, int *D, double *prC,
		   double *prN, double *prA
		   ){
  int i, j;
  double dtemp;
  int *vitemp = intArray(n_grp);

  bprobitMixedGibbs(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, 0, 
		    delta0, A0R, tau0s[3], T0R, 1);
  
  /* Compute probabilities of R = Robs */ 
  for (j = 0; j < n_grp; j++) vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    dtemp = 0;
    if (AT) { /* always-takers */
      for (j = 3; j < n_fixedR; j++)
	dtemp += Xr[i][j]*delta[j];
      if (random) 
	for (j = 2; j < n_randomR; j++)
	  dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
      else
	for (j = 0; j < n_randomR; j++)
	  dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
      if (random) {
	if (Z[i] == 0)  
	  prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	else
	  prC[i] = R[i]*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 0, 0);
	prA[i] = R[i]*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 1, 0) +
	  (1-R[i])*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 0, 0);
      } else {
	if (Z[i] == 0)
	  prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	else
	  prC[i] = R[i]*pnorm(dtemp+delta[0], 0, 1, 1, 0) +
	    (1-R[i])*pnorm(dtemp+delta[0], 0, 1, 0, 0);
	prA[i] = R[i]*pnorm(dtemp+delta[2], 0, 1, 1, 0) +
	  (1-R[i])*pnorm(dtemp+delta[2], 0, 1, 0, 0);
      }
      prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	(1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
    } else { /* no always-takers */
      for (j = 2; j < n_fixedR; j++)
	dtemp += Xr[i][j]*delta[j];
      if (random) 
	for (j = 1; j < n_randomR; j++)
	  dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
      else
	for (j = 0; j < n_randomR; j++)
	  dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
      if (random) {
	if (Z[i] == 0) 
	  prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) + 
	    (1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	else
	  prC[i] = R[i]*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 1, 0) + 
	    (1-R[i])*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 0, 0);
      } else {
	if (Z[i] == 0)
	  prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) + 
	    (1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	else
	  prC[i] = R[i]*pnorm(dtemp+delta[0], 0, 1, 1, 0) + 
	    (1-R[i])*pnorm(dtemp+delta[0], 0, 1, 0, 0);
      }
      prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	(1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
    } 
    vitemp[grp[i]]++;
  }
  
  free(vitemp);
}


/* 
   Compliance model
*/

void CompMixed(int logitC, int AT, int *C, double **Xc, double ***Zc,
	       int *grp, double *betaC, double ***xiC, double ***Psi,
	       int n_samp, int n_fixedC, int n_randomC , int n_grp,
	       double *beta0, double **A0C, int *tau0s, double **T0C, 
	       double *tune_fixed, double *tune_random, int *acc_fixed,
	       int *acc_random, int *A, int max_samp_grp, 
	       double *betaA, double **T0A
	       ){
  int i, j;
  int itemp;
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);

  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  int *Atemp = intArray(n_samp);
  double ***Ztemp = doubleMatrix3D(n_grp, max_samp_grp + n_randomC,
				   n_randomC +1);
  int *grp_temp = intArray(n_samp);
  
  if (logitC) 
    if (AT) 
      logitMixedMetro(C, Xc, Zc, grp, betaC, xiC, Psi, n_samp, 2,
		      n_fixedC, n_randomC, n_grp, beta0, A0C, tau0s[0],
		      T0C, tune_fixed, tune_random, 1, acc_fixed, acc_random);
    else 
      logitMixedMetro(C, Xc, Zc, grp, betaC, xiC, Psi, n_samp, 1,
		      n_fixedC, n_randomC, n_grp, beta0, A0C,
		      tau0s[0], T0C, tune_fixed, tune_random, 1,
		      acc_fixed, acc_random);
  else {
    /* complier vs. noncomplier */
    bprobitMixedGibbs(C, Xc, Zc, grp, betaC, xiC[0], Psi[0], n_samp,
		      n_fixedC, n_randomC, n_grp, 0, beta0, A0C,
		      tau0s[0], T0C, 1); 
    if (AT) {
      /* never-taker vs. always-taker */
      /* subset the data */
      itemp = 0;
      for (j = 0; j < n_grp; j++) {
	vitemp[j] = 0; vitemp1[j] = 0;
      }
      for (i = 0; i < n_samp; i++) {
	if (C[i] == 0) {
	  Atemp[itemp] = A[i]; grp_temp[itemp] = grp[i];
	  for (j = 0; j < n_fixedC; j++)
	    Xtemp[itemp][j] = Xc[i][j];
	  for (j = 0; j < n_randomC; j++)
	    Ztemp[grp[i]][vitemp1[grp[i]]][j] = Zc[grp[i]][vitemp[grp[i]]][j];
	  itemp++; vitemp1[grp[i]]++;
	}
	vitemp[grp[i]]++;
      }
      for (i = n_samp; i < n_samp + n_fixedC; i++) {
	for (j = 0; j <= n_fixedC; j++)
	  Xtemp[itemp][j] = Xc[i][j];
	itemp++;
      }
      bprobitMixedGibbs(Atemp, Xtemp, Ztemp, grp_temp, betaA, xiC[1],
			Psi[1], itemp-n_fixedC, n_fixedC, n_randomC,
			n_grp, 0, beta0, A0C, tau0s[1], T0A, 1); 
    }      
  }    

  free(vitemp);
  free(vitemp1);
  FreeMatrix(Xtemp, n_samp+n_fixedC);
  free(Atemp);
  Free3DMatrix(Ztemp, n_grp, max_samp_grp + n_randomC);
  free(grp_temp);
}

/*
  Sample Compliance Covariate
*/

void SampCompMixed(int n_grp, int n_samp, int n_fixedC, double **Xc, 
		   double *betaC, double ***Zc, int *grp, 
		   double ***xiC, int n_randomC, int AT, int logitC,
		   double *qC, double *qN, int *Z, int *D, int *R, 
		   int *RD, double *prC, double *prN, double ***Zo, 
		   double ***Zr, int *C, double **Xo, double **Xr,
		   int random, double **Xobs, double ***Zobs, 
		   double *prA, double *pA, int *A, double *betaA, 
		   double *pC, double *pN
		   ) {

  int i, j;
  int itemp;
  double dtemp, dtemp1, dtemp2;
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);

  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);
  double *meana = doubleArray(n_samp);

  itemp = 0;
  for (j = 0; j < n_grp; j++) {
    vitemp[j] = 0; vitemp1[j] = 0;
  }
  for (i = 0; i < n_samp; i++) {
    meanc[i] = 0;
    for (j = 0; j < n_fixedC; j++) 
      meanc[i] += Xc[i][j]*betaC[j];
    for (j = 0; j < n_randomC; j++)
      meanc[i] += Zc[grp[i]][vitemp[grp[i]]][j]*xiC[0][grp[i]][j];
    if (AT) { /* some always-takers */
      meana[i] = 0;
      for (j = 0; j < n_randomC; j++)
	meana[i] += Zc[grp[i]][vitemp[grp[i]]][j]*xiC[1][grp[i]][j];
      if (logitC) { /* if logistic regression is used */
	for (j = 0; j < n_fixedC; j++) 
	  meana[i] += Xc[i][j]*betaC[j+n_fixedC];
	qC[i] = exp(meanc[i])/(1 + exp(meanc[i]) + exp(meana[i]));
	qN[i] = 1/(1 + exp(meanc[i]) + exp(meana[i]));
      } else { /* double probit regressions */
	for (j = 0; j < n_fixedC; j++) 
	  meana[i] += Xc[i][j]*betaA[j];
	qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	qN[i] = (1-qC[i])*pnorm(meana[i], 0, 1, 0, 0);
      }
      if (RD[i] == 0) { /* units with missing treatment status */
	if (R[i] == 1) {
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	  dtemp1 = (qC[i]*pC[i]*prC[i] + qN[i]*pN[i]*prN[i]) /
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	} else {
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+qN[i]*prN[i]+(1-qC[i]-qN[i])*prA[i]);
	  dtemp1 = (qC[i]*prC[i] + qN[i]*prN[i]) / 
		   (qC[i]*prC[i]+qN[i]*prN[i]+(1-qC[i]-qN[i])*prA[i]);
	}
	dtemp2 = unif_rand();
	if (dtemp2 < dtemp) { /* compliers */
	  C[i] = 1; A[i] = 0; D[i] = Z[i];
	  Xo[i][1-Z[i]] = 1; Xr[i][1-Z[i]] = 1;
	  Xo[i][Z[i]] = 0; Xr[i][Z[i]] = 0;
	  Xo[i][2] = 0; Xr[i][2] = 0;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	  }
	} else if (dtemp2 < dtemp1) { /* never-takers */
	  C[i] = 0; A[i] = 0; D[i] = 0;
	  Xo[i][0] = 0; Xr[i][0] = 0;
	  Xo[i][1] = 0; Xr[i][1] = 0;
	  Xo[i][2] = 0; Xr[i][2] = 0;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	  }
	} else { /* always-takers */
	  if (logitC)
	    C[i] = 2;
	  else
	    C[i] = 0; 
	  A[i] = 1; D[i] = 1;
	  Xo[i][0] = 0; Xr[i][0] = 0; 
	  Xo[i][1] = 0; Xr[i][1] = 0; 
	  Xo[i][2] = 1; Xr[i][2] = 1;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zo[grp[i]][vitemp[grp[i]]][1] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 1;
	  }
	} 	
      } else if ((Z[i] == 0) && (D[i] == 0)) {
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i]/(qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]);
	else
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+qN[i]*prN[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	  }
	} else {
	  C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	  }
	}  
      } else if ((Z[i] == 1) && (D[i] == 1)){
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i]/(qC[i]*pC[i]*prC[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	else
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i]-qN[i])*prA[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; Xo[i][0] = 1; Xr[i][0] = 1; 
	  A[i] = 0; Xo[i][2] = 0; Xr[i][2] = 0; 
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	  } 
	} else {
	  if (logitC)
	    C[i] = 2;
	  else
	    C[i] = 0; 
	  A[i] = 1; Xo[i][0] = 0; Xr[i][0] = 0; 
	  Xo[i][2] = 1; Xr[i][2] = 1;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zo[grp[i]][vitemp[grp[i]]][1] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 1;
	  }
	}
      }  
      if (R[i] == 1) {
	Xobs[itemp][0] = Xo[i][0];
	Xobs[itemp][1] = Xo[i][1];
	Xobs[itemp][2] = Xo[i][2];
	if (random) {
	  Zobs[grp[i]][vitemp1[grp[i]]][0] = Zo[grp[i]][vitemp[grp[i]]][0]; 
	  Zobs[grp[i]][vitemp1[grp[i]]][1] = Zo[grp[i]][vitemp[grp[i]]][1]; 
	}
      }
    } else { /* no always-takers */
      if ((Z[i] == 0) || (RD[i] == 0)){
	if (logitC)
	  qC[i] = 1/(1+exp(-meanc[i]));
	else
	  qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	if (R[i] == 1)
	  dtemp = qC[i]*pC[i]*prC[i]/(qC[i]*pC[i]*prC[i]+(1-qC[i])*pN[i]*prN[i]);
	else
	  dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i])*prN[i]);
	if (unif_rand() < dtemp) {
	  C[i] = 1; D[i] = Z[i];
	  Xo[i][1-Z[i]] = 1; Xo[i][Z[i]] = 0; 
	  Xr[i][1-Z[i]] = 1; Xr[i][Z[i]] = 0;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	  }
	} else {
	  C[i] = 0;  D[i] = 0;
	  Xo[i][0] = 0; Xr[i][0] = 0;
	  Xo[i][1] = 0; Xr[i][1] = 0;
	  if (random) {
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	  }
	}
      }
      if (R[i] == 1) {
	Xobs[itemp][0] = Xo[i][0];
	Xobs[itemp][1] = Xo[i][1];
	if (random) {
	  Zobs[grp[i]][vitemp1[grp[i]]][0] = Zo[grp[i]][vitemp[grp[i]]][0];
	  Zobs[grp[i]][vitemp1[grp[i]]][1] = Zo[grp[i]][vitemp[grp[i]]][1];
	}
      } 
    }
    if (R[i] == 1) {
      itemp++; vitemp1[grp[i]]++;
    }
    vitemp[grp[i]]++;
  }
  
  free(vitemp);
  free(vitemp1);
  free(meana);
  free(meanc);
}

/* Calculating univariate QoI */
void uniQoIcalMixed(int Insample, int n_grp, double *ITT, double *Y1barC,
		    double *Y0barC, int n_samp, int **n_comp, int **n_never, int **n_always,
		    double *p_comp, double *p_never, double *CACE, double *YbarN,
		    double *YbarA, int AT) {
  int j;
  
  for (j = 0; j < (n_grp+1); j++) {	  
    p_comp[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			  n_comp[j][1]+n_never[j][1]+n_always[j][1]);  
    p_never[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			   n_comp[j][1]+n_never[j][1]+n_always[j][1]);  
    if (Insample) { /* insample QoI */
      ITT[j] = (Y1barC[j]-Y0barC[j]) /
	(double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
		 n_comp[j][1]+n_never[j][1]+n_always[j][1]);
      Y1barC[j] /= (double)(n_comp[j][0]+n_comp[j][1]);  
      Y0barC[j] /= (double)(n_comp[j][0]+n_comp[j][1]); 
      YbarN[j] /= (double)(n_never[j][0]+n_never[j][1]);
      if (AT)
	YbarA[j] /= (double)(n_always[j][0]+n_always[j][1]);
    } else { /* population QoI */
      ITT[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			 n_comp[j][1]+n_never[j][1]+n_always[j][1]);
      Y1barC[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			    n_comp[j][1]+n_never[j][1]+n_always[j][1]);
      Y0barC[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			    n_comp[j][1]+n_never[j][1]+n_always[j][1]);
      YbarN[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			   n_comp[j][1]+n_never[j][1]+n_always[j][1]);
      if (AT)
	YbarA[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			     n_comp[j][1]+n_never[j][1]+n_always[j][1]);
    }
    CACE[j] = Y1barC[j]-Y0barC[j]; 
  }
}

/* 
   Bayesian binary probit mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LIbprobitMixed(int *Y,         /* binary outcome variable */ 
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
		    int *grp,       /* group indicator: 0, 1, 2, ... */
		    int *Ymiss,     /* number of missing obs in Y */
		    int *AT,        /* Are there always-takers? */
		    int *Insample,  /* Insample (=1) or population QoI? */
		    int *random,    /* Want to include compliance
				       random effects? */
		    double *dXc,    /* fixed effects for compliance model */
		    double *dZc,    /* random effects for compliance model */
		    double *dXo,    /* fixed effects for outcome model */
		    double *dZo,    /* random effects for outcome model */
		    double *dXr,    /* fixed effects for response model */
		    double *dZr,    /* random effects for response model */
		    double *betaC,   /* fixed effects for compliance model */
		    double *betaA,   /* fixed effects for always-takers model */
		    double *gamma,  /* fixed effects for outcome model */
		    double *delta,  /* fixed effects for response model */
		    int *in_samp,   /* # of observations */
		    int *n_gen,     /* # of Gibbs draws */
		    int *in_grp,     /* # of groups */
		    int *max_samp_grp, /* max # of obs within group */
		    int *in_fixed,  /* # of fixed effects for
				       compliance, outcome, and response models */
		    int *in_random, /* # of random effects for
				       compliance, outcome, and response modeld */ 
		    double *dPsiC,  /* precision for compliance model */
		    double *dPsiA,  /* precision for always-takers model */
		    double *dPsiO,  /* precision for outcome model */
		    double *dPsiR,  /* precision for response model */
		    double *beta0,  /* prior mean for betaC and betaA */ 
		    double *gamma0, /* prior mean for gamma */
		    double *delta0, /* prior mean for delta */
		    double *dA0C,   /* prior precision for betaC and betaA */ 
		    double *dA0O,   /* prior precision for gamma */
		    double *dA0R,   /* prior precision for delta */
		    int *tau0s,     /* prior df for PsiC, PsiA, PsiO, PsiR */
		    double *dT0C,   /* prior scale for PsiC */
		    double *dT0A,   /* prior scale for PsiA */
		    double *dT0O,   /* prior scale for PsiO */
		    double *dT0R,   /* prior scale for PsiR */
		    double *tune_fixed, /* proposal variance */
		    double *tune_random, /* proposal variance */
		    int *logitC,    /* Use logistic regression for the
				       compliance model? */
		    int *param,     /* Want to keep paramters? */
		    int *burnin,    /* number of burnin */
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
		    double *sPsiC,  /* Storage for precisions of
				       random effects */
		    double *sPsiA,
		    double *sPsiO,
		    double *sPsiR, 
		    double *QoI     /* Storage of quantities of
				       interest */		 
		    ) {
   /** counters **/
  int n_samp = *in_samp; 
  int n_grp = *in_grp;
  int n_fixedC = in_fixed[0]; int n_randomC = in_random[0];
  int n_fixedO = in_fixed[1]; int n_randomO = in_random[1];
  int n_fixedR = in_fixed[2]; int n_randomR = in_random[2];
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  int *grp_obs = intArray(n_obs);

  /*** observed Y ***/
  int *Yobs = intArray(n_obs);

  /* covariates for fixed effects in the compliance model */
  double **Xc = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  /* covariates for random effects */
  double ***Zc = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				n_randomC +1);

  /* covariates for fixed effects in the outcome model */
  double **Xo = doubleMatrix(n_samp+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				n_randomO +1);

  /* covariates for fixed effecs in the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zobs = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				  n_randomO +1);

  /* covariates for fixed effects in the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_fixedR, n_fixedR+1);    
  /* covariates for random effects */
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_randomR,
				n_randomR +1);

  /*** model parameters ***/
  /* random effects */
  double ***xiC = doubleMatrix3D(2, n_grp, n_randomC);
  double **xiO = doubleMatrix(n_grp, n_randomO);
  double **xiR = doubleMatrix(n_grp, n_randomR);

  /* covariances for random effects */
  double ***Psi = doubleMatrix3D(2, n_randomC, n_randomC);
  double **PsiO = doubleMatrix(n_randomO, n_randomO);
  double **PsiR = doubleMatrix(n_randomR, n_randomR);

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

  /* probability of being Y = 1 for always-taker */
  double *pA = doubleArray(n_samp);

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* qoi by groups */
  double *ITT = doubleArray(n_grp+1);
  double *CACE = doubleArray(n_grp+1);
  double *Y1barC = doubleArray(n_grp+1);
  double *Y0barC = doubleArray(n_grp+1);
  double *YbarN = doubleArray(n_grp+1); 
  double *YbarA = doubleArray(n_grp+1);
  int **n_comp = intMatrix(n_grp+1, 2);          /* number of compliers */
  int **n_never = intMatrix(n_grp+1, 2);
  int **n_always = intMatrix(n_grp+1, 2);
  double *p_comp = doubleArray(n_grp+1); 
  double *p_never = doubleArray(n_grp+1); /* prob. of being a particular type */


  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  double dtemp, dtemp1;
  int progress = 1;
  int keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv, itempPO, itempPC, itempPA, itempPR;

  /*** get random seed **/
  GetRNGstate();

  /*** Data prep etc. ***/
  PrepMixed(dXc, dZc, dXo, dZo, dXr, dZr, Xc, Xo, Xr, Xobs, Zc, Zo,
	    Zr, n_samp, n_grp, n_obs, n_miss, n_fixedC, n_fixedO,
	    n_fixedR, n_randomC, n_randomO, n_randomR, Zobs, R,
	    grp, grp_obs, xiC, xiO, xiR, dPsiC, dPsiA, dPsiO, dPsiR,
	    Psi, PsiO, PsiR, dT0C, T0C, dT0A, T0A, dT0O, T0O, dT0R,
	    T0R, dA0C, A0C, dA0O, A0O, dA0R, A0R, *logitC, pC, pN,
	    pA, prC, prN, prA, *AT, beta0, gamma0, delta0, 1);

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0;   
  itempPO = 0; itempPA = 0; itempPC = 0; itempPR = 0;
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      ResponseMixed(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, delta0, A0R, tau0s, T0R,
		    *AT, *random, Z, D, prC,prN, prA);
    
    /** Step 2: COMPLIANCE MODEL **/
    CompMixed(*logitC, *AT, C, Xc, Zc, grp, betaC, xiC, Psi, n_samp,
	      n_fixedC, n_randomC , n_grp, beta0, A0C, tau0s, T0C, 
	      tune_fixed, tune_random, acc_fixed, acc_random, A, 
	      *max_samp_grp, betaA, T0A);

    /** Step 3: SAMPLE COMPLIANCE COVARIATE **/
    SampCompMixed(n_grp, n_samp, n_fixedC, Xc, betaC, Zc, grp, 
		  xiC, n_randomC, *AT, *logitC, qC, qN, Z, D, R, RD, 
		  prC, prN, Zo, Zr, C, Xo, Xr, *random, Xobs, Zobs, 
		  prA, pA, A, betaA, pC, pN);

    /** Step 4: OUTCOME MODEL **/
    bprobitMixedGibbs(Yobs, Xobs, Zobs, grp_obs, gamma, xiO, PsiO,
		      n_obs, n_fixedO, n_randomO, n_grp, 0,
		      gamma0, A0O, tau0s[2], T0O, 1); 
    
    /** Compute probabilities of Y = Yobs **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 2; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (*random) {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 0, 0);
	      pA[i] = Y[i]*pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 0, 0);
	    } else {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pA[i] = Y[i]*pnorm(meano[i]+gamma[2], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	    }
	    pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	      (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	  }
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 1; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else 
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (*random)
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 1, 0) + 
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 0, 0);
	    else 
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0) + 
		(1-Y[i])*pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	    pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	      (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	  }
      }
      vitemp[grp[i]]++;
    }
    
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	for (j = 0; j < (n_grp+1); j++) {
	  n_comp[j][0] = 0; n_comp[j][1] = 0;
	  n_never[j][0] = 0; n_never[j][1] = 0;
	  n_always[j][0] = 0; n_always[j][1] = 0;
	  p_comp[j] = 0; p_never[j] = 0; ITT[j] = 0;
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp[grp[i]] += qC[i]; p_comp[n_grp] += qC[i];
	  p_never[grp[i]] += qN[i]; p_never[n_grp] += qN[i];
	  /* counting */
	  if (C[i] == 1) { 
	    if (Z[i] == 1) {
	      n_comp[grp[i]][1]++; n_comp[n_grp][1]++;
	    } else {
	      n_comp[grp[i]][0]++; n_comp[n_grp][0]++;
	    }
	  } else if (A[i] == 1) {
	    if (Z[i] == 1) {
	      n_always[grp[i]][1]++; n_always[n_grp][1]++;
	    } else {
	      n_always[grp[i]][0]++; n_always[n_grp][0]++;
	    }
	  } else {
	    if (Z[i] == 1) {
	      n_never[grp[i]][1]++; n_never[n_grp][1]++;
	    } else {
	      n_never[grp[i]][0]++; n_never[n_grp][0]++;
	    }
	  }
	  /* insample QoI */
	  if (*Insample) { 
	    if (C[i] == 1) { /* compliers */
	      if (*random) {
		dtemp = (double)((meano[i]+gamma[0]+xiO[grp[i]][0]+norm_rand()) > 0);
		dtemp1 = (double)((meano[i]+gamma[1]+xiO[grp[i]][0]+norm_rand()) > 0);
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
	      Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	      Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    } else if (A[i] == 1) { /* always-takers */ 
	      if (R[i] == 1)
		dtemp = (double)Y[i];
	      else if (*random) 
		dtemp = (double)((meano[i]+gamma[2]+xiO[grp[i]][1]+norm_rand()) > 0);
	      else
		dtemp = (double)((meano[i]+gamma[2]+norm_rand()) > 0);
	      YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    } else { /* never-takers */
	      if (R[i] == 1)
		dtemp = (double)Y[i];
	      else
		dtemp = (double)((meano[i]+norm_rand()) > 0);
	      YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    if (*random) {
	      dtemp = pnorm(meano[i]+gamma[0]+xiO[grp[i]][0], 0, 1, 1, 0);
	      dtemp1 = pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 1, 0); 
	    } else {
	      dtemp = pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
	      dtemp1 = pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	    }
	    Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	    Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    ITT[grp[i]] += ((dtemp-dtemp1)*qC[i]);
	    ITT[n_grp] += ((dtemp-dtemp1)*qC[i]);
	    /* always-takers */
	    if (*random)
	      dtemp = pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0);
	    else
	      dtemp = pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	    YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    /* never-takers */
	    dtemp = pnorm(meano[i], 0, 1, 1, 0);
	    YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	  }
	}
	
	uniQoIcalMixed(*Insample, n_grp, ITT, Y1barC, Y0barC, n_samp, n_comp,
		       n_never, n_always, p_comp, p_never, CACE, YbarN,
		       YbarA, *AT);

	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_comp[j]; 	  
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_never[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = YbarN[j];
	if (*AT)
	  for (j = 0; j < (n_grp+1); j++)
	    QoI[itempQ++] = YbarA[j];
	
	if (*param) {
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  for (i = 0; i < n_randomC; i++)
	    for (j = i; j < n_randomC; j++)
	      sPsiC[itempPC++] = Psi[0][i][j];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    } else {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	    }
	    for (i = 0; i < n_randomC; i++)
	      for (j = i; j < n_randomC; j++)
		sPsiA[itempPA++] = Psi[1][i][j];
	  }
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  for (i = 0; i < n_randomO; i++)
	    for (j = i; j < n_randomO; j++)
	      sPsiO[itempPO++] = PsiO[i][j];
	  if (n_miss > 0) {
	    for (j = 0; j < n_fixedR; j++)
	      coefR[itempR++] = delta[j];
	    for (i = 0; i < n_randomR; i++)
	      for (j = i; j < n_randomR; j++)
		sPsiR[itempPR++] = PsiR[i][j];
	  }
	}
	keep = 1;
      } else
	keep++;
    }

    if (*verbose) {
      if (main_loop == itempP) {
	Rprintf("%3d percent done.\n", progress*10);
       	if (*logitC) {
	  Rprintf("  Current Acceptance Ratio for fixed effects:");
	  if (*AT)
	    for (j = 0; j < n_fixedC*2; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_fixedC; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  Rprintf("\n");
	  Rprintf("  Current Acceptance Ratio for random effects:");
	  if (*AT)
	    for (j = 0; j < 2; j++)
	      Rprintf("%10g", (double)acc_random[j]/(double)main_loop);
	  else
	    Rprintf("%10g", (double)acc_random[0]/(double)main_loop);
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
  free(grp_obs);
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_fixedC);
  Free3DMatrix(Zc, n_grp, *max_samp_grp + n_randomC);
  FreeMatrix(Xo, n_samp+n_fixedO);
  Free3DMatrix(Zo, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xobs, n_obs+n_fixedO);
  Free3DMatrix(Zobs, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xr, n_samp+n_fixedR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_randomR);
  Free3DMatrix(xiC, 2, n_grp);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiR, n_grp);
  Free3DMatrix(Psi, 2, n_randomC);
  FreeMatrix(PsiO, n_randomO);
  FreeMatrix(PsiR, n_randomR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  free(ITT);
  free(CACE);
  free(Y1barC);
  free(Y0barC);
  free(YbarN);
  free(YbarA);
  FreeintMatrix(n_comp, n_grp+1);
  FreeintMatrix(n_never, n_grp+1);
  FreeintMatrix(n_always, n_grp+1);
  free(p_comp);
  free(p_never);
  free(vitemp);
  free(acc_fixed);
  free(acc_random);

} /* end of LIbprobitMixed */



/* 
   Bayesian Normal mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LINormalMixed(double *Y,      /* Gaussian outcome variable */ 
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
		   int *grp,       /* group indicator: 0, 1, 2, ... */
		   int *Ymiss,     /* number of missing obs in Y */
		   int *AT,        /* Are there always-takers? */
		   int *Insample,  /* Insample (=1) or population QoI? */
		   int *random,    /* Want to include compliance
				      random effects? */
		   double *dXc,    /* fixed effects for compliance model */
		   double *dZc,    /* random effects for compliance model */
		   double *dXo,    /* fixed effects for outcome model */
		   double *dZo,    /* random effects for outcome model */
		   double *dXr,    /* fixed effects for response model */
		   double *dZr,    /* random effects for response model */
		   double *betaC,  /* fixed effects for compliance model */
		   double *betaA,  /* fixed effects for always-takers model */
		   double *gamma,  /* fixed effects for outcome model */
		   double *delta,  /* fixed effects for response model */
		   double *sig2,   /* variance parameter for outcome
				      model */
		   int *in_samp,   /* # of observations */
		   int *n_gen,     /* # of Gibbs draws */
		   int *in_grp,    /* # of groups */
		   int *max_samp_grp, /* max # of obs within group */
		   int *in_fixed,  /* # of fixed effects for
				      compliance, outcome, and response models */
		   int *in_random, /* # of random effects for
				      compliance, outcome, and response modeld */ 
		   double *dPsiC,  /* precision for compliance model */
		   double *dPsiA,  /* precision for always-takers model */
		   double *dPsiO,  /* precision for outcome model */
		   double *dPsiR,  /* precision for response model */
		   double *beta0,  /* prior mean for betaC and betaA */ 
		   double *gamma0, /* prior mean for gamma */
		   double *delta0, /* prior mean for delta */
		   double *dA0C,   /* prior precision for betaC and betaA */ 
		   double *dA0O,   /* prior precision for gamma */
		   double *dA0R,   /* prior precision for delta */
		   int *nu0,       /* prior df for sig2 */
		   double *s0,     /* prior scale for sig2 */
		   int *tau0s,     /* prior df for PsiC, PsiA, PsiO, PsiR */
		   double *dT0C,   /* prior scale for PsiC */
		   double *dT0A,   /* prior scale for PsiA */
		   double *dT0O,   /* prior scale for PsiO */
		   double *dT0R,   /* prior scale for PsiR */
		   double *tune_fixed, /* proposal variance */
		   double *tune_random, /* proposal variance */
		   int *logitC,    /* Use logistic regression for the
				      compliance model? */
		   int *param,     /* Want to keep paramters? */
		   int *burnin,    /* number of burnin */
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
		   double *ssig2,  /* Storage for outcome variance */
		   double *sPsiC,  /* Storage for precisions of
				      random effects */
		   double *sPsiA,
		   double *sPsiO,
		   double *sPsiR, 
		   double *QoI     /* Storage of quantities of
				      interest */		 
		   ) {
   /** counters **/
  int n_samp = *in_samp; 
  int n_grp = *in_grp;
  int n_fixedC = in_fixed[0]; int n_randomC = in_random[0];
  int n_fixedO = in_fixed[1]; int n_randomO = in_random[1];
  int n_fixedR = in_fixed[2]; int n_randomR = in_random[2];
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  int *grp_obs = intArray(n_obs);

  /*** observed Y ***/
  double *Yobs = doubleArray(n_obs);

  /* covariates for fixed effects in the compliance model */
  double **Xc = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  /* covariates for random effects */
  double ***Zc = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				n_randomC +1);

  /* covariates for fixed effects in the outcome model */
  double **Xo = doubleMatrix(n_samp+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				n_randomO +1);

  /* covariates for fixed effecs in the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zobs = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				  n_randomO +1);

  /* covariates for fixed effects in the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_fixedR, n_fixedR+1);    
  /* covariates for random effects */
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_randomR,
				n_randomR +1);

  /*** model parameters ***/
  /* random effects */
  double ***xiC = doubleMatrix3D(2, n_grp, n_randomC);
  double **xiO = doubleMatrix(n_grp, n_randomO);
  double **xiR = doubleMatrix(n_grp, n_randomR);

  /* covariances for random effects */
  double ***Psi = doubleMatrix3D(2, n_randomC, n_randomC);
  double **PsiO = doubleMatrix(n_randomO, n_randomO);
  double **PsiR = doubleMatrix(n_randomR, n_randomR);

  /* density of Y = Yobs for complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);

  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);

  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);

  /* density of Y = Yobs for always-taker */
  double *pA = doubleArray(n_samp);

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* qoi by groups */
  double *ITT = doubleArray(n_grp+1);
  double *CACE = doubleArray(n_grp+1);
  double *Y1barC = doubleArray(n_grp+1);
  double *Y0barC = doubleArray(n_grp+1);
  double *YbarN = doubleArray(n_grp+1); 
  double *YbarA = doubleArray(n_grp+1);
  int **n_comp = intMatrix(n_grp+1, 2);          /* number of compliers */
  int **n_never = intMatrix(n_grp+1, 2);
  int **n_always = intMatrix(n_grp+1, 2);
  double *p_comp = doubleArray(n_grp+1); 
  double *p_never = doubleArray(n_grp+1); /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int progress; progress = 1;
  int keep; keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv, itempS;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Data prep etc. ***/
  PrepMixed(dXc, dZc, dXo, dZo, dXr, dZr, Xc, Xo, Xr, Xobs, Zc, Zo,
	    Zr, n_samp, n_grp, n_obs, n_miss, n_fixedC, n_fixedO,
	    n_fixedR, n_randomC, n_randomO, n_randomR, Zobs, R,
	    grp, grp_obs, xiC, xiO, xiR, dPsiC, dPsiA, dPsiO, dPsiR,
	    Psi, PsiO, PsiR, dT0C, T0C, dT0A, T0A, dT0O, T0O, dT0R,
	    T0R, dA0C, A0C, dA0O, A0O, dA0R, A0R, *logitC, pC, pN,
	    pA, prC, prN, prA, *AT, beta0, gamma0, delta0, 1);

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) {
      Yobs[itemp] = Y[i];
      Xobs[itemp++][n_fixedO] = Y[i];
    }

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0; itempS = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      ResponseMixed(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, delta0, A0R, tau0s, T0R,
		    *AT, *random, Z, D, prC,prN, prA);
    
    /** Step 2: COMPLIANCE MODEL **/
    CompMixed(*logitC, *AT, C, Xc, Zc, grp, betaC, xiC, Psi, n_samp,
	      n_fixedC, n_randomC , n_grp, beta0, A0C, tau0s, T0C, 
	      tune_fixed, tune_random, acc_fixed, acc_random, A, 
	      *max_samp_grp, betaA, T0A);

    /** Step 3: SAMPLE COMPLIANCE COVARIATE **/
    SampCompMixed(n_grp, n_samp, n_fixedC, Xc, betaC, Zc, grp, 
		  xiC, n_randomC, *AT, *logitC, qC, qN, Z, D, R, RD, 
		  prC, prN, Zo, Zr, C, Xo, Xr, *random, Xobs, Zobs, 
		  prA, pA, A, betaA, pC, pN);

    /** Step 4: OUTCOME MODEL **/
    bNormalMixedGibbs(Yobs, Xobs, Zobs, grp_obs, gamma, 
		      xiO, sig2, PsiO, n_obs, n_fixedO, 
		      n_randomO, n_grp, 0, gamma0, A0O, 
		      0, *nu0, *s0, tau0s[2], T0O, 1); 
    
    /** Compute probabilities of Y = Yobs **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 2; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (*random) {
	      pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], sqrt(*sig2), 0);
	      pA[i] = dnorm(Y[i], meano[i]+gamma[2]+xiO[grp[i]][1], sqrt(*sig2), 0);
	    } else {
	      pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]], sqrt(*sig2), 0);
	      pA[i] = dnorm(Y[i], meano[i]+gamma[2], sqrt(*sig2), 0);
	    }
	    pN[i] = dnorm(Y[i], meano[i], sqrt(*sig2), 0);
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 1; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (*random) 
	      pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], sqrt(*sig2), 0);
	    else 
	      pC[i] = dnorm(Y[i], meano[i]+gamma[1-Z[i]], sqrt(*sig2), 0);
	    pN[i] = dnorm(Y[i], meano[i], sqrt(*sig2), 0);
	  } 
      }
      vitemp[grp[i]]++;
    }
    
     /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	for (j = 0; j < (n_grp+1); j++) {
	  n_comp[j][0] = 0; n_comp[j][1] = 0;
	  n_never[j][0] = 0; n_never[j][1] = 0;
	  n_always[j][0] = 0; n_always[j][1] = 0;
	  p_comp[j] = 0; p_never[j] = 0; ITT[j] = 0;
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp[grp[i]] += qC[i]; p_comp[n_grp] += qC[i];
	  p_never[grp[i]] += qN[i]; p_never[n_grp] += qN[i];
	  /* counting */
	  if (C[i] == 1) { 
	    if (Z[i] == 1) {
	      n_comp[grp[i]][1]++; n_comp[n_grp][1]++;
	    } else {
	      n_comp[grp[i]][0]++; n_comp[n_grp][0]++;
	    }
	  } else if (A[i] == 1) {
	    if (Z[i] == 1) {
	      n_always[grp[i]][1]++; n_always[n_grp][1]++;
	    } else {
	      n_always[grp[i]][0]++; n_always[n_grp][0]++;
	    }
	  } else {
	    if (Z[i] == 1) {
	      n_never[grp[i]][1]++; n_never[n_grp][1]++;
	    } else {
	      n_never[grp[i]][0]++; n_never[n_grp][0]++;
	    }
	  }
	  /* insample QoI */
	  if (*Insample) { 
	    if (C[i] == 1) { /* compliers */
	      if (*random) {
		dtemp = rnorm(meano[i]+gamma[0]+xiO[grp[i]][0], sqrt(*sig2));
		dtemp1 = rnorm(meano[i]+gamma[1]+xiO[grp[i]][0], sqrt(*sig2));
	      } else { 
		dtemp = rnorm(meano[i]+gamma[0], sqrt(*sig2));
		dtemp1 = rnorm(meano[i]+gamma[1], sqrt(*sig2));
	      }
	      if (R[i] == 1) {
		if (Z[i] == 1) {
		  dtemp = Y[i]; 
		} else { 
		  dtemp1 = Y[i];
		}
	      }
	      Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	      Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    } else if (A[i] == 1) { /* always-takers */ 
	      if (R[i] == 1)
		dtemp = Y[i];
	      else if (*random) 
		dtemp = rnorm(meano[i]+gamma[2]+xiO[grp[i]][1], sqrt(*sig2));
	      else
		dtemp = rnorm(meano[i]+gamma[2], sqrt(*sig2));
	      YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    } else { /* never-takers */
	      if (R[i] == 1)
		dtemp = Y[i];
	      else
		dtemp = rnorm(meano[i], sqrt(*sig2));
	      YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    if (*random) {
	      dtemp = (meano[i]+gamma[0]+xiO[grp[i]][0]);
	      dtemp1 = (meano[i]+gamma[1]+xiO[grp[i]][0]); 
	    } else {
	      dtemp = (meano[i]+gamma[0]);
	      dtemp1 = (meano[i]+gamma[1]); 
	    }
	    Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	    Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    ITT[grp[i]] += (dtemp-dtemp1)*qC[i];
	    ITT[n_grp] += (dtemp-dtemp1)*qC[i];
	    /* always-takers */
	    if (*random)
	      dtemp = (meano[i]+gamma[2]+xiO[grp[i]][1]);
	    else
	      dtemp = (meano[i]+gamma[2]);
	    YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    /* never-takers */
	    dtemp = meano[i];
	    YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	  }
	}
	
	uniQoIcalMixed(*Insample, n_grp, ITT, Y1barC, Y0barC, n_samp, n_comp,
		       n_never, n_always, p_comp, p_never, CACE, YbarN,
		       YbarA, *AT);

	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_comp[j]; 	  
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_never[j];
	for (j = 0; j < (n_grp+1); j++) 
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = YbarN[j];
	if (*AT)
	  for (j = 0; j < (n_grp+1); j++)
	    QoI[itempQ++] = YbarA[j];

	if (*param) {
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  ssig2[itempS++] = sig2[0];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    } else {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
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
	  Rprintf("  Current Acceptance Ratio for fixed effects:");
	  if (*AT)
	    for (j = 0; j < n_fixedC*2; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_fixedC; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  Rprintf("\n");
	  Rprintf("  Current Acceptance Ratio for random effects:");
	  if (*AT)
	    for (j = 0; j < 2; j++)
	      Rprintf("%10g", (double)acc_random[j]/(double)main_loop);
	  else
	    Rprintf("%10g", (double)acc_random[0]/(double)main_loop);
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
  free(grp_obs);
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_fixedC);
  Free3DMatrix(Zc, n_grp, *max_samp_grp + n_randomC);
  FreeMatrix(Xo, n_samp+n_fixedO);
  Free3DMatrix(Zo, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xobs, n_obs+n_fixedO);
  Free3DMatrix(Zobs, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xr, n_samp+n_fixedR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_randomR);
  Free3DMatrix(xiC, 2, n_grp);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiR, n_grp);
  Free3DMatrix(Psi, 2, n_randomC);
  FreeMatrix(PsiO, n_randomO);
  FreeMatrix(PsiR, n_randomR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  free(ITT);
  free(CACE);
  free(Y1barC);
  free(Y0barC);
  free(YbarN);
  free(YbarA);
  FreeintMatrix(n_comp, n_grp+1);
  FreeintMatrix(n_never, n_grp+1);
  FreeintMatrix(n_always, n_grp+1);
  free(p_comp);
  free(p_never);
  free(vitemp);
  free(acc_fixed);
  free(acc_random);

} /* end of LINormalMixed */


/* 
   Bayesian binary probit mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LIboprobitMixed(int *Y,         /* binary outcome variable */ 
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
		     int *grp,       /* group indicator: 0, 1, 2, ... */
		     int *Ymiss,     /* number of missing obs in Y */
		     int *AT,        /* Are there always-takers? */
		     int *Insample,  /* Insample (=1) or population QoI? */
		     int *random,    /* random compliance status */
		     double *dXc,    /* fixed effects for compliance model */
		     double *dZc,    /* random effects for compliance model */
		     double *dXo,    /* fixed effects for outcome model */
		     double *dZo,    /* random effects for outcome model */
		     double *dXr,    /* fixed effects for response model */
		     double *dZr,    /* random effects for response model */
		     double *betaC,  /* fixed effects for compliance model */
		     double *betaA,  /* fixed effects for always-takers model */
		     double *gamma,  /* fixed effects for outcome
					model */
		     double *tau,    /* cutpoints */
		     double *delta,  /* fixed effects for response model */
		     int *in_samp,   /* # of observations */
		     int *in_cat,     /* # of categories */
		     int *n_gen,     /* # of Gibbs draws */
		     int *in_grp,    /* # of groups */
		     int *max_samp_grp, /* max # of obs within group */
		     int *in_fixed,  /* # of fixed effects for
					compliance, outcome, and response models */
		     int *in_random, /* # of random effects for
					compliance, outcome, and response modeld */ 
		     double *dPsiC,  /* precision for compliance model */
		     double *dPsiA,  /* precision for always-takers model */
		     double *dPsiO,  /* precision for outcome model */
		     double *dPsiR,  /* precision for response model */
		     double *beta0,  /* prior mean for betaC and betaA */ 
		     double *gamma0, /* prior mean for gamma */
		     double *delta0, /* prior mean for delta */
		     double *dA0C,   /* prior precision for betaC and betaA */ 
		     double *dA0O,   /* prior precision for gamma */
		     double *dA0R,   /* prior precision for delta */
		     int *tau0s,     /* prior df for PsiC, PsiA, PsiO, PsiR */
		     double *dT0C,   /* prior scale for PsiC */
		     double *dT0A,   /* prior scale for PsiA */
		     double *dT0O,   /* prior scale for PsiO */
		     double *dT0R,   /* prior scale for PsiR */
		     double *tune_fixed,  /* proposal variance */
		     double *tune_random, /* proposal variance */
		     double *tune_tau,    /* proposal variance */
		     int *mh,        /* Use MH step for taus? */
		     int *logitC,    /* Use logistic regression for the
					compliance model? */
		     int *param,     /* Want to keep paramters? */
		     int *burnin,    /* number of burnin */
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
		     double *tauO,
		     double *sPsiC,  /* Storage for precisions of
					random effects */
		     double *sPsiA,
		     double *sPsiO,
		     double *sPsiR, 
		     double *QoI     /* Storage of quantities of
					interest */		 
		     ) {
   /** counters **/
  int n_samp = *in_samp; 
  int n_grp = *in_grp;
  int n_cat = *in_cat;
  int n_fixedC = in_fixed[0]; int n_randomC = in_random[0];
  int n_fixedO = in_fixed[1]; int n_randomO = in_random[1];
  int n_fixedR = in_fixed[2]; int n_randomR = in_random[2];
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  int *grp_obs = intArray(n_obs);

  /*** observed Y ***/
  int *Yobs = intArray(n_obs);

  /* covariates for fixed effects in the compliance model */
  double **Xc = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  /* covariates for random effects */
  double ***Zc = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				n_randomC +1);

  /* covariates for fixed effects in the outcome model */
  double **Xo = doubleMatrix(n_samp+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				n_randomO +1);

  /* covariates for fixed effecs in the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zobs = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				  n_randomO +1);

  /* covariates for fixed effects in the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_fixedR, n_fixedR+1);    
  /* covariates for random effects */
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_randomR,
				n_randomR +1);

  /*** model parameters ***/
  /* random effects */
  double ***xiC = doubleMatrix3D(2, n_grp, n_randomC);
  double **xiO = doubleMatrix(n_grp, n_randomO);
  double **xiR = doubleMatrix(n_grp, n_randomR);

  /* covariances for random effects */
  double ***Psi = doubleMatrix3D(2, n_randomC, n_randomC);
  double **PsiO = doubleMatrix(n_randomO, n_randomO);
  double **PsiR = doubleMatrix(n_randomR, n_randomR);

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

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
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* quantities of interest: ITT, CACE  */
  double **ITT = doubleMatrix(n_grp+1, n_cat-1);
  double **CACE = doubleMatrix(n_grp+1, n_cat-1);
  double **Y1barC = doubleMatrix(n_grp+1, n_cat-1); 
  double **Y0barC = doubleMatrix(n_grp+1, n_cat-1); 
  double **YbarN = doubleMatrix(n_grp+1, n_cat-1);
  double **YbarA = doubleMatrix(n_grp+1, n_cat-1);
  int **n_comp = intMatrix(n_grp+1, 2);          /* number of compliers */
  int **n_never = intMatrix(n_grp+1, 2);
  int **n_always = intMatrix(n_grp+1, 2);
  double *p_comp = doubleArray(n_grp+1); 
  double *p_never = doubleArray(n_grp+1); /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int progress; progress = 1;
  int keep; keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);    /* number of acceptance */
  int *acc_tau = intArray(0);
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, k, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR, itempT;
  int itempAv, itempCv, itempOv, itempRv;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Data prep etc. ***/
  PrepMixed(dXc, dZc, dXo, dZo, dXr, dZr, Xc, Xo, Xr, Xobs, Zc, Zo,
	    Zr, n_samp, n_grp, n_obs, n_miss, n_fixedC, n_fixedO,
	    n_fixedR, n_randomC, n_randomO, n_randomR, Zobs, R,
	    grp, grp_obs, xiC, xiO, xiR, dPsiC, dPsiA, dPsiO, dPsiR,
	    Psi, PsiO, PsiR, dT0C, T0C, dT0A, T0A, dT0O, T0O, dT0R,
	    T0R, dA0C, A0C, dA0O, A0O, dA0R, A0R, *logitC, pC, pN,
	    pA, prC, prN, prA, *AT, beta0, gamma0, delta0, 1);

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0; itempT = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0; acc_tau[0] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      ResponseMixed(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, delta0, A0R, tau0s, T0R,
		    *AT, *random, Z, D, prC,prN, prA);
    
    /** Step 2: COMPLIANCE MODEL **/
    CompMixed(*logitC, *AT, C, Xc, Zc, grp, betaC, xiC, Psi, n_samp,
	      n_fixedC, n_randomC , n_grp, beta0, A0C, tau0s, T0C, 
	      tune_fixed, tune_random, acc_fixed, acc_random, A, 
	      *max_samp_grp, betaA, T0A);

    /** Step 3: SAMPLE COMPLIANCE COVARIATE **/
    SampCompMixed(n_grp, n_samp, n_fixedC, Xc, betaC, Zc, grp, 
		  xiC, n_randomC, *AT, *logitC, qC, qN, Z, D, R, RD, 
		  prC, prN, Zo, Zr, C, Xo, Xr, *random, Xobs, Zobs, 
		  prA, pA, A, betaA, pC, pN);

    /** Step 4: OUTCOME MODEL **/
    if (*mh && (main_loop == 1)) {
      itemp = 0;
      for (j = 0; j < n_grp; j++)
	vitemp[j] = 0;
      for (i = 0; i < n_samp; i++){
	if (R[i] == 1) {
	  dtemp = 0; dtemp1 = 0;
	  for (j = 0; j < n_fixedO; j++)
	    dtemp += Xo[i][j]*gamma[j];
	  for (j = 0; j < n_randomO; j++)
	    dtemp1 += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	  dtemp += dtemp1;
	  if (Y[i] == 0)
	    Xobs[itemp++][n_fixedO] = 
	      TruncNorm(dtemp-1000,0,dtemp,1,0) - dtemp1;
	  else
	    Xobs[itemp++][n_fixedO] =
	      TruncNorm(tau[Y[i]-1],tau[Y[i]],dtemp,1,0) - dtemp1;
	}
	vitemp[grp[i]]++;
      }
    }
    
    boprobitMixedMCMC(Yobs, Xobs, Zobs, grp_obs, gamma, xiO, tau, 
		      PsiO, n_obs, n_cat, n_fixedO, n_randomO, n_grp,
		      0, gamma0, A0O, tau0s[2], T0O, *mh, tune_tau,
		      acc_tau, 1);  
    
    /** Compute probabilities of Y = Yobs **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 2; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (*random) {
	      if (Y[i] == 0) {
		pC[i] = pnorm(tau[0], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0);
		pA[i] = pnorm(tau[0], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0);
	      } else {
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0);
		pA[i] = pnorm(tau[Y[i]], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0);
	      }
	    } else {
	      if (Y[i] == 0) {
		pC[i] = pnorm(tau[0], meano[i]+gamma[1-Z[i]], 1, 1, 0);
		pA[i] = pnorm(tau[0], meano[i]+gamma[2], 1, 1, 0);
	      } else {
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1-Z[i]], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1-Z[i]], 1, 1, 0);
		pA[i] = pnorm(tau[Y[i]], meano[i]+gamma[2], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[2], 1, 1, 0);
	      }
	    }
	    if (Y[i] == 0)
	      pN[i] = pnorm(tau[0], meano[i], 1, 1, 0);
	    else
	      pN[i] = pnorm(tau[Y[i]], meano[i], 1, 1, 0) -
		pnorm(tau[Y[i]-1], meano[i], 1, 1, 0);
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 1; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (*random)
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 1, 1, 0);
	    else
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1-Z[i]], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1-Z[i]], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1-Z[i]], 1, 1, 0);
	    if (Y[i] == 0)
	      pN[i] = pnorm(tau[0], meano[i], 1, 1, 0);
	    else
	      pN[i] = pnorm(tau[Y[i]], meano[i], 1, 1, 0)
		- pnorm(tau[Y[i]-1], meano[i], 1, 1, 0); 
	  }
      }
      vitemp[grp[i]]++;
    }
    
   /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	for (j = 0; j < (n_grp+1); j++) {
	  n_comp[j][0] = 0; n_comp[j][1] = 0;
	  n_never[j][0] = 0; n_never[j][1] = 0;
	  n_always[j][0] = 0; n_always[j][1] = 0;
	  p_comp[j] = 0; p_never[j] = 0;
	  for (k = 0; k < (n_cat-1); k++) {
	    Y1barC[j][k] = 0; Y0barC[j][k] = 0; 
	    YbarN[j][k] = 0; YbarA[j][k] = 0; ITT[j][k] = 0;
	  }
	}
	for (i = 0; i < n_samp; i++){
	  p_comp[grp[i]] += qC[i]; p_comp[n_grp] += qC[i];
	  p_never[grp[i]] += qN[i]; p_never[n_grp] += qN[i];
	  /* counting */
	  if (C[i] == 1) { 
	    if (Z[i] == 1) {
	      n_comp[grp[i]][1]++; n_comp[n_grp][1]++;
	    } else {
	      n_comp[grp[i]][0]++; n_comp[n_grp][0]++;
	    }
	  } else if (A[i] == 1) {
	    if (Z[i] == 1) {
	      n_always[grp[i]][1]++; n_always[n_grp][1]++;
	    } else {
	      n_always[grp[i]][0]++; n_always[n_grp][0]++;
	    }
	  } else {
	    if (Z[i] == 1) {
	      n_never[grp[i]][1]++; n_never[n_grp][1]++;
	    } else {
	      n_never[grp[i]][0]++; n_never[n_grp][0]++;
	    }
	  }
	  /* insample QoI */
	  for (j = 1; j < n_cat; j++) {
	    if (*Insample) { 
	      if (C[i] == 1) { /* compliers */
		if (*random) {
		  dtemp = (double)(unif_rand() < 
				   (pnorm(tau[j], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0)));
		  dtemp1 = (double)(unif_rand() < 
				    (pnorm(tau[j], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0) 
				     -pnorm(tau[j-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0)));
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
		Y1barC[grp[i]][j-1] += dtemp; Y0barC[grp[i]][j-1] += dtemp1;
		Y1barC[n_grp][j-1] += dtemp; Y0barC[n_grp][j-1] += dtemp1;
	      } else if (A[i] == 1) { /* always-takers */ 
		  if (R[i] == 1)
		    dtemp = (double)(Y[i] == j);
		  else if (*random)
		    dtemp = (double)(unif_rand() < 
				     (pnorm(tau[j], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0) 
				      -pnorm(tau[j-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0)));	
		  else
		    dtemp = (double)(unif_rand() < 
				     (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
				      -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0)));	
		YbarA[grp[i]][j-1] += dtemp; YbarA[n_grp][j-1] += dtemp;
	      } else { /* never-takers */
		  if (R[i] == 1)
		    dtemp = (double)(Y[i] == j);
		  else 
		    dtemp = (double)(unif_rand() < 
				     (pnorm(tau[j], meano[i], 1, 1, 0) 
				      -pnorm(tau[j-1], meano[i], 1, 1, 0)));		
		  YbarN[grp[i]][j-1] += dtemp; YbarN[n_grp][j-1] += dtemp;
	      } 
	    } else { /* population QoI */
	      /* compliers */
	      if (*random) {
		dtemp = (pnorm(tau[j], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0) 
			 -pnorm(tau[j-1], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0));
		dtemp1 = (pnorm(tau[j], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0) 
			  -pnorm(tau[j-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0));
	      } else {
		dtemp = (pnorm(tau[j], meano[i]+gamma[0], 1, 1, 0) 
			 -pnorm(tau[j-1], meano[i]+gamma[0], 1, 1, 0));
		dtemp1 = (pnorm(tau[j], meano[i]+gamma[1], 1, 1, 0) 
			  -pnorm(tau[j-1], meano[i]+gamma[1], 1, 1, 0));
	      }
	      Y1barC[grp[i]][j-1] += dtemp; Y0barC[grp[i]][j-1] += dtemp1;
	      Y1barC[n_grp][j-1] += dtemp; Y0barC[n_grp][j-1] += dtemp1;
	      ITT[grp[i]][j-1] += (dtemp-dtemp1)*qC[i];
	      ITT[n_grp][j-1] += (dtemp-dtemp1)*qC[i];
	      /* always-takers */
	      if (*random)
		dtemp = (pnorm(tau[j], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0) 
			 -pnorm(tau[j-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0));
	      else
		dtemp += (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
			  -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0));
	      YbarA[grp[i]][j-1] += dtemp; YbarA[n_grp][j-1] += dtemp;
	      /* never-takers */
	      dtemp = (pnorm(tau[j], meano[i], 1, 1, 0) 
		       -pnorm(tau[j-1], meano[i], 1, 1, 0));
	      YbarN[grp[i]][j-1] += dtemp; YbarN[n_grp][j-1] += dtemp;
	    }
	  }
	}

	if (*Insample) {
	  for (j = 0; j < (n_grp+1); j++) {
	    for (k = 0; k < (n_cat-1); k++) {
	      ITT[j][k] = (Y1barC[j][k]-Y0barC[j][k]) /
		(double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
			 n_comp[j][1]+n_never[j][1]+n_always[j][1]);
	      Y1barC[j][k] /= (double)(n_comp[j][0]+n_comp[j][1]);  
	      Y0barC[j][k] /= (double)(n_comp[j][0]+n_comp[j][1]); 
	      CACE[j][k] = Y1barC[j][k]-Y0barC[j][k];    
	      YbarN[j][k] /= (double)(n_never[j][0]+n_never[j][1]);
	      if (*AT)
		YbarA[j][k] /= (double)(n_always[j][0]+n_always[j][1]);
	    }
	    p_comp[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
				  n_comp[j][1]+n_never[j][1]+n_always[j][1]);  
	    p_never[j] /= (double)(n_comp[j][0]+n_never[j][0]+n_always[j][0] + 
				   n_comp[j][1]+n_never[j][1]+n_always[j][1]);  
	  }
	} else {
	  for (j = 0; j < (n_grp+1); j++) {
	    for (k = 0; k < (n_cat-1); k++) {
	      ITT[j][k] /= (double)n_samp;
	      Y1barC[j][k] /= (double)n_samp;  
	      Y0barC[j][k] /= (double)n_samp;
	      CACE[j][k] = Y1barC[j][k]-Y0barC[j][k];
	      YbarN[j][k] /= (double)n_samp;
	      if (*AT)
		YbarA[j][k] /= (double)n_samp;
	    }
	    p_comp[j] /= (double)n_samp;
	    p_never[j] /= (double)n_samp;
	  }
	}
	
	for (j = 0; j < (n_grp+1); j++)
	  for (k = 0; k < (n_cat-1); k++) 
	    QoI[itempQ++] = ITT[j][k];   
	for (j = 0; j < (n_grp+1); j++)
	  for (k = 0; k < (n_cat-1); k++) 
	    QoI[itempQ++] = CACE[j][k];   
	for (j = 0; j < (n_grp+1); j++)
	  for (k = 0; k < (n_cat-1); k++) 
	    QoI[itempQ++] = Y1barC[j][k];
	for (j = 0; j < (n_grp+1); j++)
	  for (k = 0; k < (n_cat-1); k++) 
	    QoI[itempQ++] = Y0barC[j][k];
	for (j = 0; j < (n_grp+1); j++)
	  for (k = 0; k < (n_cat-1); k++) 
	    QoI[itempQ++] = YbarN[j][k];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_comp[j]; 	  
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_never[j];
	if (*AT)
	  for (j = 0; j < (n_grp+1); j++)
	    for (k = 0; k < (n_cat-1); k++) 
	      QoI[itempQ++] = YbarA[j][k];

	if (*param) {
	  for (j = 0; j < (n_cat-1); j++)
	    tauO[itempT++] = tau[j];
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    } else {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
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
	  Rprintf("  Current Acceptance Ratio for fixed effects:");
	  if (*AT)
	    for (j = 0; j < n_fixedC*2; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_fixedC; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  Rprintf("\n");
	  Rprintf("  Current Acceptance Ratio for random effects:");
	  if (*AT)
	    for (j = 0; j < 2; j++)
	      Rprintf("%10g", (double)acc_random[j]/(double)main_loop);
	  else
	    Rprintf("%10g", (double)acc_random[0]/(double)main_loop);
	  Rprintf("\n");
	} 
	Rprintf("  Current Acceptance Ratio for cut points:");
	Rprintf("%10g\n", (double)acc_tau[0]/(double)main_loop);
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
  free(grp_obs);
  free(Yobs);
  FreeMatrix(Xc, n_samp+n_fixedC);
  Free3DMatrix(Zc, n_grp, *max_samp_grp + n_randomC);
  FreeMatrix(Xo, n_samp+n_fixedO);
  Free3DMatrix(Zo, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xobs, n_obs+n_fixedO);
  Free3DMatrix(Zobs, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xr, n_samp+n_fixedR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_randomR);
  Free3DMatrix(xiC, 2, n_grp);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiR, n_grp);
  Free3DMatrix(Psi, 2, n_randomC);
  FreeMatrix(PsiO, n_randomO);
  FreeMatrix(PsiR, n_randomR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  FreeMatrix(ITT, n_grp+1);
  FreeMatrix(CACE, n_grp+1);
  FreeMatrix(Y1barC, n_grp+1);
  FreeMatrix(Y0barC, n_grp+1);
  FreeMatrix(YbarN, n_grp+1);
  FreeMatrix(YbarA, n_grp+1);
  FreeintMatrix(n_comp, n_grp+1);
  FreeintMatrix(n_never, n_grp+1);
  FreeintMatrix(n_always, n_grp+1);
  free(p_comp);
  free(p_never);
  free(vitemp);
  free(acc_fixed);
  free(acc_random);
  free(acc_tau);

} /* end of LIboprobitMixed */


/* 
   Bayesian Normal mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LINegBinMixed(int *Y,         /* count outcome variable */ 
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
		   int *grp,       /* group indicator: 0, 1, 2, ... */
		   int *Ymiss,     /* number of missing obs in Y */
		   int *AT,        /* Are there always-takers? */
		   int *Insample,  /* Insample (=1) or population QoI? */
		   int *random,    /* Want to include compliance
				      random effects? */
		   double *dXc,    /* fixed effects for compliance model */
		   double *dZc,    /* random effects for compliance model */
		   double *dXo,    /* fixed effects for outcome model */
		   double *dZo,    /* random effects for outcome model */
		   double *dXr,    /* fixed effects for response model */
		   double *dZr,    /* random effects for response model */
		   double *betaC,  /* fixed effects for compliance model */
		   double *betaA,  /* fixed effects for always-takers model */
		   double *gamma,  /* fixed effects for outcome model */
		   double *delta,  /* fixed effects for response model */
		   double *sig2,   /* dispersion parameter for outcome model */
		   int *in_samp,   /* # of observations and # of groups */
		   int *n_gen,     /* # of Gibbs draws */
		   int *max_samp_grp, /* max # of obs within group */
		   int *in_fixed,  /* # of fixed effects for
				      compliance, outcome, and response models */
		   int *in_random, /* # of random effects for
				      compliance, outcome, and response modeld */ 
		   double *dPsiC,  /* precision for compliance model */
		   double *dPsiA,  /* precision for always-takers model */
		   double *dPsiO,  /* precision for outcome model */
		   double *dPsiR,  /* precision for response model */
		   double *beta0,  /* prior mean for betaC and betaA */ 
		   double *gamma0, /* prior mean for gamma */
		   double *delta0, /* prior mean for delta */
		   double *dA0C,   /* prior precision for betaC and betaA */ 
		   double *dA0O,   /* prior precision for gamma */
		   double *dA0R,   /* prior precision for delta */
		   double *a0,     /* prior shape for sig2 */
		   double *b0,     /* prior scale for sig2 */
		   int *tau0s,     /* prior df for PsiC, PsiA, PsiO, PsiR */
		   double *dT0C,   /* prior scale for PsiC */
		   double *dT0A,   /* prior scale for PsiA */
		   double *dT0O,   /* prior scale for PsiO */
		   double *dT0R,   /* prior scale for PsiR */
		   double *tune_fixed, /* proposal variance */
		   double *tune_random, /* proposal variance */
		   double *varb,   /* proposal variance for beta */
		   double *varg,   /* proposal variance for gamma */
		   double *vars,   /* proposal variance for sig2 */
		   int *logitC,    /* Use logistic regression for the
				      compliance model? */
		   int *param,     /* Want to keep paramters? */
		   int *burnin,    /* number of burnin */
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
		   double *ssig2,  /* Storage for outcome variance */
		   double *sPsiC,  /* Storage for precisions of
				      random effects */
		   double *sPsiA,
		   double *sPsiO,
		   double *sPsiR, 
		   double *QoI     /* Storage of quantities of
				      interest */		 
		   ) {
   /** conuters **/
  int n_samp = in_samp[0]; 
  int n_grp = in_samp[1];
  int n_fixedC = in_fixed[0]; int n_randomC = in_random[0];
  int n_fixedO = in_fixed[1]; int n_randomO = in_random[1];
  int n_fixedR = in_fixed[2]; int n_randomR = in_random[2];
  int n_miss = *Ymiss;
  int n_obs = n_samp - n_miss;

  /*** data ***/
  int *grp_obs = intArray(n_obs);

  /*** observed Y ***/
  int *Yobs = intArray(n_obs);
  int **Ygrp = intMatrix(n_grp, *max_samp_grp);

  /* covariates for fixed effects in the compliance model */
  double **Xc = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  /* covariates for random effects */
  double ***Zc = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				n_randomC +1);

  /* covariates for fixed effects in the outcome model */
  double **Xo = doubleMatrix(n_samp, n_fixedO);    
  /* covariates for random effects */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp, n_randomO);

  /* covariates for fixed effecs in the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs, n_fixedO);    
  /* covariates for random effects */
  double ***Zobs = doubleMatrix3D(n_grp, *max_samp_grp, n_randomO);

  /* covariates for fixed effects in the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_fixedR, n_fixedR+1);    
  /* covariates for random effects */
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_randomR,
				n_randomR +1);

  /*** model parameters ***/
  /* random effects */
  double ***xiC = doubleMatrix3D(2, n_grp, n_randomC);
  double **xiO = doubleMatrix(n_grp, n_randomO);
  double **xiR = doubleMatrix(n_grp, n_randomR);

  /* covariances for random effects */
  double ***Psi = doubleMatrix3D(2, n_randomC, n_randomC);
  double **PsiO = doubleMatrix(n_randomO, n_randomO);
  double **PsiR = doubleMatrix(n_randomR, n_randomR);

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);

  /* density of Y = Yobs for a complier */
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
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* qoi by groups */
  double *ITT = doubleArray(n_grp+1);
  double *CACE = doubleArray(n_grp+1);
  double *Y1barC = doubleArray(n_grp+1);
  double *Y0barC = doubleArray(n_grp+1);
  double *YbarN = doubleArray(n_grp+1); 
  double *YbarA = doubleArray(n_grp+1);
  int **n_comp = intMatrix(n_grp+1, 2);          /* number of compliers */
  int **n_never = intMatrix(n_grp+1, 2);
  int **n_always = intMatrix(n_grp+1, 2);
  double *p_comp = doubleArray(n_grp+1); 
  double *p_never = doubleArray(n_grp+1); /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int progress; progress = 1;
  int keep; keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int *counter = intArray(2); 
  int **counterg = intMatrix(n_grp, 2);
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv, itempS;
  double dtemp, dtemp1;

  /*** get random seed **/
  GetRNGstate();

  /*** Data prep etc. ***/
  PrepMixed(dXc, dZc, dXo, dZo, dXr, dZr, Xc, Xo, Xr, Xobs, Zc, Zo,
	    Zr, n_samp, n_grp, n_obs, n_miss, n_fixedC, n_fixedO,
	    n_fixedR, n_randomC, n_randomO, n_randomR, Zobs, R,
	    grp, grp_obs, xiC, xiO, xiR, dPsiC, dPsiA, dPsiO, dPsiR,
	    Psi, PsiO, PsiR, dT0C, T0C, dT0A, T0A, dT0O, T0O, dT0R,
	    T0R, dA0C, A0C, dA0O, A0O, dA0R, A0R, *logitC, pC, pN,
	    pA, prC, prN, prA, *AT, beta0, gamma0, delta0, 0);

  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) {
      Yobs[itemp++] = Y[i];
      Ygrp[grp[i]][vitemp[grp[i]]++] = Y[i];
    }

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0; itempS = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0; 
  counter[0] = 0; counter[1] = 0;
  for (j = 0; j < n_grp; j++) {
    counterg[j][0] = 0; counterg[j][1] = 0;
  }
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      ResponseMixed(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, delta0, A0R, tau0s, T0R,
		    *AT, *random, Z, D, prC,prN, prA);

    /** Step 2: COMPLIANCE MODEL **/
    CompMixed(*logitC, *AT, C, Xc, Zc, grp, betaC, xiC, Psi, n_samp,
	      n_fixedC, n_randomC , n_grp, beta0, A0C, tau0s, T0C, 
	      tune_fixed, tune_random, acc_fixed, acc_random, A, 
	      *max_samp_grp, betaA, T0A);

    /** Step 3: SAMPLE COMPLIANCE COVARIATE **/
    SampCompMixed(n_grp, n_samp, n_fixedC, Xc, betaC, Zc, grp, 
		  xiC, n_randomC, *AT, *logitC, qC, qN, Z, D, R, RD,
		  prC, prN, Zo, Zr, C, Xo, Xr, *random, Xobs, Zobs, 
		  prA, pA, A, betaA, pC, pN);

    /** Step 4: OUTCOME MODEL **/
    bnegbinMixedMCMC(Yobs, Ygrp, Xobs, Zobs, grp_obs, gamma, 
		     xiO, sig2, PsiO, n_obs, n_fixedO, 
		     n_randomO, n_grp, *max_samp_grp, gamma0, A0O, 
		     *a0, *b0, tau0s[2], T0O, varb, *vars, varg, 
		     counter, counterg, 1); 
    
    /** Compute probabilities of Y = Yobs **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 2; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (*random) {
	      pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0]), *sig2, 0);
	      pA[i] = dnegbin(Y[i], exp(meano[i]+gamma[2]+xiO[grp[i]][1]), *sig2, 0);
	    } else {
	      pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]), *sig2, 0);
	      pA[i] = dnegbin(Y[i], exp(meano[i]+gamma[2]), *sig2, 0);
	    }
	    pN[i] = dnegbin(Y[i], exp(meano[i]), *sig2, 0);
	  } 
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (*random)
	  for (j = 1; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	else
	  for (j = 0; j < n_randomO; j++)
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (*random) 
	      pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0]), *sig2, 0);
	    else 
	      pC[i] = dnegbin(Y[i], exp(meano[i]+gamma[1-Z[i]]), *sig2, 0);
	    pN[i] = dnegbin(Y[i], exp(meano[i]), *sig2, 0);
	  } 
      }
      vitemp[grp[i]++];
    }
 
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	for (j = 0; j < (n_grp+1); j++) {
	  n_comp[j][0] = 0; n_comp[j][1] = 0;
	  n_never[j][0] = 0; n_never[j][1] = 0;
	  n_always[j][0] = 0; n_always[j][1] = 0;
	  p_comp[j] = 0; p_never[j] = 0; ITT[j] = 0;
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp[grp[i]] += qC[i]; p_comp[n_grp] += qC[i];
	  p_never[grp[i]] += qN[i]; p_never[n_grp] += qN[i];
	  /* counting */
	  if (C[i] == 1) { 
	    if (Z[i] == 1) {
	      n_comp[grp[i]][1]++; n_comp[n_grp][1]++;
	    } else {
	      n_comp[grp[i]][0]++; n_comp[n_grp][0]++;
	    }
	  } else if (A[i] == 1) {
	    if (Z[i] == 1) {
	      n_always[grp[i]][1]++; n_always[n_grp][1]++;
	    } else {
	      n_always[grp[i]][0]++; n_always[n_grp][0]++;
	    }
	  } else {
	    if (Z[i] == 1) {
	      n_never[grp[i]][1]++; n_never[n_grp][1]++;
	    } else {
	      n_never[grp[i]][0]++; n_never[n_grp][0]++;
	    }
	  }
	  /* insample QoI */
	  if (*Insample) { 
	    if (C[i] == 1) { /* compliers */
	      if (*random) {
		dtemp = rnegbin(exp(meano[i]+gamma[0]+xiO[grp[i]][0]), *sig2);
		dtemp1 = rnegbin(exp(meano[i]+gamma[1]+xiO[grp[i]][0]), *sig2);
	      } else { 
		dtemp = rnegbin(exp(meano[i]+gamma[0]), *sig2);
		dtemp1 = rnegbin(exp(meano[i]+gamma[1]), *sig2);
	      }
	      if (R[i] == 1) {
		if (Z[i] == 1) { 
		  dtemp = (double)Y[i]; 
		} else { 
		  dtemp1 = (double)Y[i];
		}
	      }
	      Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	      Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    } else if (A[i] == 1) { /* always-takers */ 
	      if (R[i] == 1)
		dtemp = (double)Y[i];
	      else if (*random) 
		dtemp = rnegbin(exp(meano[i]+gamma[2]+xiO[grp[i]][1]), *sig2);
	      else
		dtemp = rnegbin(exp(meano[i]+gamma[2]), *sig2);
	      YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    } else { /* never-takers */
	      if (R[i] == 1)
		dtemp = (double)Y[i];
	      else
		dtemp = rnegbin(exp(meano[i]), *sig2);
	      YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    if (*random) {
	      dtemp = exp(meano[i]+gamma[0]+xiO[grp[i]][0]);
	      dtemp1 = exp(meano[i]+gamma[1]+xiO[grp[i]][0]);
	    } else {
	      dtemp = exp(meano[i]+gamma[0]);
	      dtemp1 = exp(meano[i]+gamma[1]);
	    }
	    Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	    Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    ITT[grp[i]] += (dtemp-dtemp1)*qC[i];
	    ITT[n_grp] += (dtemp-dtemp1)*qC[i];
	    /* always-takers */
	    if (*random)
	      dtemp = exp(meano[i]+gamma[2]+xiO[grp[i]][1]);
	    else
	      dtemp = exp(meano[i]+gamma[2]);
	    YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    /* never-takers */
	    dtemp = exp(meano[i]);
	    YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	  }
	}
	
	uniQoIcalMixed(*Insample, n_grp, ITT, Y1barC, Y0barC, n_samp, n_comp,
		       n_never, n_always, p_comp, p_never, CACE, YbarN,
		       YbarA, *AT);
   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_comp[j]; 	  
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_never[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = YbarN[j];
	if (*AT)
	  for (j = 0; j < (n_grp+1); j++)
	    QoI[itempQ++] = YbarA[j];

	if (*param) {
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  ssig2[itempS++] = sig2[0];
	  if (*AT) {
	    if (*logitC) {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    } else {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
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
	  Rprintf("  Acceptance ratio for fixed effects in the compliance model:");
	  if (*AT)
	    for (j = 0; j < n_fixedC*2; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_fixedC; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  Rprintf("\n");
	  Rprintf("  Acceptance ratio for random effectsin the compliance model:");
	  if (*AT)
	    for (j = 0; j < 2; j++)
	      Rprintf("%10g", (double)acc_random[j]/(double)main_loop);
	  else
	    Rprintf("%10g", (double)acc_random[0]/(double)main_loop);
	  Rprintf("\n");
	}
	Rprintf("  Acceptance ratio for the outcome model:");
	Rprintf("%10g%10g%10g\n", (double)counter[0]/(double)main_loop, 
		(double)counter[1]/(double)main_loop, 
		(double)counterg[0][0]/(double)main_loop);
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
  free(grp_obs);
  free(Yobs);
  FreeintMatrix(Ygrp, n_grp);
  FreeMatrix(Xc, n_samp+n_fixedC);
  Free3DMatrix(Zc, n_grp, *max_samp_grp + n_randomC);
  FreeMatrix(Xo, n_samp);
  Free3DMatrix(Zo, n_grp, *max_samp_grp);
  FreeMatrix(Xobs, n_obs);
  Free3DMatrix(Zobs, n_grp, *max_samp_grp);
  FreeMatrix(Xr, n_samp+n_fixedR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_randomR);
  Free3DMatrix(xiC, 2, n_grp);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiR, n_grp);
  Free3DMatrix(Psi, 2, n_randomC);
  FreeMatrix(PsiO, n_randomO);
  FreeMatrix(PsiR, n_randomR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  free(ITT);
  free(CACE);
  free(Y1barC);
  free(Y0barC);
  free(YbarN);
  free(YbarA);
  FreeintMatrix(n_comp, n_grp+1);
  FreeintMatrix(n_never, n_grp+1);
  FreeintMatrix(n_always, n_grp+1);
  free(p_comp);
  free(p_never);
  free(vitemp);
  free(acc_fixed);
  free(acc_random);
  free(counter);
  FreeintMatrix(counterg, n_grp);

} /* end of LINegBinMixed */


/* 
   Bayesian twopart mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LItwopartMixed(int *Y,         /* indicator for Y > 0 */
		    double *Y1,     /* Gaussian outcome variable for
				       Y > 0 */ 
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
		    int *grp,       /* group indicator: 0, 1, 2, ... */
		    int *Ymiss,     /* number of missing obs in Y and
				       Are there always takers? */
		    int *Insample,  /* Insample (=1) or population QoI? */
		    int *random,    /* Want to include compliance
				       random effects? */
		    double *dXc,    /* fixed effects for compliance model */
		    double *dZc,    /* random effects for compliance model */
		    double *dXo,    /* fixed effects for outcome model */
		    double *dZo,    /* random effects for outcome model */
		    double *dXr,    /* fixed effects for response model */
		    double *dZr,    /* random effects for response model */
		    double *betaC,  /* fixed effects for compliance model */
		    double *betaA,  /* fixed effects for always-takers model */
		    double *gamma,  /* fixed effects for outcome model
				       1 */
		    double *delta,  /* fixed effects for response model */
		    double *sig2,   /* variance parameter for outcome
				       model */
		    int *in_samp,   /* # of observations (a vector
				       where the first element is total number of obs and
				       the second element is the number of obs with Y > 0
				    */
		    int *n_gen,     /* # of Gibbs draws */
		    int *in_grp,    /* # of groups */
		    int *max_samp_grp, /* max # of obs within group */
		    int *in_fixed,  /* # of fixed effects for
				       compliance, outcome, and response models */
		    int *in_random, /* # of random effects for
				       compliance, outcome, and response modeld */ 
		    double *dPsiC,  /* precision for compliance model */
		    double *dPsiA,  /* precision for always-takers model */
		    double *dPsiO,  /* precision for outcome model */
		    double *dPsiR,  /* precision for response model */
		    double *beta0,  /* prior mean for betaC and betaA */ 
		    double *gamma0, /* prior mean for gamma */
		    double *delta0, /* prior mean for delta */
		    double *dA0C,   /* prior precision for betaC and betaA */ 
		    double *dA0O,   /* prior precision for gamma */
		    double *dA0R,   /* prior precision for delta */
		    int *nu0,       /* prior df for sig2 */
		    double *s0,     /* prior scale for sig2 */
		    int *tau0s,     /* prior df for PsiC, PsiA, PsiO, PsiR */
		    double *dT0C,   /* prior scale for PsiC */
		    double *dT0A,   /* prior scale for PsiA */
		    double *dT0O,   /* prior scale for PsiO */
		    double *dT0R,   /* prior scale for PsiR */
		    double *tune_fixed, /* proposal variance */
		    double *tune_random, /* proposal variance */
		    int *logitC,    /* Use logistic regression for the
				       compliance model? */
		    int *param,     /* Want to keep paramters? */
		    int *burnin,    /* number of burnin */
		    int *iKeep,     /* keep ?th draws */
		    int *verbose,   /* print out messages */
		    double *coefC,  /* Storage for coefficients of the
				       compliance model */
		    double *coefA,  /* Storage for coefficients of the
				       always-takers model */
		    double *coefO,  /* Storage for coefficients of the
				       outcome model */
		    double *coefO1,  /* Storage for coefficients of the
					outcome model */
		    double *coefR,  /* Storage for coefficients of the
				       response model */	
		    double *ssig2,  /* Storage for outcome variance */
		    double *sPsiC,  /* Storage for precisions of
				       random effects */
		    double *sPsiA,
		    double *sPsiO,
		    double *sPsiO1,
		    double *sPsiR, 
		    double *QoI     /* Storage of quantities of
				       interest */		 
		    ) {
  /** counters **/
  int n_samp = in_samp[0]; 
  int n_samp1 = in_samp[1]; 
  int n_grp = *in_grp;
  int n_fixedC = in_fixed[0]; int n_randomC = in_random[0];
  int n_fixedO = in_fixed[1]; int n_randomO = in_random[1];
  int n_fixedR = in_fixed[2]; int n_randomR = in_random[2];
  int n_miss = Ymiss[0];
  int AT = Ymiss[1];
  int n_obs = n_samp - n_miss;

  double *gamma1 = doubleArray(n_fixedO); /* fixed effects for outcome model
					     2 */

  /*** data ***/
  int *grp_obs = intArray(n_obs);
  int *grp_obs1 = intArray(n_samp1);

  /*** observed Y ***/
  int *Yobs = intArray(n_obs);
  double *Yobs1 = doubleArray(n_samp1);

  /* covariates for fixed effects in the compliance model */
  double **Xc = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  /* covariates for random effects */
  double ***Zc = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				n_randomC +1);

  /* covariates for fixed effects in the outcome model */
  double **Xo = doubleMatrix(n_samp+n_fixedO, n_fixedO+1);    
  /* covariates for random effects */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				n_randomO +1);

  /* covariates for fixed effecs in the outcome model: only units with observed Y */     
  double **Xobs = doubleMatrix(n_obs+n_fixedO, n_fixedO+1);    
  double **Xobs1 = doubleMatrix(n_samp1+n_fixedO, n_fixedO+1);    

  /* covariates for random effects */
  double ***Zobs = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				  n_randomO +1);
  double ***Zobs1 = doubleMatrix3D(n_grp, *max_samp_grp + n_randomO,
				   n_randomO +1);

  /* covariates for fixed effects in the response model: includes all obs */     
  double **Xr = doubleMatrix(n_samp+n_fixedR, n_fixedR+1);    
  /* covariates for random effects */
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_randomR,
				n_randomR +1);

  /*** model parameters ***/
  /* random effects */
  double ***xiC = doubleMatrix3D(2, n_grp, n_randomC);
  double **xiO = doubleMatrix(n_grp, n_randomO);
  double **xiO1 = doubleMatrix(n_grp, n_randomO);
  double **xiR = doubleMatrix(n_grp, n_randomR);

  /* covariances for random effects */
  double ***Psi = doubleMatrix3D(2, n_randomC, n_randomC);
  double **PsiO = doubleMatrix(n_randomO, n_randomO);
  double **PsiO1 = doubleMatrix(n_randomO, n_randomO);
  double **PsiR = doubleMatrix(n_randomR, n_randomR);

  /* density of Y = Yobs for complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);

  /* probability of R = 1 */
  double *prC = doubleArray(n_samp);
  double *prN = doubleArray(n_samp);
  double *prA = doubleArray(n_samp);

  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);

  /* density of Y = Yobs for always-taker */
  double *pA = doubleArray(n_samp);

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);
  double *meano1 = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* qoi by groups */
  double *ITT = doubleArray(n_grp+1);
  double *CACE = doubleArray(n_grp+1);
  double *Y1barC = doubleArray(n_grp+1);
  double *Y0barC = doubleArray(n_grp+1);
  double *YbarN = doubleArray(n_grp+1); 
  double *YbarA = doubleArray(n_grp+1);
  int **n_comp = intMatrix(n_grp+1, 2);          /* number of compliers */
  int **n_never = intMatrix(n_grp+1, 2);
  int **n_always = intMatrix(n_grp+1, 2);
  double *p_comp = doubleArray(n_grp+1); 
  double *p_never = doubleArray(n_grp+1); /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int progress; progress = 1;
  int keep; keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempO1, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv, itempS;
  double dtemp, dtemp1;
  double **mtemp = doubleMatrix(n_fixedO, n_fixedO);
  int *vitemp1 = intArray(n_grp);

  /*** get random seed **/
  GetRNGstate();

  /*** Data prep etc. ***/
  PrepMixed(dXc, dZc, dXo, dZo, dXr, dZr, Xc, Xo, Xr, Xobs, Zc, Zo,
	    Zr, n_samp, n_grp, n_obs, n_miss, n_fixedC, n_fixedO,
	    n_fixedR, n_randomC, n_randomO, n_randomR, Zobs, R,
	    grp, grp_obs, xiC, xiO, xiR, dPsiC, dPsiA, dPsiO, dPsiR,
	    Psi, PsiO, PsiR, dT0C, T0C, dT0A, T0A, dT0O, T0O, dT0R,
	    T0R, dA0C, A0C, dA0O, A0O, dA0R, A0R, *logitC, pC, pN,
	    pA, prC, prN, prA, AT, beta0, gamma0, delta0, 1);

  itemp = 0;
  for (j = 0; j < n_grp; j++) {
    vitemp[j] = 0; 
    vitemp1[j] = 0;
  }
  for (i = 0; i < n_samp; i++) {
    if ((R[i] == 1) && (Y[i] == 1)) {
      for (j = 0; j < n_fixedO; j++) 
	Xobs1[itemp][j] = Xo[i][j];
      grp_obs1[itemp] = grp[i];
      Yobs1[itemp] = log(Y1[i]);
      Xobs1[itemp++][n_fixedO] = log(Y1[i]);
      for (j = 0; j < n_randomO; j++)
	Zobs1[grp[i]][vitemp1[grp[i]]][j] = Zo[grp[i]][vitemp[grp[i]]][j];
      vitemp1[grp[i]]++;
    }
    vitemp[grp[i]]++;
  }

  dcholdc(A0O, n_fixedO, mtemp);
  for (i = 0; i < n_fixedO; i++) {
    Xobs1[itemp+i][n_fixedO] = 0;
    for (j = 0; j < n_fixedO; j++) {
      Xobs1[itemp+i][n_fixedO] += mtemp[i][j]*gamma0[j];
      Xobs1[itemp+i][j] = mtemp[i][j];
    }
  }

  for (j = 0; j < n_randomO; j++) {
    for (i = 0; i < n_randomO; i++)
      PsiO1[j][i] = PsiO[j][i];
    for (i = 0; i < n_grp; i++)
      xiO1[i][j] = xiO[i][j];
  }

  itemp = 0;
  for (i = 0; i < n_samp; i++) 
    if (R[i] == 1) 
      Yobs[itemp++] = Y[i];

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempO1 = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0; itempS = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0)
      ResponseMixed(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
		    n_fixedR, n_randomR, n_grp, delta0, A0R, tau0s, T0R,
		    AT, *random, Z, D, prC,prN, prA);
    
    /** Step 2: COMPLIANCE MODEL **/
    CompMixed(*logitC, AT, C, Xc, Zc, grp, betaC, xiC, Psi, n_samp,
	      n_fixedC, n_randomC , n_grp, beta0, A0C, tau0s, T0C, 
	      tune_fixed, tune_random, acc_fixed, acc_random, A, 
	      *max_samp_grp, betaA, T0A);

    /** Step 3: SAMPLE COMPLIANCE COVARIATE **/
    SampCompMixed(n_grp, n_samp, n_fixedC, Xc, betaC, Zc, grp, 
		  xiC, n_randomC, AT, *logitC, qC, qN, Z, D, R, RD, 
		  prC, prN, Zo, Zr, C, Xo, Xr, *random, Xobs, Zobs, 
		  prA, pA, A, betaA, pC, pN);

    /** Step 4: OUTCOME MODEL **/
    bprobitMixedGibbs(Yobs, Xobs, Zobs, grp_obs, gamma, xiO, PsiO,
		      n_obs, n_fixedO, n_randomO, n_grp, 0,
		      gamma0, A0O, tau0s[2], T0O, 1);
    bNormalMixedGibbs(Yobs1, Xobs1, Zobs1, grp_obs1, gamma1, 
		      xiO1, sig2, PsiO1, n_samp1, n_fixedO, 
		      n_randomO, n_grp, 0, gamma0, A0O, 
		      0, *nu0, *s0, tau0s[2], T0O, 1); 

    /** Compute probabilities of Y = Yobs **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0; meano1[i] = 0;
      if (AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++) {
	  meano[i] += Xo[i][j]*gamma[j];
	  meano1[i] += Xo[i][j]*gamma1[j];
	}
	if (*random)
	  for (j = 2; j < n_randomO; j++) {
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	    meano1[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO1[grp[i]][j];
	  }
	else
	  for (j = 0; j < n_randomO; j++) {
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	    meano1[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO1[grp[i]][j];
	  }
	if (R[i] == 1) {
	  if ((RD[i] == 0) || (Z[i] == D[i])) {
	    if (Y[i] == 1) {
	      if (*random) {
		pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]]+xiO1[grp[i]][0], sqrt(*sig2), 0) *
		  pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 1, 0);
		pA[i] = dlnorm(Y1[i], meano1[i]+gamma1[2]+xiO1[grp[i]][1], sqrt(*sig2), 0) * 
		  pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0);
	      } else {
		pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]], sqrt(*sig2), 0) * 
		  pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0);
		pA[i] = dlnorm(Y1[i], meano1[i]+gamma1[2], sqrt(*sig2), 0) *
		  pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	      }
	      pN[i] = dlnorm(Y1[i], meano1[i], sqrt(*sig2), 0) * 
		pnorm(meano[i], 0, 1, 1, 0);
	    } else {
	      if (*random) {
		pC[i] = pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 0, 0);
		pA[i] = pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 0, 0);
	      } else {
		pC[i] = pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
		pA[i] = pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	      }
	      pN[i] = pnorm(meano[i], 0, 1, 0, 0);
	    }
	  }
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++) {
	  meano[i] += Xo[i][j]*gamma[j];
	  meano1[i] += Xo[i][j]*gamma1[j];
	}
	if (*random)
	  for (j = 1; j < n_randomO; j++) {
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	    meano1[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO1[grp[i]][j];
	  }
	else
	  for (j = 0; j < n_randomO; j++) {
	    meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	    meano1[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO1[grp[i]][j];
	  }
	if (R[i] == 1)
	  if ((Z[i] == 0) || (RD[i] == 0)) {
	    if (Y[i] == 1) {
	      if (*random) 
		pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]]+xiO1[grp[i]][0], sqrt(*sig2), 0) * 
		  pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 1, 0);
	      else 
		pC[i] = dlnorm(Y1[i], meano1[i]+gamma1[1-Z[i]], sqrt(*sig2), 0) *
		  pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 1, 0);
	      pN[i] = dlnorm(Y1[i], meano1[i], sqrt(*sig2), 0) *
		pnorm(meano[i], 0, 1, 1, 0);
	    } else {
	      if (*random) 
		pC[i] = pnorm(meano[i]+gamma[1-Z[i]]+xiO[grp[i]][0], 0, 1, 0, 0);
	      else 
		pC[i] = pnorm(meano[i]+gamma[1-Z[i]], 0, 1, 0, 0);
	      pN[i] =  pnorm(meano[i], 0, 1, 1, 0);
	    } 
	  }
      }
      vitemp[grp[i]]++;
    }
    
    /** storing the results **/
    if (main_loop > *burnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	for (j = 0; j < (n_grp+1); j++) {
	  n_comp[j][0] = 0; n_comp[j][1] = 0;
	  n_never[j][0] = 0; n_never[j][1] = 0;
	  n_always[j][0] = 0; n_always[j][1] = 0;
	  p_comp[j] = 0; p_never[j] = 0; ITT[j] = 0;
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp[grp[i]] += qC[i]; p_comp[n_grp] += qC[i];
	  p_never[grp[i]] += qN[i]; p_never[n_grp] += qN[i];
	  /* counting */
	  if (C[i] == 1) { 
	    if (Z[i] == 1) {
	      n_comp[grp[i]][1]++; n_comp[n_grp][1]++;
	    } else {
	      n_comp[grp[i]][0]++; n_comp[n_grp][0]++;
	    }
	  } else if (A[i] == 1) {
	    if (Z[i] == 1) {
	      n_always[grp[i]][1]++; n_always[n_grp][1]++;
	    } else {
	      n_always[grp[i]][0]++; n_always[n_grp][0]++;
	    }
	  } else {
	    if (Z[i] == 1) {
	      n_never[grp[i]][1]++; n_never[n_grp][1]++;
	    } else {
	      n_never[grp[i]][0]++; n_never[n_grp][0]++;
	    }
	  }
	  /* insample QoI */
	  if (*Insample) { 
	    if (C[i] == 1) { /* compliers */
	      if (*random) {
		dtemp = rlnorm(meano1[i]+gamma1[0]+xiO1[grp[i]][0], sqrt(*sig2)) * 
		  ((meano[i]+gamma[0]+xiO[grp[i]][0]+norm_rand()) > 0);
		dtemp1 = rlnorm(meano1[i]+gamma1[1]+xiO1[grp[i]][0], sqrt(*sig2)) * 
		  ((meano[i]+gamma[1]+xiO[grp[i]][0]+norm_rand()) > 0);
	      } else { 
		dtemp = rlnorm(meano1[i]+gamma1[0], sqrt(*sig2)) * 
		  ((meano[i]+gamma[0]+norm_rand()) > 0);
		dtemp1 = rlnorm(meano1[i]+gamma1[1], sqrt(*sig2)) * 
		  ((meano[i]+gamma[1]+norm_rand()) > 0);
	      }
	      if (R[i] == 1) {
		if (Z[i] == 1) { 
		  dtemp = Y1[i]; 
		} else { 
		  dtemp1 = Y1[i];
		}
	      }
	      Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	      Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    } else if (A[i] == 1) { /* always-takers */ 
	      if (R[i] == 1)
		dtemp = Y1[i];
	      else if (*random) 
		dtemp = rlnorm(meano1[i]+gamma1[2]+xiO1[grp[i]][1], sqrt(*sig2)) * 
		  ((meano[i]+gamma[2]+xiO[grp[i]][1]+norm_rand()) > 0);
	      else
		dtemp = rlnorm(meano1[i]+gamma1[2], sqrt(*sig2)) * 
		  ((meano[i]+gamma[2]+norm_rand()) > 0);
	      YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    } else { /* never-takers */
	      if (R[i] == 1)
		dtemp = Y1[i];
	      else
		dtemp = rlnorm(meano1[i], sqrt(*sig2)) *
		  ((meano[i]+norm_rand()) > 0);
	      YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	    } 
	  } else { /* population QoI */
	    /* compliers */
	    if (*random) {
	      dtemp = exp(meano1[i]+gamma1[0]+xiO1[grp[i]][0]+0.5*sig2[0]) *
		pnorm(meano[i]+gamma[0]+xiO[grp[i]][0], 0, 1, 1, 0);
	      dtemp1 = exp(meano1[i]+gamma1[1]+xiO1[grp[i]][0]+0.5*sig2[0]) *  
		pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 1, 0);
	    } else {
	      dtemp = exp(meano1[i]+gamma1[0]+0.5*sig2[0]) *
		pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
	      dtemp1 = exp(meano1[i]+gamma1[1]+0.5*sig2[0]) *
		pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	    }
	    Y1barC[grp[i]] += dtemp; Y0barC[grp[i]] += dtemp1;
	    Y1barC[n_grp] += dtemp; Y0barC[n_grp] += dtemp1;
	    ITT[grp[i]] += (dtemp-dtemp1)*qC[i];
	    ITT[n_grp] += (dtemp-dtemp1)*qC[i];
	    /* always-takers */
	    if (*random)
	      dtemp = exp(meano1[i]+gamma1[2]+xiO1[grp[i]][1]+0.5*sig2[0]) *
		pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0);
	    else
	      dtemp = exp(meano1[i]+gamma1[2]+0.5*sig2[0]) *
		pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	    YbarA[grp[i]] += dtemp; YbarA[n_grp] += dtemp;
	    /* never-takers */
	    dtemp = exp(meano1[i]+0.5*sig2[0]);
	    YbarN[grp[i]] += dtemp; YbarN[n_grp] += dtemp;
	  }
	}
	
	uniQoIcalMixed(*Insample, n_grp, ITT, Y1barC, Y0barC, n_samp, n_comp,
		       n_never, n_always, p_comp, p_never, CACE, YbarN,
		       YbarA, AT);
	
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_comp[j]; 	  
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = p_never[j];
	for (j = 0; j < (n_grp+1); j++) 
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (n_grp+1); j++)
	  QoI[itempQ++] = YbarN[j];
	if (AT)
	  for (j = 0; j < (n_grp+1); j++)
	    QoI[itempQ++] = YbarA[j];
	
	if (*param) {
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  ssig2[itempS++] = sig2[0];
	  if (AT) {
	    if (*logitC) {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    } else {
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	    }
	  }
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  for (j = 0; j < n_fixedO; j++)
	    coefO1[itempO1++] = gamma1[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
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
	  Rprintf("  Current Acceptance Ratio for fixed effects:");
	  if (AT)
	    for (j = 0; j < n_fixedC*2; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  else
	    for (j = 0; j < n_fixedC; j++)
	      Rprintf("%10g", (double)acc_fixed[j]/(double)main_loop);
	  Rprintf("\n");
	  Rprintf("  Current Acceptance Ratio for random effects:");
	  if (AT)
	    for (j = 0; j < 2; j++)
	      Rprintf("%10g", (double)acc_random[j]/(double)main_loop);
	  else
	    Rprintf("%10g", (double)acc_random[0]/(double)main_loop);
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
  free(gamma1);
  free(grp_obs);
  free(grp_obs1);
  free(Yobs);
  free(Yobs1);
  FreeMatrix(Xc, n_samp+n_fixedC);
  Free3DMatrix(Zc, n_grp, *max_samp_grp + n_randomC);
  FreeMatrix(Xo, n_samp+n_fixedO);
  Free3DMatrix(Zo, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xobs, n_obs+n_fixedO);
  FreeMatrix(Xobs1, n_samp1+n_fixedO);
  Free3DMatrix(Zobs, n_grp, *max_samp_grp + n_randomO);
  Free3DMatrix(Zobs1, n_grp, *max_samp_grp + n_randomO);
  FreeMatrix(Xr, n_samp+n_fixedR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_randomR);
  Free3DMatrix(xiC, 2, n_grp);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiO1, n_grp);
  FreeMatrix(xiR, n_grp);
  Free3DMatrix(Psi, 2, n_randomC);
  FreeMatrix(PsiO, n_randomO);
  FreeMatrix(PsiO1, n_randomO);
  FreeMatrix(PsiR, n_randomR);
  free(meano);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  free(ITT);
  free(CACE);
  free(Y1barC);
  free(Y0barC);
  free(YbarN);
  free(YbarA);
  FreeintMatrix(n_comp, n_grp+1);
  FreeintMatrix(n_never, n_grp+1);
  FreeintMatrix(n_always, n_grp+1);
  free(p_comp);
  free(p_never);
  free(vitemp);
  free(acc_fixed);
  free(acc_random);
  FreeMatrix(mtemp, n_fixedO);
  free(vitemp1);

} /* end of LItwopartMixed */
