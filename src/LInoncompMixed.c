/****

     This file contains models for clustered randomized experiments
     (Frangakis, Rubin, and Zhou Biostatistics) with noncompliance
     under the assumptions of latent ignorability of Frangakis and
     Rubin (1999, Biometrika) .

****/

#include <stddef.h>
#include <stdio.h>      
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"


/* 
   Bayesian binary probit mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LIbprobitMixed(int *Y,         /* binary outcome variable */ 
		    int *R,         /* missingness indicator for Y */
		    int *Z,         /* treatment assignment */
		    int *D,         /* treatment status */ 
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

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);
  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);

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
  double *meana = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  int *Atemp = intArray(n_samp);
  double ***Ztemp = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				   n_randomC +1);
  int *grp_temp = intArray(n_samp);

  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);
  int progress = 1;
  int keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, k, l, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv;
  double dtemp;
  double **mtempC = doubleMatrix(n_fixedC, n_fixedC); 
  double **mtempO = doubleMatrix(n_fixedO, n_fixedO); 
  double **mtempR = doubleMatrix(n_fixedR, n_fixedR); 

  /*** get random seed **/
  GetRNGstate();

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
      Yobs[itemp] = Y[i];
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

  /* prior for fixed effects as additional data points */ 
  itemp = 0; 
  if ((*logitC == 1) && (*AT == 1))
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

  if (*logitC != 1) {
    dcholdc(A0C, n_fixedC, mtempC);
    for (i = 0; i < n_fixedC; i++) {
      Xc[n_samp+i][n_fixedC]=0;
      for (j = 0; j < n_fixedC; j++) {
	Xc[n_samp+i][n_fixedC] += mtempC[i][j]*beta0[j];
	Xc[n_samp+i][j] = mtempC[i][j];
      }
    }
  }

  dcholdc(A0O, n_fixedO, mtempO);
  for (i = 0; i < n_fixedO; i++) {
    Xobs[n_obs+i][n_fixedO]=0;
    for (j = 0; j < n_fixedO; j++) {
      Xobs[n_obs+i][n_fixedO] += mtempO[i][j]*gamma0[j];
      Xobs[n_obs+i][j] = mtempO[i][j];
    }
  }
  
  dcholdc(A0R, n_fixedR, mtempR);
  for (i = 0; i < n_fixedR; i++) {
    Xr[n_samp+i][n_fixedR]=0;
    for (j = 0; j < n_fixedR; j++) {
      Xr[n_samp+i][n_fixedR] += mtempR[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtempR[i][j];
    }
  }
  
  /*** starting values for probabilities ***/
  for (i = 0; i < n_samp; i++) {
    pC[i] = unif_rand(); 
    pN[i] = unif_rand(); 
    pA[i] = unif_rand();
    if (n_miss > 0) {
      prC[i] = unif_rand(); 
      prN[i] = unif_rand(); 
      prA[i] = unif_rand();
    } else { /* no missing values in Y */
      prC[i] = 1;
      prN[i] = 1;
      prA[i] = 1;
    }
  }

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0) {
      bprobitMixedGibbs(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
			n_fixedR, n_randomR, n_grp, 0, 
			delta0, A0R, tau0s[3], T0R, 1);
 
      /* Compute probabilities of R = Robs */ 
      for (j = 0; j < n_grp; j++) vitemp[j] = 0;
      for (i = 0; i < n_samp; i++) {
	dtemp = 0;
	if (*AT) { /* always-takers */
	  for (j = 3; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  if (*random) 
	    for (j = 2; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  else
	    for (j = 0; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    if (*random) 
	      prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	     else 
	       prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) +
		 (1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    if (*random) {
	      prC[i] = R[i]*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 0, 0);
	      prA[i] = R[i]*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 0, 0);
	    } else {
	      prC[i] = R[i]*pnorm(dtemp+delta[0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[0], 0, 1, 0, 0);
	      prA[i] = R[i]*pnorm(dtemp+delta[2], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[2], 0, 1, 0, 0);
	    }
	  }
	} else { /* no always-takers */
	  for (j = 2; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  if (*random) 
	    for (j = 1; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  else
	    for (j = 0; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if (Z[i] == 0) {
	    if (*random)
	      prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) + 
		(1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	    else
	      prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) + 
		(1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  }
	} 
	vitemp[grp[i]]++;
      }
    }

    /** Step 2: COMPLIANCE MODEL **/
    if (*logitC) 
      if (*AT) 
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
      if (*AT) {
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

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
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
      if (*AT) { /* some always-takers */
	meana[i] = 0;
	for (j = 0; j < n_randomC; j++)
	  meana[i] += Zc[grp[i]][vitemp[grp[i]]][j]*xiC[1][grp[i]][j];
	if (*logitC) { /* if logistic regression is used */
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
	if ((Z[i] == 0) && (D[i] == 0)){
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]);
	  else 
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+qN[i]*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1;
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    }
	  } else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0; 
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0; 
	    }
	  }  
	}
	if ((Z[i] == 1) && (D[i] == 1)){
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	  else
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i]-qN[i])*prA[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][0] = 1; Xr[i][0] = 1; 
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    A[i] = 0; Xo[i][2] = 0; Xr[i][2] = 0; 
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 1; Xobs[itemp][2] = 0; 
	      if (*random) {
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1; 
		Zobs[grp[i]][vitemp1[grp[i]]][1] = 0;
	      }
	    }
	  } else {
	    if (*logitC)
	      C[i] = 2;
	    else
	      C[i] = 0; 
	    A[i] = 1; Xo[i][0] = 0; Xr[i][0] = 0; 
	    Xo[i][2] = 1; Xr[i][2] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0; 
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0; 
	      Zo[grp[i]][vitemp[grp[i]]][1] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][1] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 0; 
	      Xobs[itemp][2] = 1; 
	      if (*random) {
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
		Zobs[grp[i]][vitemp1[grp[i]]][1] = 1;
	      } 
	    }
	  }  
	}
      } else { /* no always-takers */
	if (Z[i] == 0){
	  if (*logitC)
	    qC[i] = 1/(1+exp(-meanc[i]));
	  else
	    qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+(1-qC[i])*pN[i]*prN[i]);
	  else
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i])*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1; 
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    } 
	  } else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0;
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
	    } 
	  }
	}
      }
      if (R[i] == 1) {
	itemp++; vitemp1[grp[i]]++;
      }
      vitemp[grp[i]]++;
    }

    /** Step 4: OUTCOME MODEL **/
    bprobitMixedGibbs(Yobs, Xobs, Zobs, grp_obs, gamma, xiO, PsiO,
		      n_obs, n_fixedO, n_randomO, n_grp, 0,
		      gamma0, A0O, tau0s[2], T0O, 1); 
    
    /** Compute probabilities of Y = 1 **/
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
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    if (*random)
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 0, 0);
	    else
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[1], 0, 1, 0, 0);
	    pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	      (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    if (*random) {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[0]+xiO[grp[i]][0], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[0]+xiO[grp[i]][0], 0, 1, 0, 0);
	      pA[i] = Y[i]*pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 0, 0);
	    } else {
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[0], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[0], 0, 1, 0, 0);
	      pA[i] = Y[i]*pnorm(meano[i]+gamma[2], 0, 1, 1, 0) +
		(1-Y[i])*pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	    }
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
	  if (Z[i] == 0) {
	    if (*random)
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 1, 0) + 
		(1-Y[i])*pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 0, 0);
	    else 
	      pC[i] = Y[i]*pnorm(meano[i]+gamma[1], 0, 1, 1, 0) + 
		(1-Y[i])*pnorm(meano[i]+gamma[1], 0, 1, 0, 0);
	    pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	      (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
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
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i];
	  p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	    if (*Insample) { /* insample QoI */
	      if (Z[i] == 1) {
		if (R[i] == 1)
		  Y1barC += (double)Y[i];
		else if (*random)
		  Y1barC += (double)((meano[i]+gamma[0]+xiO[grp[i]][0]+norm_rand()) > 0);
		else
		  Y1barC += (double)((meano[i]+gamma[0]+norm_rand()) > 0);
	      } else {
		if (R[i] == 1)
		  Y0barC += (double)Y[i];
		else if (*random)
		  Y0barC += (double)((meano[i]+gamma[1]+xiO[grp[i]][0]+norm_rand()) > 0);
		else
		  Y0barC += (double)((meano[i]+gamma[1]+norm_rand()) > 0);
	      }
	    } else { /* population QoI */
	      if (*random) {
		Y1barC += pnorm(meano[i]+gamma[0]+xiO[grp[i]][0], 0, 1, 1, 0);
		Y0barC += pnorm(meano[i]+gamma[1]+xiO[grp[i]][0], 0, 1, 1, 0); 
	      } else {
		Y1barC += pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
		Y0barC += pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	      }
	    }
	  } else { /* Estimates for always-takers and never-takers */
	    if (A[i] == 1) {
	      if (Z[i] == 1)
		n_always[1]++;
	      else
		n_always[0]++;
	      if (*Insample)
		if (R[i] == 1)
		  YbarA += (double)Y[i];
		else if (*random)
		  YbarA += (double)((meano[i]+gamma[2]+xiO[grp[i]][1]+norm_rand()) > 0);
		else
		  YbarA += (double)((meano[i]+gamma[2]+norm_rand()) > 0);
	      else if (*random)
		YbarA += pnorm(meano[i]+gamma[2]+xiO[grp[i]][1], 0, 1, 1, 0);
	      else
		YbarA += pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	    } else {
	      if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	      if (*Insample)
		if (R[i] == 1)
		  YbarN += (double)Y[i];
		else 
		  YbarN += (double)((meano[i]+norm_rand()) > 0);
	      else 
		YbarN += pnorm(meano[i], 0, 1, 1, 0);
	    }
	  }
	}
	if (*Insample) { 
	  ITT = Y1barC/(double)(n_comp[1]+n_never[1]+n_always[1]) -
	    Y0barC/(double)(n_comp[0]+n_never[0]+n_always[0]);
	  Y1barC /= (double)n_comp[1];  
	  Y0barC /= (double)n_comp[0]; 
	  p_comp = (double)(n_comp[0]+n_comp[1])/(double)n_samp;
	  p_never = (double)(n_never[0]+n_never[1])/(double)n_samp;
	} else {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
					   a complier */ 
	  p_never /= (double)n_samp; /* Prob. of being a never-taker */
	}
	CACE = Y1barC-Y0barC;    /* CACE */
	YbarN /= (double)(n_never[0]+n_never[1]);
	if (*AT)
	  YbarA /= (double)(n_always[0]+n_always[1]);
	
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
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  if (*AT)
	    if (*logitC)
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    else
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      }
      else
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
  free(meanc);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  free(meana);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  FreeMatrix(Xtemp, n_samp+n_fixedC);
  free(Atemp);
  Free3DMatrix(Ztemp, n_grp, *max_samp_grp + n_randomC);
  free(grp_temp);
  free(n_comp);
  free(n_never);
  free(n_always);
  free(vitemp);
  free(vitemp1);
  free(acc_fixed);
  free(acc_random);
  FreeMatrix(mtempC, n_fixedC);
  FreeMatrix(mtempO, n_fixedO);
  FreeMatrix(mtempR, n_fixedR); 

} /* end of LIbprobitMixed */



/* 
   Bayesian Normal mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LINormalMixed(double *Y,      /* Gaussian outcome variable */ 
		   int *R,         /* missingness indicator for Y */
		   int *Z,         /* treatment assignment */
		   int *D,         /* treatment status */ 
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

  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);
  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);

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
  double *meana = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  int *Atemp = intArray(n_samp);
  double ***Ztemp = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				   n_randomC +1);
  int *grp_temp = intArray(n_samp);

  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);
  int progress = 1;
  int keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);      /* number of acceptance */
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, k, l, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR;
  int itempAv, itempCv, itempOv, itempRv, itempS;
  double dtemp;
  double **mtempC = doubleMatrix(n_fixedC, n_fixedC); 
  double **mtempO = doubleMatrix(n_fixedO, n_fixedO); 
  double **mtempR = doubleMatrix(n_fixedR, n_fixedR); 

  /*** get random seed **/
  GetRNGstate();

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
      Yobs[itemp] = Y[i];
      Xobs[itemp][n_fixedO] = Y[i];
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

  /* prior for fixed effects as additional data points */ 
  itemp = 0; 
  if ((*logitC == 1) && (*AT == 1))
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

  if (*logitC != 1) {
    dcholdc(A0C, n_fixedC, mtempC);
    for (i = 0; i < n_fixedC; i++) {
      Xc[n_samp+i][n_fixedC]=0;
      for (j = 0; j < n_fixedC; j++) {
	Xc[n_samp+i][n_fixedC] += mtempC[i][j]*beta0[j];
	Xc[n_samp+i][j] = mtempC[i][j];
      }
    }
  }

  dcholdc(A0O, n_fixedO, mtempO);
  for (i = 0; i < n_fixedO; i++) {
    Xobs[n_obs+i][n_fixedO]=0;
    for (j = 0; j < n_fixedO; j++) {
      Xobs[n_obs+i][n_fixedO] += mtempO[i][j]*gamma0[j];
      Xobs[n_obs+i][j] = mtempO[i][j];
    }
  }
  
  dcholdc(A0R, n_fixedR, mtempR);
  for (i = 0; i < n_fixedR; i++) {
    Xr[n_samp+i][n_fixedR]=0;
    for (j = 0; j < n_fixedR; j++) {
      Xr[n_samp+i][n_fixedR] += mtempR[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtempR[i][j];
    }
  }
  
  /*** starting values for probabilities ***/
  for (i = 0; i < n_samp; i++) {
    if (R[i] == 1) {
      pC[i] = unif_rand(); 
      pN[i] = unif_rand(); 
      pA[i] = unif_rand();
    } else {
      pC[i] = 1;
      pN[i] = 1;
      pA[i] = 1;
    }
    if (n_miss > 0) {
      prC[i] = unif_rand(); 
      prN[i] = unif_rand(); 
      prA[i] = unif_rand();
    } else { /* no missing values in Y */
      prC[i] = 1;
      prN[i] = 1;
      prA[i] = 1;
    }
  }

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0; itempS = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){
    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0) {
      bprobitMixedGibbs(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
			n_fixedR, n_randomR, n_grp, 0, 
			delta0, A0R, tau0s[3], T0R, 1);
 
      /* Compute probabilities of R = Robs */ 
      for (j = 0; j < n_grp; j++) vitemp[j] = 0;
      for (i = 0; i < n_samp; i++) {
	dtemp = 0;
	if (*AT) { /* always-takers */
	  for (j = 3; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  for (j = 2; j < n_randomR; j++)
	    dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    prC[i] = R[i]*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 0, 0);
	    prA[i] = R[i]*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 0, 0);
	  }
	} else { /* no always-takers */
	  for (j = 2; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  for (j = 1; j < n_randomR; j++)
	    dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if (Z[i] == 0) {
	    prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) + 
	      (1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  }
	} 
	vitemp[grp[i]]++;
      }
    }

    /** Step 2: COMPLIANCE MODEL **/
    if (*logitC) 
      if (*AT) 
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
			n_fixedC, n_randomC, n_grp, 0, 
			beta0, A0C, tau0s[0], T0C, 1);
      if (*AT) {
	/* never-taker vs. always-taker */
	/* subset the data */
	itemp = 0;
	for (j = 0; j < n_grp; j++) {
	  vitemp[j] = 0; vitemp1[j] = 0;
	}
	for (i = 0; i < n_samp; i++) {
	  if (C[i] != 1) {
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

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
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
      if (*AT) { /* some always-takers */
	meana[i] = 0;
	for (j = 0; j < n_randomC; j++)
	  meana[i] += Zc[grp[i]][vitemp[grp[i]]][j]*xiC[1][grp[i]][j];
	if (*logitC) { /* if logistic regression is used */
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
	if ((Z[i] == 0) && (D[i] == 0)){
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1; 
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1; Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    }
	  } else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0; 
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0; Zobs[grp[i]][vitemp1[grp[i]]][0] = 0; 
	    }
	  }  
	}
	if ((Z[i] == 1) && (D[i] == 1)){
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][0] = 1; Xr[i][0] = 1; 
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    A[i] = 0; Xo[i][2] = 0; Xr[i][2] = 0; 
	    Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 1; Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	      Xobs[itemp][2] = 0; Zobs[grp[i]][vitemp1[grp[i]]][1] = 0;
	    }
	  }
	  else {
	    if (*logitC)
	      C[i] = 2;
	    else
	      C[i] = 0; 
	    A[i] = 1; Xo[i][0] = 0; Xr[i][0] = 0; 
	    Xo[i][2] = 1; Xr[i][2] = 1; 
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0; 
	    Zo[grp[i]][vitemp[grp[i]]][1] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][1] = 1;
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 0; Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
	      Xobs[itemp][2] = 1; Zobs[grp[i]][vitemp1[grp[i]]][1] = 1;
	    } 
	  }  
	}
      } else { /* no always-takers */
	if (Z[i] == 0){
	  if (*logitC)
	    qC[i] = 1/(1+exp(-meanc[i]));
	  else
	    qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	  dtemp = qC[i]*pC[i]*prC[i] / 
	    (qC[i]*pC[i]*prC[i]+(1-qC[i])*pN[i]*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	    Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1; Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    } 
	  }
	  else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	    Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	    Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0; Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
	    } 
	  }
	}
      }
      if (R[i] == 1) {
	itemp++; vitemp1[grp[i]]++;
      }
      vitemp[grp[i]]++;
    }

    /** Step 4: OUTCOME MODEL **/
    bNormalMixedGibbs(Yobs, Xobs, Zobs, grp_obs, gamma, 
		      xiO, sig2, PsiO, n_obs, n_fixedO, 
		      n_randomO, n_grp, 0, gamma0, A0O, 
		      0, *nu0, *s0, tau0s[2], T0O, 1); 
    
    /** Compute probabilities of Y = 1 **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      meano[i] = 0;
      if (*AT) { /* always-takers */
	for (j = 3; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	for (j = 2; j < n_randomO; j++)
	  meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1) {
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    pC[i] = dnorm(Y[i], meano[i]+gamma[1]+xiO[grp[i]][0], sqrt(*sig2), 0);
	    pN[i] = dnorm(Y[i], meano[i], sqrt(*sig2), 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    pC[i] = dnorm(Y[i], meano[i]+gamma[0]+xiO[grp[i]][0], sqrt(*sig2), 0);
	    pA[i] = dnorm(Y[i], meano[i]+gamma[2]+xiO[grp[i]][1], sqrt(*sig2), 0);
	  }
	}
      } else { /* no always-takers */
	for (j = 2; j < n_fixedO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	for (j = 1; j < n_randomO; j++)
	  meano[i] += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	if (R[i] == 1)
	  if (Z[i] == 0) {
	    pC[i] = dnorm(Y[i], meano[i]+gamma[1]+xiO[grp[i]][0], sqrt(*sig2), 0);
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
	p_comp = 0; p_never = 0; 
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i];
	  p_never += qN[i];
	  if (C[i] == 1) { /* ITT effects */
	    if (Z[i] == 1) 
	      n_comp[1]++;
	    else
	      n_comp[0]++;
	    if (*Insample) { /* insample QoI */
	      if (Z[i] == 1) 
		if (R[i] == 1)
		  Y1barC += Y[i];
		else 
		  Y1barC += rnorm(meano[i]+gamma[0]+xiO[grp[i]][0], sqrt(*sig2));
	      else 
		if (R[i] == 1)
		  Y0barC += Y[i];
		else
		  Y0barC += rnorm(meano[i]+gamma[1]+xiO[grp[i]][0], sqrt(*sig2));
	    } else { /* population QoI */
	      Y1barC += (meano[i]+gamma[0]+xiO[grp[i]][0]);
	      Y0barC += (meano[i]+gamma[1]+xiO[grp[i]][0]); 
	    }
	  } else { /* Estimates for always-takers and never-takers */
	    if (A[i] == 1) {
	      if (Z[i] == 1)
		n_always[1]++;
	      else
		n_always[0]++;
	      if (*Insample)
		if (R[i] == 1)
		  YbarA += Y[i];
		else 
		  YbarA += rnorm(meano[i]+gamma[2]+xiO[grp[i]][1], sqrt(*sig2));
	      else 
		YbarA += (meano[i]+gamma[2]+xiO[grp[i]][1]);
	    } else {
	      if (Z[i] == 1)
		n_never[1]++;
	      else
		n_never[0]++;
	      if (*Insample)
		if (R[i] == 1)
		  YbarN += Y[i];
		else 
		  YbarN += rnorm(meano[i], sqrt(*sig2));
	      else 
		YbarN += meano[i];
	    }
	  }
	}
	if (*Insample) { 
	  ITT = Y1barC/(double)(n_comp[1]+n_never[1]+n_always[1]) -
	    Y0barC/(double)(n_comp[0]+n_never[0]+n_always[0]);
	  Y1barC /= (double)n_comp[1];  
	  Y0barC /= (double)n_comp[0]; 
	  p_comp = (double)(n_comp[0]+n_comp[1])/(double)n_samp;
	  p_never = (double)(n_never[0]+n_never[1])/(double)n_samp;
	} else {
	  ITT = (Y1barC-Y0barC)/(double)n_samp;     /* ITT effect */
	  Y1barC /= (double)(n_comp[0]+n_comp[1]);  
	  Y0barC /= (double)(n_comp[0]+n_comp[1]); 
	  p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
					   a complier */ 
	  p_never /= (double)n_samp; /* Prob. of being a never-taker */
	}
	CACE = Y1barC-Y0barC;    /* CACE */
	YbarN /= (double)(n_never[0]+n_never[1]);
	if (*AT)
	  YbarA /= (double)(n_always[0]+n_always[1]);
	
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
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  ssig2[itempS++] = sig2[0];
	  if (*AT)
	    if (*logitC)
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    else
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      }
      else
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
  free(meanc);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  free(meana);
  free(n_comp);
  free(n_never);
  free(n_always);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  FreeMatrix(Xtemp, n_samp+n_fixedC);
  free(Atemp);
  Free3DMatrix(Ztemp, n_grp, *max_samp_grp + n_randomC);
  free(grp_temp);
  free(vitemp);
  free(vitemp1);
  free(acc_fixed);
  free(acc_random);
  FreeMatrix(mtempC, n_fixedC);
  FreeMatrix(mtempO, n_fixedO);
  FreeMatrix(mtempR, n_fixedR); 

} /* end of LINormalMixed */


/* 
   Bayesian binary probit mixed effects model for clustered randomized
   experiments with noncompliance; Latent Ignorability assumption for
   subsequent missing outcomes 
*/

void LIboprobitMixed(int *Y,         /* binary outcome variable */ 
		     int *R,         /* missingness indicator for Y */
		     int *Z,         /* treatment assignment */
		     int *D,         /* treatment status */ 
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
  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);

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
  double *meana = doubleArray(n_samp);

  /* prior precision matrices */
  double **A0C = doubleMatrix(n_fixedC*2, n_fixedC*2);
  double **A0O = doubleMatrix(n_fixedO, n_fixedO);
  double **A0R = doubleMatrix(n_fixedR, n_fixedR);
  
  /* prior scale for Psi's */
  double **T0C = doubleMatrix(n_randomC, n_randomC);
  double **T0A = doubleMatrix(n_randomC, n_randomC);
  double **T0O = doubleMatrix(n_randomO, n_randomO);
  double **T0R = doubleMatrix(n_randomR, n_randomR);

  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_fixedC, n_fixedC+1);
  int *Atemp = intArray(n_samp);
  double ***Ztemp = doubleMatrix3D(n_grp, *max_samp_grp + n_randomC,
				   n_randomC +1);
  int *grp_temp = intArray(n_samp);

  /* quantities of interest: ITT, CACE  */
  double *ITT = doubleArray(n_cat-1);
  double *CACE = doubleArray(n_cat-1);
  double *Y1barC = doubleArray(n_cat-1); 
  double *Y0barC = doubleArray(n_cat-1); 
  double *YbarN = doubleArray(n_cat-1);
  double *YbarA = doubleArray(n_cat-1);
  int *n_comp = intArray(2);          /* number of compliers */
  int *n_never = intArray(2);
  int *n_always = intArray(2);
  double p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int *vitemp = intArray(n_grp);
  int *vitemp1 = intArray(n_grp);
  int progress = 1;
  int keep = 1;
  int *acc_fixed = intArray(n_fixedC*2);    /* number of acceptance */
  int *acc_tau = intArray(0);
  int *acc_random = intArray(2*n_grp);      /* number of acceptance */
  int i, j, k, l, main_loop;
  int itempP = ftrunc((double) *n_gen/10);
  int itemp, itempA, itempC, itempO, itempQ, itempR, itempT;
  int itempAv, itempCv, itempOv, itempRv;
  double dtemp, dtemp1;
  double **mtempC = doubleMatrix(n_fixedC, n_fixedC); 
  double **mtempO = doubleMatrix(n_fixedO, n_fixedO); 
  double **mtempR = doubleMatrix(n_fixedR, n_fixedR); 

  /*** get random seed **/
  GetRNGstate();

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
      Yobs[itemp] = Y[i];
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

  /* prior for fixed effects as additional data points */ 
  itemp = 0; 
  if (*logitC && *AT)
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

  if (*logitC != 1) {
    dcholdc(A0C, n_fixedC, mtempC);
    for (i = 0; i < n_fixedC; i++) {
      Xc[n_samp+i][n_fixedC]=0;
      for (j = 0; j < n_fixedC; j++) {
	Xc[n_samp+i][n_fixedC] += mtempC[i][j]*beta0[j];
	Xc[n_samp+i][j] = mtempC[i][j];
      }
    }
  }

  dcholdc(A0O, n_fixedO, mtempO);
  for (i = 0; i < n_fixedO; i++) {
    Xobs[n_obs+i][n_fixedO]=0;
    for (j = 0; j < n_fixedO; j++) {
      Xobs[n_obs+i][n_fixedO] += mtempO[i][j]*gamma0[j];
      Xobs[n_obs+i][j] = mtempO[i][j];
    }
  }
  
  dcholdc(A0R, n_fixedR, mtempR);
  for (i = 0; i < n_fixedR; i++) {
    Xr[n_samp+i][n_fixedR]=0;
    for (j = 0; j < n_fixedR; j++) {
      Xr[n_samp+i][n_fixedR] += mtempR[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtempR[i][j];
    }
  }
  
  /*** starting values for probabilities ***/
  for (i = 0; i < n_samp; i++) {
    pC[i] = unif_rand(); 
    pN[i] = unif_rand(); 
    pA[i] = unif_rand();
    if (n_miss > 0) {
      prC[i] = unif_rand(); 
      prN[i] = unif_rand(); 
      prA[i] = unif_rand();
    } else { /* no missing values in Y */
      prC[i] = 1;
      prN[i] = 1;
      prA[i] = 1;
    }
  }

  /*** Gibbs Sampler! ***/
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0; itempR = 0; itempT = 0;   
  itempAv = 0; itempCv = 0; itempOv = 0; itempRv = 0;   
  for (j = 0; j < n_fixedC*2; j++)
    acc_fixed[j] = 0;
  acc_random[0] = 0; acc_random[1] = 0; acc_tau[0] = 0;
  for (main_loop = 1; main_loop <= *n_gen; main_loop++){

    /* Step 1: RESPONSE MODEL */
    if (n_miss > 0) {
      bprobitMixedGibbs(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp, 
			n_fixedR, n_randomR, n_grp, 0, 
			delta0, A0R, tau0s[3], T0R, 1);
      
      /* Compute probabilities of R = Robs */ 
      for (j = 0; j < n_grp; j++) vitemp[j] = 0;
      for (i = 0; i < n_samp; i++) {
	dtemp = 0;
	if (*AT) { /* always-takers */
	  for (j = 3; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  if (*random)
	    for (j = 2; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  else 
	    for (j = 0; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    if (*random)
	      prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	    else
	      prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    if (*random) {
	      prC[i] = R[i]*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[0]+xiR[grp[i]][0], 0, 1, 0, 0);
	      prA[i] = R[i]*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[2]+xiR[grp[i]][1], 0, 1, 0, 0);
	    } else {
	      prC[i] = R[i]*pnorm(dtemp+delta[0], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[0], 0, 1, 0, 0);
	      prA[i] = R[i]*pnorm(dtemp+delta[2], 0, 1, 1, 0) +
		(1-R[i])*pnorm(dtemp+delta[2], 0, 1, 0, 0);
	    }
	  }
	} else { /* no always-takers */
	  for (j = 2; j < n_fixedR; j++)
	    dtemp += Xr[i][j]*delta[j];
	  if (*random)
	    for (j = 1; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  else
	    for (j = 0; j < n_randomR; j++)
	      dtemp += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  if (Z[i] == 0) {
	    if (*random) 
	      prC[i] = R[i]*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 1, 0) + 
		(1-R[i])*pnorm(dtemp+delta[1]+xiR[grp[i]][0], 0, 1, 0, 0);
	    else 
	      prC[i] = R[i]*pnorm(dtemp+delta[1], 0, 1, 1, 0) + 
		(1-R[i])*pnorm(dtemp+delta[1], 0, 1, 0, 0);
	    prN[i] = R[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-R[i])*pnorm(dtemp, 0, 1, 0, 0);
	  }
	} 
	vitemp[grp[i]]++;
      }
    }
    
    /** Step 2: COMPLIANCE MODEL **/
    if (*logitC) 
      if (*AT) 
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
      if (*AT) {
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

    /* Step 3: SAMPLE COMPLIANCE COVARITE */
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
      if (*AT) { /* some always-takers */
	meana[i] = 0;
	for (j = 0; j < n_randomC; j++)
	  meana[i] += Zc[grp[i]][vitemp[grp[i]]][j]*xiC[1][grp[i]][j];
	if (*logitC) { /* if logistic regression is used */
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
	if ((Z[i] == 0) && (D[i] == 0)){
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+qN[i]*pN[i]*prN[i]);
	  else 
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+qN[i]*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1;
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    }
	  } else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0; 
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0; 
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0; 
	    }
	  }  
	}
	if ((Z[i] == 1) && (D[i] == 1)){
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+(1-qC[i]-qN[i])*pA[i]*prA[i]);
	  else
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i]-qN[i])*prA[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][0] = 1; Xr[i][0] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    A[i] = 0; Xo[i][2] = 0; Xr[i][2] = 0;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][1] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][1] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 1; Xobs[itemp][2] = 0; 
	      if (*random) {
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
		Zobs[grp[i]][vitemp1[grp[i]]][1] = 0;
	      }
	    }
	  } else {
	    if (*logitC)
	      C[i] = 2;
	    else
	      C[i] = 0; 
	    A[i] = 1; Xo[i][0] = 0; Xr[i][0] = 0; 
	    Xo[i][2] = 1; Xr[i][2] = 1; 
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0; 
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0; 
	      Zo[grp[i]][vitemp[grp[i]]][1] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][1] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][0] = 0; 
	      Xobs[itemp][2] = 1; 
	      if (*random) {
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
		Zobs[grp[i]][vitemp1[grp[i]]][1] = 1;
	      } 
	    }
	  }  
	}
      } else { /* no always-takers */
	if (Z[i] == 0){
	  if (*logitC)
	    qC[i] = 1/(1+exp(-meanc[i]));
	  else
	    qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	  if (R[i] == 1)
	    dtemp = qC[i]*pC[i]*prC[i] / 
	      (qC[i]*pC[i]*prC[i]+(1-qC[i])*pN[i]*prN[i]);
	  else
	    dtemp = qC[i]*prC[i]/(qC[i]*prC[i]+(1-qC[i])*prN[i]);
	  if (unif_rand() < dtemp) {
	    C[i] = 1; Xo[i][1] = 1; Xr[i][1] = 1;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 1;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 1;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 1;
	      if (*random) 
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 1;
	    } 
	  } else {
	    C[i] = 0; Xo[i][1] = 0; Xr[i][1] = 0;
	    if (*random) {
	      Zo[grp[i]][vitemp[grp[i]]][0] = 0;
	      Zr[grp[i]][vitemp[grp[i]]][0] = 0;
	    }
	    if (R[i] == 1) {
	      Xobs[itemp][1] = 0; 
	      if (*random)
		Zobs[grp[i]][vitemp1[grp[i]]][0] = 0;
	    } 
	  }
	}
      }
      if (R[i] == 1) {
	itemp++; vitemp1[grp[i]]++;
      }
      vitemp[grp[i]]++;
    }

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
    
    /** Compute probabilities of Y = j **/
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
	  if ((Z[i] == 0) && (D[i] == 0)) {
	    if (*random)
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0);
	    else
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1], 1, 1, 0);
	    if (Y[i] == 0)
	      pN[i] = pnorm(tau[0], meano[i], 1, 1, 0);
	    else
	      pN[i] = pnorm(tau[Y[i]], meano[i], 1, 1, 0) -
		pnorm(tau[Y[i]-1], meano[i], 1, 1, 0);
	  } 
	  if ((Z[i] == 1) && (D[i] == 1)) {
	    if (*random) {
	      if (Y[i] == 0) {
		pC[i] = pnorm(tau[0], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0); 
		pA[i] = pnorm(tau[0], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0);
	      } else {
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0) 
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0); 
		pA[i] = pnorm(tau[Y[i]], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0);
	      }
	    } else {
	      if (Y[i] == 0) {
		pC[i] = pnorm(tau[0], meano[i]+gamma[0], 1, 1, 0); 
		pA[i] = pnorm(tau[0], meano[i]+gamma[2], 1, 1, 0);
	      } else {
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[0], 1, 1, 0) 
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[0], 1, 1, 0); 
		pA[i] = pnorm(tau[Y[i]], meano[i]+gamma[2], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[2], 1, 1, 0);
	      }
	    }
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
	  if (Z[i] == 0) {
	    if (*random)
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0);
	    else
	      if (Y[i] == 0)
		pC[i] = pnorm(tau[0], meano[i]+gamma[1], 1, 1, 0);
	      else
		pC[i] = pnorm(tau[Y[i]], meano[i]+gamma[1], 1, 1, 0)
		  - pnorm(tau[Y[i]-1], meano[i]+gamma[1], 1, 1, 0);
	    if (Y[i] == 0)
	      pN[i] = pnorm(tau[0], meano[i], 1, 1, 0);
	    else
	      pN[i] = pnorm(tau[Y[i]], meano[i], 1, 1, 0)
		- pnorm(tau[Y[i]-1], meano[i], 1, 1, 0); 
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
	for (j = 0; j < (n_cat-1); j++) {
	  Y1barC[j] = 0; Y0barC[j] = 0; YbarN[j] = 0; YbarA[j] = 0;
	}
	for (i = 0; i < n_samp; i++){
	  p_comp += qC[i];
	  p_never += qN[i];
	  if (Z[i] == 1) 
	    if (C[i] == 1)
	      n_comp[1]++;
	    else if (A[i] == 1)
	      n_always[1]++;
	    else
	      n_never[1]++;
	  else
	    if (C[i] == 1)
	      n_comp[0]++;
	    else if (A[i] == 1)
	      n_always[0]++;
	    else 
	      n_never[0]++;
	  for (j = 1; j < n_cat; j++) {
	    if (C[i] == 1) { /* ITT effects */
	      if (*Insample) { /* insample QoI */
		if (Z[i] == 1) {
		  if (R[i] == 1)
		    Y1barC[j-1] += (double)(Y[i] == j);
		  else 
		    if (*random)
		      Y1barC[j-1] += (double)(unif_rand() < 
					      (pnorm(tau[j], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0) 
					       -pnorm(tau[j-1], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0)));
		    else
		      Y1barC[j-1] += (double)(unif_rand() < 
					      (pnorm(tau[j], meano[i]+gamma[0], 1, 1, 0) 
					       -pnorm(tau[j-1], meano[i]+gamma[0], 1, 1, 0)));
		} else {
		  if (R[i] == 1)
		    Y0barC[j-1] += (double)(Y[i] == j);
		  else
		    if (*random)
		      Y0barC[j-1] += (double)(unif_rand() < 
					      (pnorm(tau[j], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0) 
					       -pnorm(tau[j-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0)));
		    else
		      Y0barC[j-1] += (double)(unif_rand() < 
					      (pnorm(tau[j], meano[i]+gamma[1], 1, 1, 0) 
					       -pnorm(tau[j-1], meano[i]+gamma[1], 1, 1, 0)));		
		}
	      } else { /* population QoI */
		if (Z[i] == 1)
		  if (*random)
		    Y1barC[j-1] += (pnorm(tau[j], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[0]+xiO[grp[i]][0], 1, 1, 0));
		  else
		    Y1barC[j-1] += (pnorm(tau[j], meano[i]+gamma[0], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[0], 1, 1, 0));
		else
		  if (*random)
		    Y0barC[j-1] += (pnorm(tau[j], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[1]+xiO[grp[i]][0], 1, 1, 0));
		  else
		    Y0barC[j-1] += (pnorm(tau[j], meano[i]+gamma[1], 1, 1, 0) 
				    -pnorm(tau[j-1], meano[i]+gamma[1], 1, 1, 0));
	      }
	    } else { /* Estimates for always-takers and never-takers */
	      if (A[i] == 1) {
		if (*Insample)
		  if (R[i] == 1)
		    YbarA[j-1] += (double)(Y[i] == j);
		  else
		    if (*random)
		      YbarA[j-1] += (double)(unif_rand() < 
					     (pnorm(tau[j], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0) 
					      -pnorm(tau[j-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0)));	
		    else
		      YbarA[j-1] += (double)(unif_rand() < 
					     (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
					      -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0)));	
		else 
		  if (*random)
		    YbarA[j-1] += (pnorm(tau[j], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0) 
				   -pnorm(tau[j-1], meano[i]+gamma[2]+xiO[grp[i]][1], 1, 1, 0));
		  else
		    YbarA[j-1] += (pnorm(tau[j], meano[i]+gamma[2], 1, 1, 0) 
				   -pnorm(tau[j-1], meano[i]+gamma[2], 1, 1, 0));
	      } else {
		if (*Insample)
		  if (R[i] == 1)
		    YbarN[j-1] += (double)(Y[i] == j);
		  else 
		    YbarN[j-1] += (double)(unif_rand() < 
					   (pnorm(tau[j], meano[i], 1, 1, 0) 
					    -pnorm(tau[j-1], meano[i], 1, 1, 0)));		
		else 
		  YbarN[j-1] += (pnorm(tau[j], meano[i], 1, 1, 0) 
				 -pnorm(tau[j-1], meano[i], 1, 1, 0));
	      }
	    }
	  }
	}
	
	if (*Insample) { 
	  for (j = 0; j < (n_cat-1); j++) {
	    ITT[j] = Y1barC[j]/(double)(n_comp[1]+n_never[1]+n_always[1]) -
	      Y0barC[j]/(double)(n_comp[0]+n_never[0]+n_always[0]);
	    Y1barC[j] /= (double)n_comp[1];  
	    Y0barC[j] /= (double)n_comp[0]; 
	  }
	  p_comp = (double)(n_comp[0]+n_comp[1])/(double)n_samp;
	  p_never = (double)(n_never[0]+n_never[1])/(double)n_samp;
	} else {
	  for (j = 0; j < (n_cat-1); j++) {
	    ITT[j] = (Y1barC[j]-Y0barC[j])/(double)n_samp;     /* ITT effect */
	    Y1barC[j] /= (double)(n_comp[0]+n_comp[1]);  
	    Y0barC[j] /= (double)(n_comp[0]+n_comp[1]);
	  } 
	  p_comp /= (double)n_samp;  /* ITT effect on D; Prob. of being
					a complier */ 
	  p_never /= (double)n_samp; /* Prob. of being a never-taker */
	}
	for (j = 0; j < (n_cat-1); j++) {
	  CACE[j] = Y1barC[j]-Y0barC[j];    /* CACE */
	  YbarN[j] /= (double)(n_never[0]+n_never[1]);
	  if (*AT)
	    YbarA[j] /= (double)(n_always[0]+n_always[1]);
	}
	
	for (j = 0; j < (n_cat-1); j++) 
	  QoI[itempQ++] = ITT[j];   
	for (j = 0; j < (n_cat-1); j++) 
	  QoI[itempQ++] = CACE[j];   
	for (j = 0; j < (n_cat-1); j++) 
	  QoI[itempQ++] = Y1barC[j];
	for (j = 0; j < (n_cat-1); j++) 
	  QoI[itempQ++] = Y0barC[j];
	for (j = 0; j < (n_cat-1); j++) 
	  QoI[itempQ++] = YbarN[j];
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	if (*AT)
	  for (j = 0; j < (n_cat-1); j++) 
	    QoI[itempQ++] = YbarA[j];
	
	if (*param) {
	  for (j = 0; j < (n_cat-1); j++)
	    tauO[itempT++] = tau[j];
	  for (j = 0; j < n_fixedC; j++)
	    coefC[itempC++] = betaC[j];
	  if (*AT)
	    if (*logitC)
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaC[j+n_fixedC];
	    else
	      for (j = 0; j < n_fixedC; j++)
		coefA[itempA++] = betaA[j];
	  for (j = 0; j < n_fixedO; j++)
	    coefO[itempO++] = gamma[j];
	  if (n_miss > 0) 
	    for (j = 0; j < n_fixedR; j++)
	      coefR[itempR++] = delta[j];
	}
	keep = 1;
      }
      else
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
  free(meanc);
  free(pC);
  free(pN);
  free(prC);
  free(prN);
  free(prA);
  free(qC);
  free(qN);
  free(pA);
  free(meana);
  free(n_comp);
  free(n_never);
  free(n_always);
  FreeMatrix(A0C, n_fixedC*2);
  FreeMatrix(A0O, n_fixedO);
  FreeMatrix(A0R, n_fixedR);
  FreeMatrix(T0C, n_randomC);
  FreeMatrix(T0A, n_randomC);
  FreeMatrix(T0O, n_randomO);
  FreeMatrix(T0R, n_randomR);
  FreeMatrix(Xtemp, n_samp+n_fixedC);
  free(Atemp);
  Free3DMatrix(Ztemp, n_grp, *max_samp_grp + n_randomC);
  free(grp_temp);
  free(vitemp);
  free(vitemp1);
  free(acc_fixed);
  free(acc_random);
  free(acc_tau);
  FreeMatrix(mtempC, n_fixedC);
  FreeMatrix(mtempO, n_fixedO);
  FreeMatrix(mtempR, n_fixedR); 

} /* end of LIboprobitMixed */
