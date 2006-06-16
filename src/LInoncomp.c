#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

/* Bayesian binary probit for randomized experiments with
   noncompliance; Latent Ignorability assumption for subsequent
   missing outcomes */

void LIbprobit(int *Y,         /* binary outcome variable */ 
	       int *R,         /* missingness indicator for Y */
	       int *Z,         /* treatment assignment */
	       int *D,         /* treatment status */ 
	       int *C,         /* compliance status; complier = 1,
			          noncomplier = 0 */
	       int *A,         /* always-takers; always-taker = 1, others
			          = 0 */
	       int *AT,        /* Are there always-takers? */
	       double *dXc,    /* model matrix for compliance model */
	       double *dXo,    /* model matrix for outcome model */
	       double *betaC,  /* coefficients for compliance model */
	       double *betaA,  /* coefficients for always-takers model */
	       double *gamma,  /* coefficients for outcome model */
	       int *in_samp,   /* number of observations */
	       int *in_gen,    /* number of Gibbs draws */
	       int *in_covC,   /* number of covariates for compliance
				  model */ 
	       int *in_covO,   /* number of covariates for outcome
				  model */
	       double *beta0,  /* prior mean for betaC and betaA */ 
	       double *gamma0, /* prior mean for gamma */
	       double *dA0C,   /* prior precision for betaC and betaA */ 
	       double *dA0O,   /* prior precision for gamma */ 
	       int *param,     /* Want to keep paramters? */
	       int *mda,       /* Want to use marginal data
				  augmentation? */
	       int *iBurnin,   /* number of burnin */
	       int *iKeep,     /* keep ?th draws */
	       int *verbose,   /* print out messages */
	       double *coefC,  /* Storage for coefficients of the
				  compliance model */
	       double *coefA,  /* Storage for coefficients of the
				  always-takers model */
	       double *coefO,  /* Storage for coefficients of the
				  outcome model */
	       double *QoI     /* Storage of quantities of interest */
	       ) {
  /** counters **/
  int n_samp = *in_samp;
  int n_gen = *in_gen;
  int n_covC = *in_covC;
  int n_covO = *in_covO;

  /*** data ***/
  /* covariates for the compliance model */
  double **Xc = doubleMatrix(n_samp+n_covC, n_covC+1);
  /* covariates for the outcome model */     
  double **Xo = doubleMatrix(n_samp+n_covO, n_covO+1);    
  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_covC, n_covC+1);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp); 
  double qC, qA;
  /* probability of being a always-taker */
  if (AT)
    double *pA = doubleArray(n_samp); 
  double **A0C = doubleMatrix(n_covC, n_covC);
  double **A0O = doubleMatrix(n_covO, n_covO);

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, l, main_loop;  
  int itemp, itempP = ftrunc((double) n_gen/10);
  double dtemp, ndraw, cdraw;
  double *vtemp = doubleArray(n_samp);
  double **mtempC = doubleMatrix(n_covC, n_covC); 
  double **mtempO = doubleMatrix(n_covO, n_covO); 

  /*** get random seed **/
  GetRNGstate();

  /*** read the data ***/
  itemp = 0;
  for (j = 0; j < n_covC; j++)
    for (i = 0; i < n_samp; i++)
      Xc[i][j] = dXc[itemp++];

  itemp = 0;
  for (j = 0; j < n_covO; j++)
    for (i = 0; i < n_samp; i++)
      Xo[i][j] = dXo[itemp++];
  
  /*** read the prior and it as additional data points ***/ 
  itemp = 0;
  for (k = 0; k < n_covC; k++)
    for (j = 0; j < n_covC; j++)
      A0C[j][k] = dA0C[itemp++];

  itemp = 0;
  for (k = 0; k < n_covO; k++)
    for (j = 0; j < n_covO; j++)
      A0O[j][k] = dAo[itemp++];

  dcholdc(A0C, n_covC, mtempC);
  for(i = 0; i < n_covC; i++) {
    Xc[n_samp+i][n_covC]=0;
    for(j = 0; j < n_covC; j++) {
      Xc[n_samp+i][n_covC] += mtempC[i][j]*beta0[j];
      Xc[n_samp+i][j] = mtempC[i][j];
    }
  }

  dcholdc(A0O, n_covO, mtempO);
  for(i = 0; i < n_covO; i++) {
    Xo[n_samp+i][n_covO]=0;
    for(j = 0; j < n_covO; j++) {
      Xo[n_samp+i][n_covO] += mtempO[i][j]*gamma0[j];
      Xo[n_samp+i][j] = mtempO[i][j];
    }
  }

  /*** starting values for probabilities ***/
  for (i = 0; i < n_samp; i++) {
    pC[i] = unif_rand(); 
    pN[i] = unif_rand(); 
    if (AT)
      pA[i] = unif_rand();
  }

  /*** Gibbs Sampler! ***/
  itemp=0;     
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /* Sample complier status for control group */
    if (AT) { /* some always-takers */
      for(i = 0; i < n_samp; i++) {
	dtemp = 0;
	for(j = 0; j < n_covC; j++) 
	  dtemp += Xc[i][j]*betaC[j];
	qC = pnorm(dtemp, 0, 1, 1, 0);
	dtemp = 0;
	for(j = 0; j < n_covC; j++) 
	  dtemp += Xc[i][j]*betaA[j];
	qA = (1-qC)*pnorm(dtemp, 0, 1, 1, 0);
	if (Z[i] == 1 & D[i] == 1){
	  if (unif_rand() < (qC*pC[i]/(qC*pC[i]+qA*pA[i]))) { 
	    C[i] = 1; 
	    Xo[i][0] = 1;
	    Xo[i][2] = 0;
	  }
	  else {
	    C[i] = 0; 
	    Xo[i][0] = 0;
	    Xo[i][2] = 1;
	  }  
	}
	if (Z[i] == 0 & D[i] == 0){
	  if (unif_rand() < (qC*pC[i]/(qC*pC[i]+(1-qC-qA)*pN[i]))) { 
	    C[i] = 1; 
	    Xo[i][1] = 1;
	  }
	  else {
	    C[i] = 0; 
	    Xo[i][1] = 0;
	  }  
	}
      }
    } else { /* no always-takers */
      for(i = 0; i < n_samp; i++)
	if(Z[i] == 0){
	  dtemp = 0;
	  for(j = 0; j < n_covC; j++) 
	    dtemp += Xc[i][j]*betaC[j];
	  qC = pnorm(dtemp, 0, 1, 1, 0);
	  if (unif_rand() < (qC*pC[i]/(qC*pC[i]+(1-qC)*pN[i]))) { 
	    C[i] = 1; 
	    Xo[i][1] = 1; 
	  }
	  else {
	    C[i] = 0; 
	    Xo[i][1] = 0;
	  }
	}
    }

    /** COMPLIANCE MODEL **/    
    bprobitGibbs(C, Xc, betaC, n_samp, n_covC, 0, beta0, A0C, *mda, 1);

    /** ALWAYS-TAKERS MODEL **/
    if (AT) {
      /* subset the data */
      itemp = 0;
      for (i = 0; i < n_samp; i++)
	if (C[i] == 0)
	  for (j = 0; j < n_covC; j++)
	    Xtemp[itemp++][j] = X[i][j];

    }
      

    /** OUTCOME MODEL **/
    bprobitGibbs(Y, Xo, gamma, n_samp, n_covO, 0, gamma0, A0O, *mda, 1);

    /** Compute probabilities of Y = 1 **/ 
    if (AT) {
      


    } else { /* no always-takers */
      for(i = 0; i < n_samp; i++){
	if (Z[i] == 0) {
	  dtemp = 0;
	  for(j = 2; j < n_covO; j++)
	    dtemp += Xo[i][j]*gamma[j];
	  if (Y[i] == 0){
	    pC[i] = Y[i]*pnorm(dtemp+gamma[1], 0, 1, 1, 0) + 
	      (1-Y[i])*pnorm(dtemp+gamma[1], 0, 1, 0, 0);
	    pN[i] = Y[i]*pnorm(dtemp, 0, 1, 1, 0) +
	      (1-Y[i])*pnorm(dtemp, 0, 1, 0, 0);
	  }
	}
      } 
    }
    
    /** Compute quantities of interest **/
    n_comp = 0; n_compC = 0; n_ncompC = 0; base[0] = 0; base[1] = 0; 
    for (i = 0; i <= Ymax; i++)
      ITTc[i] = 0;
    if (*smooth) {
      for(i = 0; i < n11; i++){
	if(C[i] == 1) {
	  n_comp++;
	  if (Z[i] == 0) {
	    n_compC++;
	    base[0] += (double)Y[i];
	  }
	  pcmean = vtemp[i];
	  pnmean = vtemp[i]-treat[i]+gamma[0];
	  ndraw = rnorm(pnmean, 1);
	  cdraw = rnorm(pcmean, 1);
	  if (*insample && R[i]==0) 
	    dtemp = (double)(Y[i]==0) - (double)(ndraw < 0);
	  else
	    dtemp = pnorm(0, pcmean, 1, 1, 0) - pnorm(0, pnmean, 1, 1, 0);
	  ITTc[0] += dtemp;
	  if (Ymax == 1) { /* binary probit */
	    if (*insample && R[i]==0) 
	      dtemp = (double)Y[i] - (double)(ndraw > 0);
	    else
	      dtemp = pnorm(0, pcmean, 1, 0, 0) - pnorm(0, pnmean, 1, 0, 0);
	    ITTc[1] += dtemp;
	  }
	  else { /* ordered probit */
	    if (*insample && R[i]==0) 
	      dtemp = (double)(Y[i]==Ymax) - (double)(ndraw > tau[Ymax-1]);
	    else
	      dtemp = pnorm(tau[Ymax-1], pcmean, 1, 0, 0) -
		pnorm(tau[Ymax-1], pnmean, 1, 0, 0);
	    ITTc[Ymax] += dtemp; 
	    for (j = 1; j < Ymax; j++) {
	      if (*insample && R[i]==0)
		  dtemp = (double)(Y[i]==j) - (double)(ndraw < tau[j] &&
						       ndraw > tau[j-1]);
	      else
		dtemp = (pnorm(tau[j], pcmean, 1, 1, 0) - 
			 pnorm(tau[j-1], pcmean, 1, 1, 0)) 
		  - (pnorm(tau[j], pnmean, 1, 1, 0) - 
		     pnorm(tau[j-1], pnmean, 1, 1, 0));
	      ITTc[j] += dtemp;
	    }
	  }
	}
	else
	  if (Z[i] == 0) {
	    n_ncompC++;
	    base[1] += (double)Y[i];
	  } 
      }
    }
    else {
      for(i = 0; i < n_samp; i++){
	if(C[i] == 1) {
	  n_comp++;
	  if (Z[i] == 1) {
	    pcmean = vtemp[i];
	    pnmean = vtemp[i]-gamma[0]+gamma[1];
	  }
	  else {
	    n_compC++;
	    base[0] += (double)Y[i];
	    pcmean = vtemp[i]+gamma[0]-gamma[1];
	    pnmean = vtemp[i];
	  }
	  ndraw = rnorm(pnmean, 1);
	  cdraw = rnorm(pcmean, 1);
	  if (*insample && R[i]==0) {
	    if (Z[i] == 1)
	      dtemp = (double)(Y[i]==0) - (double)(ndraw < 0);
	    else
	      dtemp = (double)(cdraw < 0) - (double)(Y[i]==0);
	  }
	  else 
	    dtemp = pnorm(0, pcmean, 1, 1, 0) - pnorm(0, pnmean, 1, 1, 0);
	  ITTc[0] += dtemp;
	  if (Ymax == 1) { /* binary probit */
	    if (*insample && R[i]==0) {
	      if (Z[i] == 1)
		dtemp = (double)Y[i] - (double)(ndraw > 0);
	      else
		dtemp = (double)(cdraw > 0) - (double)Y[i];
	    }
	    else
	      dtemp = pnorm(0, pcmean, 1, 0, 0) - pnorm(0, pnmean, 1, 0, 0);
	    ITTc[1] += dtemp;
	  }
	  else { /* ordered probit */
	    if (*insample && R[i]==0) {
	      if (Z[i] == 1)
		dtemp = (double)(Y[i]==Ymax) - (double)(ndraw > tau[Ymax-1]);
	      else
		dtemp = (double)(cdraw > tau[Ymax-1]) - (double)(Y[i]==Ymax);
	    }
	    else 
	      dtemp = pnorm(tau[Ymax-1], pcmean, 1, 0, 0) -
		pnorm(tau[Ymax-1], pnmean, 1, 0, 0);
	    ITTc[Ymax] += dtemp; 
	    for (j = 1; j < Ymax; j++) {
	      if (*insample && R[i]==0) {
		if (Z[i] == 1)
		  dtemp = (double)(Y[i]==j) - (double)(ndraw < tau[j] && ndraw > tau[j-1]);
		else
		  dtemp = (pnorm(tau[j], pcmean, 1, 1, 0) - 
			   pnorm(tau[j-1], pcmean, 1, 1, 0)) - (double)(Y[i]==j);
	      }
	      else
		dtemp = (pnorm(tau[j], pcmean, 1, 1, 0) - 
			 pnorm(tau[j-1], pcmean, 1, 1, 0)) 
		  - (pnorm(tau[j], pnmean, 1, 1, 0) - 
		     pnorm(tau[j-1], pnmean, 1, 1, 0));
	      ITTc[j] += dtemp;
	    }
	  }
	}
	else
	  if (Z[i] == 0) {
	    n_ncompC++;
	    base[1] += (double)Y[i];
	  }
      } 
    }
    
    /** storing the results **/
    if (main_loop > *iBurnin) {
      if (keep == *iKeep) {
	pdStore[itemp++]=(double)n_comp/(double)n_samp;
	if (Ymax == 1) {
	  pdStore[itemp++]=ITTc[1]/(double)n_comp;
	  pdStore[itemp++]=ITTc[1]/(double)n_samp;
	  pdStore[itemp++] = base[0]/(double)n_compC;
	  pdStore[itemp++] = base[1]/(double)n_ncompC;
	  pdStore[itemp++] = (base[0]+base[1])/(double)(n_compC+n_ncompC);
	}
	else {
	  for (i = 0; i <= Ymax; i++) 
	    pdStore[itemp++]=ITTc[i]/(double)n_comp;
	  for (i = 0; i <= Ymax; i++) 
	    pdStore[itemp++]=ITTc[i]/(double)n_samp;
	}
	if (*param) {
	  for(i = 0; i < n_cov; i++) 
	    pdStore[itemp++]=beta[i];
	  if (*smooth) {
	    for(i = 0; i < n_covoX; i++)
	      pdStore[itemp++]=gamma[i];
	    for(i = 0; i < n11; i++)
	      pdStore[itemp++]=treat[i];
	  }
	  else
	    for(i = 0; i < n_covo; i++)
	      pdStore[itemp++]=gamma[i];
	  if (Ymax > 1)
	    for (i = 0; i < Ymax; i++)
	      pdStore[itemp++]=tau[i];
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
    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  FreeMatrix(X, n_samp+n_cov);
  FreeMatrix(Xo, n_samp+n_covo);
  free(W);
  free(beta);
  free(gamma);
  free(q);
  free(pc);
  free(pn);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(SSo, n_covo+1);
  free(meanb);
  free(meano);
  free(meanr);
  FreeMatrix(V, n_cov);
  FreeMatrix(Vo, n_covo);
  FreeMatrix(Vr, 3);
  FreeMatrix(A, n_cov);
  FreeMatrix(Ao, n_covo);
  free(tau);
  free(taumax);
  free(taumin);
  free(ITTc);
  free(vtemp);
  free(treat);
  free(base);
  FreeMatrix(mtemp, n_cov);
  FreeMatrix(mtempo, n_covo);

} /* main */

