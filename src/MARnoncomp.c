#include <stdio.h>      
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void MARprobit(int *Y, /* binary outcome variable */ 
	       int *Ymiss, /* missingness indicator for Y */
	       int *iYmax,  /* maximum value of Y; 0,1,...,Ymax */
	       int *Z, /* treatment assignment */
	       int *D, /* treatment status */ 
	       int *C, /* compliance status */
	       double *dX, double *dXo, /* covariates */
	       double *dBeta, double *dGamma,    /* coefficients */
	       int *iNsamp, int *iNgen, int *iNcov, int *iNcovo,
	       int *iNcovoX, int *iN11, 
	       /* counters */
	       double *beta0, double *gamma0, double *dA, double *dAo, /*prior */
	       int *insample, /* 1: insample inference, 2: conditional inference */
	       int *smooth,   
	       int *param, int *mda, int *iBurnin, 
	       int *iKeep, int *verbose, /* options */
	       double *pdStore
	       ) {
  
  /*** counters ***/
  int n_samp = *iNsamp;   /* sample size */
  int n_gen = *iNgen;     /* number of gibbs draws */
  int n_cov = *iNcov;     /* number of covariates */
  int n_covo = *iNcovo;   /* number of all covariates for outcome model */
  int n_covoX = *iNcovoX; /* number of covariates excluding smooth
			     terms */
  int n11 = *iN11;        /* number of compliers in the treament group */
  
  /*** data ***/
  double **X;     /* covariates for the compliance model */
  double **Xo;    /* covariates for the outcome model */
  double *W;      /* latent variable */
  int Ymax = *iYmax;

  /*** model parameters ***/
  double *beta;   /* coef for compliance model */
  double *gamma;  /* coef for outcomme model */
  double *q;      /* some parameters for sampling C */
  double *pc; 
  double *pn;
  double pcmean;
  double pnmean;
  double **SS;    /* matrix folders for SWEEP */
  double **SSo; 
  double **SSr;
  double *meanb;  /* means for beta and gamma */
  double *meano;
  double *meanr;
  double **V;     /* variances for beta and gamma */
  double **Vo;
  double **Vr;
  double **A;
  double **Ao;
  double *tau;    /* thresholds: tau_0, ..., tau_{Ymax-1} */
  double *taumax; /* corresponding max and min for tau */
  double *taumin; /* tau_0 is fixed to 0 */
  double *treat;  /* smooth function for treat */

  /*** quantities of interest ***/
  int n_comp, n_compC, n_ncompC; 
  double *ITTc;
  double *base;

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, main_loop;  
  int itemp, itempP = ftrunc((double) n_gen/10);
  double dtemp, ndraw, cdraw;
  double *vtemp;
  double **mtemp, **mtempo;

  /*** marginal data augmentation ***/
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;

  /*** get random seed **/
  GetRNGstate();


  /*** define vectors and matricies **/
  X = doubleMatrix(n_samp+n_cov, n_cov+1);
  Xo = doubleMatrix(n_samp+n_covo, n_covo+1);
  W = doubleArray(n_samp);
  tau = doubleArray(Ymax);
  taumax = doubleArray(Ymax);
  taumin = doubleArray(Ymax);
  SS = doubleMatrix(n_cov+1, n_cov+1);
  SSo = doubleMatrix(n_covo+1, n_covo+1);
  SSr = doubleMatrix(4, 4);
  V = doubleMatrix(n_cov, n_cov);
  Vo = doubleMatrix(n_covo, n_covo);
  Vr = doubleMatrix(3, 3);
  beta = doubleArray(n_cov); 
  gamma = doubleArray(n_covo); 
  meanb = doubleArray(n_cov); 
  meano = doubleArray(n_covo); 
  meanr = doubleArray(3); 
  q = doubleArray(n_samp); 
  pc = doubleArray(n_samp); 
  pn = doubleArray(n_samp); 
  A = doubleMatrix(n_cov, n_cov);
  Ao = doubleMatrix(n_covo, n_covo);
  vtemp = doubleArray(n_samp);
  mtemp = doubleMatrix(n_cov, n_cov);
  mtempo = doubleMatrix(n_covo, n_covo);
  ITTc = doubleArray(Ymax+1);
  treat = doubleArray(n11);
  base = doubleArray(2);

  /*** read the data ***/
  itemp = 0;
  for (j =0; j < n_cov; j++)
    for (i = 0; i < n_samp; i++)
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j =0; j < n_covo; j++)
    for (i = 0; i < n_samp; i++)
      Xo[i][j] = dXo[itemp++];
  
  /*** read the prior and it as additional data points ***/ 
  itemp = 0;
  for (k = 0; k < n_cov; k++)
    for (j = 0; j < n_cov; j++)
      A[j][k] = dA[itemp++];

  itemp = 0;
  for (k = 0; k < n_covo; k++)
    for (j = 0; j < n_covo; j++)
      Ao[j][k] = dAo[itemp++];

  dcholdc(A, n_cov, mtemp);
  for(i = 0; i < n_cov; i++) {
    X[n_samp+i][n_cov]=0;
    for(j = 0; j < n_cov; j++) {
      X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
      X[n_samp+i][j] = mtemp[i][j];
    }
  }

  dcholdc(Ao, n_covo, mtempo);
  for(i = 0; i < n_covo; i++) {
    Xo[n_samp+i][n_covo]=0;
    for(j = 0; j < n_covo; j++) {
      Xo[n_samp+i][n_covo] += mtempo[i][j]*gamma0[j];
      Xo[n_samp+i][j] = mtempo[i][j];
    }
  }

  /*** starting values ***/
  for (i = 0; i < n_cov; i++) 
    beta[i] = dBeta[i];
  for (i = 0; i < n_covo; i++)
    gamma[i] = dGamma[i];
  
  if (Ymax > 1) {
    tau[0] = 0.0;
    taumax[0] = 0.0;
    taumin[0] = 0.0;
    for (i = 1; i < Ymax; i++)
      tau[i] = tau[i-1]+2/(double)(Ymax-1);
  }
  for (i = 0; i < n_samp; i++) {
    pc[i] = unif_rand(); 
    pn[i] = unif_rand();
  }

  /*** Gibbs Sampler! ***/
  itemp=0;     
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /** COMPLIANCE MODEL **/    
    if (*mda) sig2 = s0/rchisq((double)nu0);
    /* Draw complier status for control group */
    for(i = 0; i < n_samp; i++){
      dtemp = 0;
      for(j = 0; j < n_cov; j++) 
	dtemp += X[i][j]*beta[j];
      if(Z[i] == 0){
	q[i] = pnorm(dtemp, 0, 1, 1, 0);
	if(unif_rand() < (q[i]*pc[i]/(q[i]*pc[i]+(1-q[i])*pn[i]))) { 
	  C[i] = 1; Xo[i][1] = 1; 
	}
	else {
	  C[i] = 0; Xo[i][1] = 0;
	}
      }
      /* Sample W */
      if(C[i]==0) 
	W[i] = TruncNorm(dtemp-100,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+100,dtemp,1,0);
      X[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0; i < n_samp+n_cov; i++)
      for(j = 0; j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);
    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      meanb[j] = SS[j][n_cov];
    if (*mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k] = -SS[j][k]*sig2;
    rMVN(beta, meanb, V, n_cov);

    /* rescale the parameters */
    if(*mda) {
      for (i = 0; i < n_cov; i++) beta[i] /= sqrt(sig2);
    }

    /** OUTCOME MODEL **/
    /* Sample W */
    if (Ymax > 1) { /* tau_0=0, tau_1, ... */
      for (j = 1; j < (Ymax - 1); j++) {
	taumax[j] = tau[j+1];
	taumin[j] = tau[j-1];
      }
      taumax[Ymax-1] = tau[Ymax-1]+100;
      taumin[Ymax-1] = tau[Ymax-2];
    }
    if (*mda) sig2 = s0/rchisq((double)nu0);
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_covo; j++) dtemp += Xo[i][j]*gamma[j];
      if (Ymiss[i] == 1) {
	W[i] = dtemp + norm_rand();
	if (Ymax == 1) { /* binary probit */
	  if (W[i] > 0) Y[i] = 1;
	  else Y[i] = 0;
	}
	else { /* ordered probit */
	  if (W[i] >= tau[Ymax-1])
	    Y[i] = Ymax;
	  else {
	    j = 0;
	    while (W[i] > tau[j] && j < Ymax) j++;
	    Y[i] = j;
	  }
	}
      }
      else {
	if(Ymax == 1) { /* binary probit */
	  if(Y[i] == 0) W[i] = TruncNorm(dtemp-100,0,dtemp,1,0);
	  else W[i] = TruncNorm(0,dtemp+100,dtemp,1,0);
	}
	else {         /* ordered probit */
	  if (Y[i] == 0) 
	    W[i] = TruncNorm(dtemp-100, 0, dtemp, 1, 0);
	  else if (Y[i] == Ymax) {
	    W[i] = TruncNorm(tau[Ymax-1], dtemp+100, dtemp, 1, 0);
	    if (W[i] < taumax[Ymax-1]) taumax[Ymax-1] = W[i];
	  }
	  else {
	    W[i] = TruncNorm(tau[Y[i]-1], tau[Y[i]], dtemp, 1, 0);
	    if (W[i] > taumin[Y[i]]) taumin[Y[i]] = W[i];
	    if (W[i] < taumax[Y[i]-1]) taumax[Y[i]-1] = W[i];
	  }
	}
      }
      Xo[i][n_covo] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }
    /* draw tau */
    if (Ymax > 1) 
      for (j = 1; j < Ymax; j++) 
	tau[j] = runif(taumin[j], taumax[j])*sqrt(sig2);
    /* SS matrix */
    for(j = 0; j <= n_covo; j++)
      for(k = 0; k <= n_covo; k++)
	SSo[j][k]=0;
    for(i = 0;i < n_samp+n_covo; i++)
      for(j = 0;j <= n_covo; j++)
	for(k = 0; k <= n_covo; k++) 
	  SSo[j][k] += Xo[i][j]*Xo[i][k];
    /* SWEEP SS matrix */
    for(j = 0; j < n_covo; j++)
      SWP(SSo, j, n_covo+1);

    /* draw gamma */    
    for(j = 0; j < n_covo; j++)
      meano[j] = SSo[j][n_covo];
    if (*mda) 
      sig2=(SSo[n_covo][n_covo]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_covo; j++)
      for(k = 0; k < n_covo; k++) Vo[j][k]=-SSo[j][k]*sig2;
    rMVN(gamma, meano, Vo, n_covo); 
    
    /* rescaling the parameters */
    if(*mda) {
      for (i = 0; i < n_covo; i++) gamma[i] /= sqrt(sig2);
      if (Ymax > 1)
	for (i = 1; i < Ymax; i++)
	  tau[i] /= sqrt(sig2);
    }

    /* computing smooth terms */
    if (*smooth) {
      for (i = 0; i < n11; i++) {
	treat[i] = 0;
	for (j = n_covoX; j < n_covo; j++)
	  treat[i] += Xo[i][j]*gamma[j]; 
      }
    }

    /** Compute probabilities **/ 
    for(i = 0; i < n_samp; i++){
      vtemp[i] = 0;
      for(j = 0; j < n_covo; j++)
	vtemp[i] += Xo[i][j]*gamma[j];
    }

    for(i = 0; i < n_samp; i++){
      if(Z[i]==0){
	if (C[i] == 1) {
	  pcmean = vtemp[i];
	  if (*smooth)
	    pnmean = vtemp[i]-gamma[0];
	  else
	    pnmean = vtemp[i]-gamma[1];
	}
	else {
	  if (*smooth)
	    pcmean = vtemp[i]+gamma[0];
	  else
	    pcmean = vtemp[i]+gamma[1];
	  pnmean = vtemp[i];
	}
	if (Y[i] == 0){
	  pc[i] = pnorm(0, pcmean, 1, 1, 0);
	  pn[i] = pnorm(0, pnmean, 1, 1, 0);
	}
	else {
	  if (Ymax == 1) { /* binary probit */
	    pc[i] = pnorm(0, pcmean, 1, 0, 0);
	    pn[i] = pnorm(0, pnmean, 1, 0, 0);
	  }
	  else { /* ordered probit */
	    if (Y[i] == Ymax) {
	      pc[i] = pnorm(tau[Ymax-1], pcmean, 1, 0, 0);
	      pn[i] = pnorm(tau[Ymax-1], pnmean, 1, 0, 0);
	    }
	    else {
	      pc[i] = pnorm(tau[Y[i]], pcmean, 1, 1, 0) -
		pnorm(tau[Y[i]-1], pcmean, 1, 1, 0);
	      pn[i] = pnorm(tau[Y[i]], pnmean, 1, 1, 0) - 
		pnorm(tau[Y[i]-1], pnmean, 1, 1, 0);
	    }
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
	  if (*insample && Ymiss[i]==0) 
	    dtemp = (double)(Y[i]==0) - (double)(ndraw < 0);
	  else
	    dtemp = pnorm(0, pcmean, 1, 1, 0) - pnorm(0, pnmean, 1, 1, 0);
	  ITTc[0] += dtemp;
	  if (Ymax == 1) { /* binary probit */
	    if (*insample && Ymiss[i]==0) 
	      dtemp = (double)Y[i] - (double)(ndraw > 0);
	    else
	      dtemp = pnorm(0, pcmean, 1, 0, 0) - pnorm(0, pnmean, 1, 0, 0);
	    ITTc[1] += dtemp;
	  }
	  else { /* ordered probit */
	    if (*insample && Ymiss[i]==0) 
	      dtemp = (double)(Y[i]==Ymax) - (double)(ndraw > tau[Ymax-1]);
	    else
	      dtemp = pnorm(tau[Ymax-1], pcmean, 1, 0, 0) -
		pnorm(tau[Ymax-1], pnmean, 1, 0, 0);
	    ITTc[Ymax] += dtemp; 
	    for (j = 1; j < Ymax; j++) {
	      if (*insample && Ymiss[i]==0)
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
	  if (*insample && Ymiss[i]==0) {
	    if (Z[i] == 1)
	      dtemp = (double)(Y[i]==0) - (double)(ndraw < 0);
	    else
	      dtemp = (double)(cdraw < 0) - (double)(Y[i]==0);
	  }
	  else 
	    dtemp = pnorm(0, pcmean, 1, 1, 0) - pnorm(0, pnmean, 1, 1, 0);
	  ITTc[0] += dtemp;
	  if (Ymax == 1) { /* binary probit */
	    if (*insample && Ymiss[i]==0) {
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
	    if (*insample && Ymiss[i]==0) {
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
	      if (*insample && Ymiss[i]==0) {
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



