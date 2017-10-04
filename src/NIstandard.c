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
   Bayesian Binary Probit with Nonignorable Missing Outcomes 
*/

void NIbprobit(int *Y,         /* binary outcome variable */ 
	       int *R,         /* recording indicator for Y */
	       double *dXo,    /* covariates */
	       double *dXr,    /* covariates */
	       double *beta,   /* coefficients */
	       double *delta,  /* coefficients */
	       int *insamp,    /* # of obs */ 
	       int *incovo,    /* # of covariates */
	       int *incovr,    /* # of covariates */
	       int *intreat,   /* # of treatments */
	       double *beta0,  /* prior mean */
	       double *delta0, /* prior mean */
	       double *dAo,    /* prior precision */
	       double *dAr,    /* prior precision */
	       int *Insample,  /* insample QoI */
	       int *param,     /* store parameters? */ 
	       int *mda,       /* marginal data augmentation? */ 
	       int *ndraws,    /* # of gibbs draws */
	       int *iBurnin,   /* # of burnin */
	       int *iKeep,     /* every ?th draws to keep */
	       int *verbose,  
	       double *coefo,  /* storage for coefficients */ 
	       double *coefr,  /* storage for coefficients */ 
	       double *ATE,     /* storage for ATE */
	       double *BASE    /* storage for baseline */
	       ) {
  
  /*** counters ***/
  int n_samp = *insamp;      /* sample size */
  int n_gen = *ndraws;       /* number of gibbs draws */
  int n_covo = *incovo;      /* number of covariates */
  int n_covr = *incovr;      /* number of covariates */
  int n_treat = *intreat;    /* number of treatments */

  /*** data ***/
  /* covariates for the response model */
  double **Xr = doubleMatrix(n_samp+n_covr, n_covr+1);
  /* covariates for the outcome model */     
  double **Xo = doubleMatrix(n_samp+n_covo, n_covo+1);

  /*** model parameters ***/
  double **Ao = doubleMatrix(n_covo, n_covo);
  double **Ar = doubleMatrix(n_covr, n_covr);
  double **mtemp1 = doubleMatrix(n_covo, n_covo);
  double **mtemp2 = doubleMatrix(n_covr, n_covr);

  /*** QoIs ***/
  double *base = doubleArray(n_treat);
  double *cATE = doubleArray(n_treat);

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, main_loop;  
  int itemp, itemp0, itemp1, itemp2, itempP = ftrunc((double) n_gen/10);
  double dtemp, pj, r0, r1;

  /*** get random seed **/
  GetRNGstate();

  /*** read the data ***/
  itemp = 0;
  for (j = 0; j < n_covo; j++)
    for (i = 0; i < n_samp; i++) 
      Xo[i][j] = dXo[itemp++];
  itemp = 0;
  for (j = 0; j < n_covr; j++)
    for (i = 0; i < n_samp; i++) 
      Xr[i][j] = dXr[itemp++];
  
  /*** read the prior and it as additional data points ***/ 
  itemp = 0;
  for (k = 0; k < n_covo; k++)
    for (j = 0; j < n_covo; j++)
      Ao[j][k] = dAo[itemp++];

  itemp = 0;
  for (k = 0; k < n_covr; k++)
    for (j = 0; j < n_covr; j++)
      Ar[j][k] = dAr[itemp++];

  dcholdc(Ao, n_covo, mtemp1);
  for(i = 0; i < n_covo; i++) {
    Xo[n_samp+i][n_covo] = 0;
    for(j = 0; j < n_covo; j++) {
      Xo[n_samp+i][n_covo] += mtemp1[i][j]*beta0[j];
      Xo[n_samp+i][j] = mtemp1[i][j];
    }
  }

  dcholdc(Ar, n_covr, mtemp2);
  for(i = 0; i < n_covr; i++) {
    Xr[n_samp+i][n_covr] = 0;
    for(j = 0; j < n_covr; j++) {
      Xr[n_samp+i][n_covr] += mtemp2[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtemp2[i][j];
    }
  }

  /*** Gibbs Sampler! ***/
  itemp = 0; itemp0 = 0; itemp1 = 0; itemp2 = 0;     
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /** Response Model: binary Probit **/    
    bprobitGibbs(R, Xr, delta, n_samp, n_covr, 0, delta0, Ar, *mda, 1);
      
    /** Outcome Model: binary probit **/
    bprobitGibbs(Y, Xo, beta, n_samp, n_covo, 0, beta0, Ao, *mda, 1);

    /** Imputing the missing data **/
    for (i = 0; i < n_samp; i++) {
      if (R[i] == 0) {
	pj = 0;
	r0 = delta[0];
	r1 = delta[1];
	for (j = 0; j < n_covo; j++) 
	  pj += Xo[i][j]*beta[j];
	for (j = 2; j < n_covr; j++) {
	  r0 += Xr[i][j]*delta[j];
	  r1 += Xr[i][j]*delta[j];
	}
	pj = pnorm(0, pj, 1, 0, 0);
	r0 = pnorm(0, r0, 1, 0, 0);
	r1 = pnorm(0, r1, 1, 0, 0);
	if (unif_rand() < ((1-r1)*pj/((1-r1)*pj+(1-r0)*(1-pj)))) {
	  Y[i] = 1;
	  Xr[i][0] = 0;
	  Xr[i][1] = 1;
	} else {
	  Y[i] = 0;
	  Xr[i][0] = 1;
	  Xr[i][1] = 0;
	} 
      }
    }
    
    /** Compute quantities of interest **/
    for (j = 0; j < n_treat; j++) 
      base[j] = 0;
    for (i = 0; i < n_samp; i++) {
      dtemp = 0; 
      for (j = n_treat; j < n_covo; j++) 
	dtemp += Xo[i][j]*beta[j];
      for (j = 0; j < n_treat; j++) {
	if (*Insample) {
	  if (Xo[i][j] == 1)
	    base[j] += (double)Y[i];
	  else
	    base[j] += (double)((dtemp+beta[j]+norm_rand()) > 0);
	} else
	  base[j] += pnorm(0, dtemp+beta[j], 1, 0, 0);
      }
    }
    for (j = 0; j < n_treat; j++) 
      base[j] /= (double)n_samp;
    
    /** Storing the results **/
    if (main_loop > *iBurnin) {
      if (keep == *iKeep) {
	for (j = 0; j < (n_treat-1); j++)
	  ATE[itemp0++] = base[j+1] - base[0];
	for (j = 0; j < n_treat; j++)
	  BASE[itemp++] = base[j];
	if (*param) {
	  for (i = 0; i < n_covo; i++) 
	    coefo[itemp1++] = beta[i];
	  for (i = 0; i < n_covr; i++) 
	    coefr[itemp2++] = delta[i];
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
  FreeMatrix(Xr, n_samp+n_covr);
  FreeMatrix(Xo, n_samp+n_covo);
  FreeMatrix(Ao, n_covo);
  FreeMatrix(Ar, n_covr);
  FreeMatrix(mtemp1, n_covo);
  FreeMatrix(mtemp2, n_covr);
  free(base);
  free(cATE);
} /* NIbprobit */




/* 
   Bayesian Binary Mixed Effects Probit with Nonignorable Missing
   Outcomes  
*/

void NIbprobitMixed(int *Y,         /* binary outcome variable */ 
		    int *R,         /* recording indicator for Y */
		    int *grp,       /* group indicator */
		    int *in_grp,    /* number of groups */
		    int *max_samp_grp, /* max # of obs within group */
		    double *dXo,    /* fixed effects covariates */
		    double *dXr,    /* fixed effects covariates */
		    double *dZo,    /* random effects covariates */
		    double *dZr,    /* random effects  covariates */
		    double *beta,   /* coefficients */
		    double *delta,  /* coefficients */
		    double *dPsio,  /* random effects variance */
		    double *dPsir,  /* random effects variance */
		    int *insamp,    /* # of obs */ 
		    int *incovo,    /* # of fixed effects */
		    int *incovr,    /* # of fixed effects */
		    int *incovoR,   /* # of random effects */
		    int *incovrR,   /* # of random effects */
		    int *intreat,   /* # of treatments */
		    double *beta0,  /* prior mean */
		    double *delta0, /* prior mean */
		    double *dAo,    /* prior precision */
		    double *dAr,    /* prior precision */
		    int *dfo,       /* prior degrees of freedom */
		    int *dfr,       /* prior degrees of freedom */
		    double *dS0o, /* prior scale */
		    double *dS0r, /* prior scale */
		    int *Insample,  /* insample QoI */
		    int *param,     /* store parameters? */ 
		    int *mda,       /* marginal data augmentation? */ 
		    int *ndraws,    /* # of gibbs draws */
		    int *iBurnin,   /* # of burnin */
		    int *iKeep,     /* every ?th draws to keep */
		    int *verbose,  
		    double *coefo,  /* storage for coefficients */ 
		    double *coefr,  /* storage for coefficients */ 
		    double *sPsiO,   /* storage for variance */
		    double *sPsiR,   /* storage for variance */
		    double *ATE,     /* storage for ATE */
		    double *BASE    /* storage for baseline */
		    ) {
  
  /*** counters ***/
  int n_samp = *insamp;      /* sample size */
  int n_gen = *ndraws;       /* number of gibbs draws */
  int n_grp = *in_grp;       /* number of groups */
  int n_covo = *incovo;      /* number of fixed effects */
  int n_covr = *incovr;      /* number of fixed effects */
  int n_covoR = *incovoR;    /* number of random effects */
  int n_covrR = *incovrR;    /* number of random effects */
  int n_treat = *intreat;    /* number of treatments */

  /*** data ***/
  /* covariates for the response model */
  double **Xr = doubleMatrix(n_samp+n_covr, n_covr+1);
  /* covariates for the outcome model */     
  double **Xo = doubleMatrix(n_samp+n_covo, n_covo+1);
  /* random effects covariates */
  double ***Zo = doubleMatrix3D(n_grp, *max_samp_grp + n_covoR,
				n_covoR + 1);
  double ***Zr = doubleMatrix3D(n_grp, *max_samp_grp + n_covrR,
				n_covrR + 1);

  /*** model parameters ***/
  double **PsiO = doubleMatrix(n_covoR, n_covoR);
  double **PsiR = doubleMatrix(n_covrR, n_covrR);
  double **xiO = doubleMatrix(n_grp, n_covoR);
  double **xiR = doubleMatrix(n_grp, n_covrR);
  double **S0o = doubleMatrix(n_covoR, n_covoR);
  double **S0r = doubleMatrix(n_covrR, n_covrR);
  double **Ao = doubleMatrix(n_covo, n_covo);
  double **Ar = doubleMatrix(n_covr, n_covr);
  double **mtemp1 = doubleMatrix(n_covo, n_covo);
  double **mtemp2 = doubleMatrix(n_covr, n_covr);

  /*** QoIs ***/
  double *base = doubleArray(n_treat);
  double *cATE = doubleArray(n_treat);

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, main_loop;  
  int itemp, itemp0, itemp1, itemp2, itemp3 = 0, itempP = ftrunc((double) n_gen/10);
  int *vitemp = intArray(n_grp);
  double dtemp, pj, r0, r1;

  /*** get random seed **/
  GetRNGstate();

  /*** fixed effects ***/
  itemp = 0;
  for (j = 0; j < n_covo; j++)
    for (i = 0; i < n_samp; i++) 
      Xo[i][j] = dXo[itemp++];
  itemp = 0;
  for (j = 0; j < n_covr; j++)
    for (i = 0; i < n_samp; i++) 
      Xr[i][j] = dXr[itemp++];
  
  /* prior */
  itemp = 0;
  for (k = 0; k < n_covo; k++)
    for (j = 0; j < n_covo; j++)
      Ao[j][k] = dAo[itemp++];

  itemp = 0;
  for (k = 0; k < n_covr; k++)
    for (j = 0; j < n_covr; j++)
      Ar[j][k] = dAr[itemp++];

  dcholdc(Ao, n_covo, mtemp1);
  for(i = 0; i < n_covo; i++) {
    Xo[n_samp+i][n_covo] = 0;
    for(j = 0; j < n_covo; j++) {
      Xo[n_samp+i][n_covo] += mtemp1[i][j]*beta0[j];
      Xo[n_samp+i][j] = mtemp1[i][j];
    }
  }

  dcholdc(Ar, n_covr, mtemp2);
  for(i = 0; i < n_covr; i++) {
    Xr[n_samp+i][n_covr] = 0;
    for(j = 0; j < n_covr; j++) {
      Xr[n_samp+i][n_covr] += mtemp2[i][j]*delta0[j];
      Xr[n_samp+i][j] = mtemp2[i][j];
    }
  }

  /* random effects */
  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    for (j = 0; j < n_covoR; j++)
      Zo[grp[i]][vitemp[grp[i]]][j] = dZo[itemp++];
    vitemp[grp[i]]++;
  }

  itemp = 0;
  for (j = 0; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    for (j = 0; j < n_covrR; j++)
      Zr[grp[i]][vitemp[grp[i]]][j] = dZr[itemp++];
    vitemp[grp[i]]++;
  }

  /* prior variance for random effects */
  itemp = 0;
  for (k = 0; k < n_covoR; k++)
    for (j = 0; j < n_covoR; j++) 
      PsiO[j][k] = dPsio[itemp++];

  itemp = 0;
  for (k = 0; k < n_covrR; k++)
    for (j = 0; j < n_covrR; j++) 
      PsiR[j][k] = dPsir[itemp++];

  itemp = 0;
  for (k = 0; k < n_covoR; k++)
    for (j = 0; j < n_grp; j++)
      xiO[j][k] = norm_rand();

  itemp = 0;
  for (k = 0; k < n_covrR; k++)
    for (j = 0; j < n_grp; j++)
      xiR[j][k] = norm_rand();

  /* hyper prior scale parameter for random effects */
  itemp = 0;
  for (k = 0; k < n_covoR; k++)
    for (j = 0; j < n_covoR; j++)
      S0o[j][k] = dS0o[itemp++];

  itemp = 0;
  for (k = 0; k < n_covrR; k++)
    for (j = 0; j < n_covrR; j++)
      S0r[j][k] = dS0r[itemp++];

  /*** Gibbs Sampler! ***/
  itemp = 0; itemp0 = 0; itemp1 = 0; itemp2 = 0;     
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /** Response Model: binary Probit **/    
    bprobitMixedGibbs(R, Xr, Zr, grp, delta, xiR, PsiR, n_samp,
		      n_covr, n_covrR, n_grp, 0, delta0, Ar, *dfr, S0r,
		      1);
      
    /** Outcome Model: binary probit **/
    bprobitMixedGibbs(Y, Xo, Zr, grp, beta, xiO, PsiO, n_samp, n_covo,
		      n_covoR, n_grp, 0, beta0, Ao, *dfo, S0o, 1);

    /** Imputing the missing data **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      if (R[i] == 0) {
	pj = 0;
	r0 = delta[0];
	r1 = delta[1];
	for (j = 0; j < n_covo; j++) 
	  pj += Xo[i][j]*beta[j];
	for (j = 2; j < n_covr; j++) {
	  r0 += Xr[i][j]*delta[j];
	  r1 += Xr[i][j]*delta[j];
	}
	for (j = 0; j < n_covoR; j++)
	  pj += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
	for (j = 0; j < n_covrR; j++) {
	  r0 += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	  r1 += Zr[grp[i]][vitemp[grp[i]]][j]*xiR[grp[i]][j];
	}
	pj = pnorm(0, pj, 1, 0, 0);
	r0 = pnorm(0, r0, 1, 0, 0);
	r1 = pnorm(0, r1, 1, 0, 0);
	if (unif_rand() < ((1-r1)*pj/((1-r1)*pj+(1-r0)*(1-pj)))) {
	  Y[i] = 1;
	  Xr[i][0] = 0;
	  Xr[i][1] = 1;
	} else {
	  Y[i] = 0;
	  Xr[i][0] = 1;
	  Xr[i][1] = 0;
	} 
      }
      vitemp[grp[i]]++;
    }
    
    /** Compute quantities of interest **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (j = 0; j < n_treat; j++) 
      base[j] = 0;
    for (i = 0; i < n_samp; i++) {
      dtemp = 0; 
      for (j = n_treat; j < n_covo; j++) 
	dtemp += Xo[i][j]*beta[j];
      for (j = 0; j < n_covoR; j++)
	dtemp += Zo[grp[i]][vitemp[grp[i]]][j]*xiO[grp[i]][j];
      for (j = 0; j < n_treat; j++) {
	if (*Insample) {
	  if (Xo[i][j] == 1)
	    base[j] += (double)Y[i];
	  else
	    base[j] += (double)((dtemp+beta[j]+norm_rand()) > 0);
	} else
	  base[j] += pnorm(0, dtemp+beta[j], 1, 0, 0);
      }
      vitemp[grp[i]]++;
    }
    for (j = 0; j < n_treat; j++) 
      base[j] /= (double)n_samp;
    
    /** Storing the results **/
    if (main_loop > *iBurnin) {
      if (keep == *iKeep) {
	for (j = 0; j < (n_treat-1); j++)
	  ATE[itemp0++] = base[j+1] - base[0];
	for (j = 0; j < n_treat; j++)
	  BASE[itemp++] = base[j];
	if (*param) {
	  for (i = 0; i < n_covo; i++) 
	    coefo[itemp1++] = beta[i];
	  for (i = 0; i < n_covr; i++) 
	    coefr[itemp2++] = delta[i];
	  for (i = 0; i < n_covoR; i++)
	    for (j = i; j < n_covoR; j++)
	      sPsiO[itemp3++] = PsiO[i][j];
	  for (i = 0; i < n_covrR; i++)
	    for (j = i; j < n_covrR; j++)
	      sPsiR[itemp3++] = PsiR[i][j];
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
  FreeMatrix(Xr, n_samp+n_covr);
  FreeMatrix(Xo, n_samp+n_covo);
  Free3DMatrix(Zo, n_grp, *max_samp_grp + n_covoR);
  Free3DMatrix(Zr, n_grp, *max_samp_grp + n_covrR);
  FreeMatrix(PsiO, n_covoR);
  FreeMatrix(PsiR, n_covrR);
  FreeMatrix(xiO, n_grp);
  FreeMatrix(xiR, n_grp);
  FreeMatrix(S0o, n_covoR);
  FreeMatrix(S0r, n_covrR);
  FreeMatrix(Ao, n_covo);
  FreeMatrix(Ar, n_covr);
  FreeMatrix(mtemp1, n_covo);
  FreeMatrix(mtemp2, n_covr);
  free(base);
  free(cATE);
  free(vitemp);
} /* NIbprobitMixed */



