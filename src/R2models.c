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

void R2bNormalReg(double *Y,       /* binary outcome variable */
		  double *dX,      /* model matrix */
		  double *beta,    /* fixed effects coefficients */
		  double *sig2,    /* variance parameter */
		  int *n_samp,     /* # of obs */ 
		  int *n_cov,      /* # of covariates */
		  int *n_gen,      /* # of gibbs draws */
		  int *pbeta,      /* 0: improper prior 
				      p(beta|X) \propto 1
				      1: proper prior for (beta)
				      p(beta|X) = normal(beta0, A0)
				   */
		  double *beta0,   /* prior mean for normal */
		  double *dA0,     /* prior precision for normal; can be
				      set to zero to induce improper prior
				      for beta alone
				   */
		  int *psig2,      /* 0: improper prior for sig2
				      p(sig2|X) \propto 1/sig2
				      1: proper prior for sig2
				      p(sigma2|X) = InvChi2(nu0, s0)
				   */
		  double *s0,      /* prior scale for InvChi2 */
		  int *nu0,        /* prior d.f. for InvChi2 */
		  int *sig2fixed,  /* 1: sig2 fixed, 0: sig2 sampled */ 
		  double *betaStore, 
		  double *sig2Store
		  ) {

  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int ibeta = 0, isig2 = 0;

  /* matrices */
  double **X = doubleMatrix(*n_samp+*n_cov, *n_cov+1);
  double **A0 = doubleMatrix(*n_cov, *n_cov);
  double **mtemp = doubleMatrix(*n_cov, *n_cov);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  /* packing the prior */
  itemp = 0; 
  for (k = 0; k < *n_cov; k++)
    for (j = 0; j < *n_cov; j++)
      A0[j][k] = dA0[itemp++];

  /* adding prior as an additional data point */
  dcholdc(A0, *n_cov, mtemp);
  for (i = 0; i < *n_samp; i++)
    X[i][*n_cov] = Y[i];
  for (i = 0; i < *n_cov; i++) {
    X[*n_samp+i][*n_cov]=0;
    for (j = 0; j < *n_cov; j++) {
      X[*n_samp+i][*n_cov] += mtemp[i][j]*beta0[j];
      X[*n_samp+i][j] = mtemp[i][j];
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    bNormalReg(X, beta, sig2, *n_samp, *n_cov, 0, *pbeta, beta0, A0,
	       *psig2, *s0, *nu0, *sig2fixed);

    /* Storing the output */
    for (j = 0; j < *n_cov; j++)
      betaStore[ibeta++] = beta[j];
    sig2Store[isig2++] = *sig2;

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  FreeMatrix(X, *n_samp+*n_cov);
  FreeMatrix(A0, *n_cov);
  FreeMatrix(mtemp, *n_cov);
}


void R2boprobit(int *Y,          /* ordinal outcome variable: 0, 1,
				    ..., J-1 */
		double *dX,      /* model matrix */
		double *beta,    /* fixed effects coefficients */
		double *tau,     /* J cut points; the first cut
				    point is set to 0 and the last one
				    set to tau_{J-1}+1000 */
		int *n_samp,     /* # of obs */ 
		int *n_cov,      /* # of covariates */
		int *n_cat,      /* # of categories: J */
		double *beta0,   /* prior mean for normal */
		double *dA0,     /* prior precision for normal: can be
				    set to zero to induce improper prior
				    for beta alone
				 */
		int *n_gen,      /* # of gibbs draws */
		int *mda,        /* marginal data augmentation? */
		int *mh,         /* Metropolis-Hasting step? */
		double *prop,    /* proposal variance for MH step */
		int *accept,     /* acceptance counter */
		double *betaStore, 
		double *tauStore
		) {
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int ibeta = 0, itau = 0; 
  double dtemp;

  /* matrices */
  double **X = doubleMatrix(*n_samp+*n_cov, *n_cov+1);
  double **A0 = doubleMatrix(*n_cov, *n_cov);
  double **mtemp = doubleMatrix(*n_cov, *n_cov);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  /* packing the prior */
  itemp = 0; 
  for (k = 0; k < *n_cov; k++)
    for (j = 0; j < *n_cov; j++)
      A0[j][k] = dA0[itemp++];

  /* adding prior as an additional data point */
  dcholdc(A0, *n_cov, mtemp);
  for (i = 0; i < *n_cov; i++) {
    X[*n_samp+i][*n_cov]=0;
    for (j = 0; j < *n_cov; j++) {
      X[*n_samp+i][*n_cov] += mtemp[i][j]*beta0[j];
      X[*n_samp+i][j] = mtemp[i][j];
    }
  }

  if (*mh) {
    for (i = 0; i < *n_samp; i++){
      dtemp = 0;
      for (j = 0; j < *n_cov; j++) 
	dtemp += X[i][j]*beta[j]; 
      if (Y[i] == 0) 
        X[i][*n_cov] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	X[i][*n_cov] = TruncNorm(tau[Y[i]-1],tau[Y[i]],dtemp,1,0);
    }
  }
  
  accept[0] = 0;
  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    boprobitMCMC(Y, X, beta, tau, *n_samp, *n_cov, *n_cat,
		 0, beta0, A0, *mda, *mh, prop, accept, 1);

    /* Storing the output */
    for (j = 0; j < *n_cov; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_cat-1; j++)
      tauStore[itau++] = tau[j];

    Rprintf("Acceptance ratio: %14g\n", (double)accept[0]/(double) main_loop);
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  FreeMatrix(X, *n_samp+*n_cov);
  FreeMatrix(A0, *n_cov);
  FreeMatrix(mtemp, *n_cov);
}



void R2logitMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		  double *dX,    /* (N x K) covariate matrix */
		  double *beta,  /* (K(J-1)) stacked coefficient vector */
		  int *n_samp,   /* # of obs */
		  int *n_dim,    /* # of categories, J-1 */
		  int *n_cov,    /* # of covariates, K */
		  double *beta0, /* (K(J-1)) prior mean vector */
		  double *dA0,   /* (K(J-1) x K(J-1)) prior precision */
		  double *Var,   /* K(J-1) proposal variances */
		  int *n_gen,     /* # of MCMC draws */
		  int *counter,  /* # of acceptance for each parameter */
		  double *store  /* Storage for beta */
		  ) {

  /* storage parameters and loop counters */
  int i, j, k, itemp, main_loop;  
  
  /* matrices */
  double **X = doubleMatrix(*n_samp, *n_cov+1);
  double **A0 = doubleMatrix(n_cov[0]*n_dim[0], n_cov[0]*n_dim[0]);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  /* packing the prior */
  itemp = 0; 
  for (k = 0; k < n_cov[0]*n_dim[0]; k++)
    for (j = 0; j < n_cov[0]*n_dim[0]; j++)
      A0[j][k] = dA0[itemp++];

  /* Gibbs Sampler! */
  itemp = 0;
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    logitMetro(Y, X, beta, *n_samp, *n_dim, *n_cov, beta0, A0,
	       Var, 1, counter);

    /* Storing the output */
    for (j = 0; j < n_dim[0]*n_cov[0]; j++)
      store[itemp++] = beta[j];

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  FreeMatrix(X, *n_samp);
  FreeMatrix(A0, *n_cov);
}



void R2bNormalMixedGibbs(double *Y,        /* outcome variable */
			 double *dX,       /* model matrix for fixed
					      effects */
			 double *dZ,       /* model matrix for random
					      effects */
			 int *grp,         /* group indicator: 0, 1, 2,... */
			 double *beta,     /* fixed effects coefficients */
			 double *dgamma,   /* random effects coefficients */
			 double *sig2,     /* variance parameter */
			 double *dPsi,     /* covariance for random
					      effects */
			 int *n_samp,      /* # of obs */ 
			 int *n_fixed,     /* # of fixed effects */
			 int *n_random,    /* # of random effects */
			 int *n_grp,       /* # of groups */
			 int *max_samp_grp, /* max # of obs within group */
			 double *beta0,    /* prior mean */
			 double *dA0,      /* prior precision */
			 int *imp,         /* do you want to use
					      improper prior for sig2
					      and Psi? */
			 int *nu0,         /* prior df for sig2 */
			 double *s0,       /* prior scale for sig2 */
			 int *tau0,        /* prior df for Psi */
			 double *dT0,      /* prior scale for Psi */
			 int *n_gen,       /* # of gibbs draws */
			 /* storage of MCMC draws */
			 double *betaStore, 
			 double *gammaStore,
			 double *sig2Store,
			 double *PsiStore
		       ) {

  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int ibeta = 0, igamma = 0, iPsi =0, isig2 = 0;

  /* matrices */
  double **X = doubleMatrix(*n_samp+*n_fixed, *n_fixed+1);
  double **gamma = doubleMatrix(*n_grp, *n_random);
  double **Psi = doubleMatrix(*n_random, *n_random);
  double **A0 = doubleMatrix(*n_fixed, *n_fixed);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double **mtemp = doubleMatrix(*n_fixed, *n_fixed);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp + *n_random, *n_random+1);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_random; j++)
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      Psi[j][k] = dPsi[itemp++];
  
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_grp; j++)
      gamma[j][k] = dgamma[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_fixed; k++)
    for (j = 0; j < *n_fixed; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  /* adding prior as an additional data point */
  dcholdc(A0, *n_fixed, mtemp);
  for (i = 0; i < *n_fixed; i++) {
    X[*n_samp+i][*n_fixed]=0;
    for (j = 0; j < *n_fixed; j++) {
      X[*n_samp+i][*n_fixed] += mtemp[i][j]*beta0[j];
      X[*n_samp+i][j] = mtemp[i][j];
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    bNormalMixedGibbs(Y, X, Zgrp, grp, beta, gamma, sig2, Psi, 
		      *n_samp, *n_fixed, *n_random, *n_grp, 
		      0, beta0, A0, *imp, *nu0, *s0, *tau0, T0, 1);

    /* Storing the output */
    for (j = 0; j < *n_fixed; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_random; j++)
      for (k = j; k < *n_random; k++)
	PsiStore[iPsi++] = Psi[j][k];
    for (j = 0; j < *n_grp; j++)
      for (k = 0; k < *n_random; k++)
	gammaStore[igamma++] = gamma[j][k];
    sig2Store[isig2++] = sig2[0];

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  FreeMatrix(X, *n_samp+*n_fixed);
  FreeMatrix(gamma, *n_grp);
  FreeMatrix(Psi, *n_random);
  FreeMatrix(A0, *n_fixed);
  FreeMatrix(T0, *n_random);
  FreeMatrix(mtemp, *n_fixed);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp + *n_random);
} /* end of normal mixed effects model */



void R2bprobitMixedGibbs(int *Y,           /* binary outcome variable */
			 double *dX,       /* model matrix for fixed
					      effects */
			 double *dZ,       /* model matrix for random
					      effects */
			 int *grp,         /* group indicator: 0, 1, 2,... */
			 double *beta,     /* fixed effects coefficients */
			 double *dPsi,     /* covariance for random
					      effects */
			 int *n_samp,      /* # of obs */ 
			 int *n_fixed,     /* # of fixed effects */
			 int *n_random,    /* # of random effects */
			 int *n_grp,       /* # of groups */
			 int *max_samp_grp, /* max # of obs within group */
			 double *beta0,    /* prior mean */
			 double *dA0,      /* prior precision */
			 int *tau0,        /* prior df */
			 double *dT0,      /* prior scale */
			 int *n_gen,       /* # of gibbs draws */
			 /* storage of MCMC draws */
			 double *betaStore, 
			 double *gammaStore,
			 double *PsiStore
		       ) {

  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int ibeta = 0, igamma = 0, iPsi =0;

  /* matrices */
  double **X = doubleMatrix(*n_samp+*n_fixed, *n_fixed+1);
  double **gamma = doubleMatrix(*n_grp, *n_random);
  double *gamma0 = doubleArray(*n_random);
  double **Psi = doubleMatrix(*n_random, *n_random);
  double **A0 = doubleMatrix(*n_fixed, *n_fixed);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double **mtemp = doubleMatrix(*n_fixed, *n_fixed);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp + *n_random, *n_random+1);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_random; j++)
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      Psi[j][k] = dPsi[itemp++];
  
  for (j = 0; j < *n_random; j++)
    gamma0[j] = 0;
  for (j = 0; j < *n_grp; j++)
    rMVN(gamma[j], gamma0, Psi, *n_random);

  itemp = 0; 
  for (k = 0; k < *n_fixed; k++)
    for (j = 0; j < *n_fixed; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  /* adding prior as an additional data point */
  dcholdc(A0, *n_fixed, mtemp);
  for (i = 0; i < *n_fixed; i++) {
    X[*n_samp+i][*n_fixed]=0;
    for (j = 0; j < *n_fixed; j++) {
      X[*n_samp+i][*n_fixed] += mtemp[i][j]*beta0[j];
      X[*n_samp+i][j] = mtemp[i][j];
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    bprobitMixedGibbs(Y, X, Zgrp, grp, beta, gamma, Psi, *n_samp,
		      *n_fixed, *n_random, *n_grp,
		      0, beta0, A0, *tau0, T0, 1);

    /* Storing the output */
    for (j = 0; j < *n_fixed; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_random; j++)
      for (k = j; k < *n_random; k++)
	PsiStore[iPsi++] = Psi[j][k];
    for (j = 0; j < *n_grp; j++)
      for (k = 0; k < *n_random; k++)
	gammaStore[igamma++] = gamma[j][k];

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  FreeMatrix(X, *n_samp+*n_fixed);
  free(gamma0);
  FreeMatrix(gamma, *n_grp);
  FreeMatrix(Psi, *n_random);
  FreeMatrix(A0, *n_fixed);
  FreeMatrix(T0, *n_random);
  FreeMatrix(mtemp, *n_fixed);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp + *n_random);
}



void R2logitMixedMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		       double *dX,    /* (N x K) covariate matrix for
					 fixed effects */
		       double *dZ,    /* covariates for random effects
					 organized by groups */
		       int *grp,      /* group indicator, 0, 1, ...,
					   G-1 */
		       double *beta,  /* (K(J-1)) stacked coefficient
					   vector for fixed effects */
		       double *dPsi,  /* LxL precision matrix for
					 random effecs for each equation */
		       int *n_samp,       /* # of obs */
		       int *n_dim,        /* # of categories, J-1 */
		       int *n_fixed,      /* # of fixed effects, K */
		       int *n_random,     /* # of random effects, L */
		       int *n_grp,        /* # of groups, G */
		       int *max_samp_grp, /* max # of obs within each
					     group */
		       double *beta0,    /* (K(J-1)) prior mean vector */
		       double *dA0,      /* (K(J-1) x K(J-1)) prior
					    precision */
		       int *tau0,        /* prior df for Psi */
		       double *dT0,      /* prior scale for Psi */
		       double *tune_fixed,  /* K(J-1) proposal variances */
		       double *tune_random, /* tuning constant for random
					       effects of each random effect */
		       int *n_gen,        /* # of MCMC draws */
		       int *acc_fixed,    /* # of acceptance for fixed effects */
		       int *acc_random,   /* # of acceptance for random
					     effects */
		       double *betaStore,
		       double *PsiStore
		       ) {

   /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int ibeta = 0, iPsi =0;

  /* matrices */
  double *gamma0 = doubleArray(*n_random);
  double **X = doubleMatrix(*n_samp, *n_fixed);
  double ***gamma = doubleMatrix3D(*n_dim, *n_grp, *n_random);
  double ***Psi = doubleMatrix3D(*n_dim, *n_random, *n_random);
  double **A0 = doubleMatrix(n_fixed[0]*n_dim[0], n_fixed[0]*n_dim[0]);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_random);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_random; j++)
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++) {
      for(i = 0; i < *n_dim; i++)
	Psi[i][j][k] = dPsi[itemp];
      itemp++;
    }

  itemp = 0;
  for (j = 0; j < *n_random; j++)
    gamma0[j] = 0;
  for (i = 0; i < *n_dim; i++)
    for (j = 0; j < *n_grp; j++)
      rMVN(gamma[i][j], gamma0, Psi[i], *n_random);

  itemp = 0; 
  for (k = 0; k < n_fixed[0]*n_dim[0]; k++)
    for (j = 0; j < n_fixed[0]*n_dim[0]; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  for (j = 0; j < n_fixed[0]*n_dim[0]; j++) 
    acc_fixed[j] = 0;
  for (j = 0; j < n_grp[0]*n_dim[0]; j++) 
    acc_random[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    logitMixedMetro(Y, X, Zgrp, grp, beta, gamma, Psi, 
		    *n_samp, *n_dim, *n_fixed, *n_random, *n_grp,
		    beta0, A0, *tau0, T0, tune_fixed, tune_random,
		    1, acc_fixed, acc_random);

    R_FlushConsole(); 
    /* Storing the output */
    for (j = 0; j < n_fixed[0]*n_dim[0]; j++)
      betaStore[ibeta++] = beta[j];
    for (i = 0; i < *n_dim; i++)
      for (j = 0; j < *n_random; j++)
	for (k = j; k < *n_random; k++)
	  PsiStore[iPsi++] = Psi[i][j][k];

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  free(gamma0);
  FreeMatrix(X, *n_samp);
  Free3DMatrix(gamma, *n_dim, *n_grp);
  Free3DMatrix(Psi, *n_dim, *n_random);
  FreeMatrix(A0, n_fixed[0]*n_dim[0]);
  FreeMatrix(T0, *n_random);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp);
} 



void R2boprobitMixedMCMC(int *Y,           /* binary outcome variable */
			 double *dX,       /* model matrix for fixed
					      effects */
			 double *dZ,       /* model matrix for random
					      effects */
			 int *grp,         /* group indicator: 0, 1, 2,... */
			 double *beta,     /* fixed effects
					      coefficients */
			 double *tau,      /* cutpoints */
			 double *dPsi,     /* covariance for random
					      effects */
			 int *n_samp,      /* # of obs */ 
			 int *n_fixed,     /* # of fixed effects */
			 int *n_random,    /* # of random effects */
			 int *n_grp,       /* # of groups */
			 int *max_samp_grp, /* max # of obs within
					       group */
			 int *n_cat,       /* # of categories */
			 double *beta0,    /* prior mean */
			 double *dA0,      /* prior precision */
			 int *tau0,        /* prior df */
			 double *dT0,      /* prior scale */
			 int *mh,          /* Metropolis-Hastings step
					      for taus? */
			 double *prop,     /* proposal variance */
			 int *accept,      /* acceptance counter */
			 int *n_gen,       /* # of gibbs draws */
			 /* storage of MCMC draws */
			 double *betaStore, 
			 double *gammaStore,
			 double *tauStore,
			 double *PsiStore
		       ) {

  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int ibeta = 0, igamma = 0, iPsi =0, itau = 0;
  double dtemp;

  /* matrices */
  double **X = doubleMatrix(*n_samp+*n_fixed, *n_fixed+1);
  double **gamma = doubleMatrix(*n_grp, *n_random);
  double *gamma0 = doubleArray(*n_random);
  double **Psi = doubleMatrix(*n_random, *n_random);
  double **A0 = doubleMatrix(*n_fixed, *n_fixed);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double **mtemp = doubleMatrix(*n_fixed, *n_fixed);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp + *n_random, *n_random+1);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_random; j++)
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      Psi[j][k] = dPsi[itemp++];
  
  for (j = 0; j < *n_random; j++)
    gamma0[j] = 0;
  for (j = 0; j < *n_grp; j++)
    rMVN(gamma[j], gamma0, Psi, *n_random);

  itemp = 0; 
  for (k = 0; k < *n_fixed; k++)
    for (j = 0; j < *n_fixed; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  /* adding prior as an additional data point */
  dcholdc(A0, *n_fixed, mtemp);
  for (i = 0; i < *n_fixed; i++) {
    X[*n_samp+i][*n_fixed]=0;
    for (j = 0; j < *n_fixed; j++) {
      X[*n_samp+i][*n_fixed] += mtemp[i][j]*beta0[j];
      X[*n_samp+i][j] = mtemp[i][j];
    }
  }

  if (*mh) {
    for (j = 0; j < *n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < *n_samp; i++){
      dtemp = 0;
      for (j = 0; j < *n_fixed; j++)
        dtemp += X[i][j]*beta[j];
      for (j = 0; j < *n_random; j++)
        dtemp += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
      vitemp[grp[i]]++;
      if (Y[i] == 0)
        X[i][*n_fixed] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else
        X[i][*n_fixed] = TruncNorm(tau[Y[i]-1],tau[Y[i]],dtemp,1,0);
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    boprobitMixedMCMC(Y, X, Zgrp, grp, beta, gamma, tau, Psi, 
		      *n_samp, *n_cat, *n_fixed, *n_random, 
		      *n_grp, 0, beta0, A0, *tau0, T0, *mh, prop,
		      accept, 1);
    
    /* Storing the output */
    for (j = 0; j < *n_fixed; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_random; j++)
      for (k = j; k < *n_random; k++)
	PsiStore[iPsi++] = Psi[j][k];
    for (j = 0; j < *n_grp; j++)
      for (k = 0; k < *n_random; k++)
	gammaStore[igamma++] = gamma[j][k];
    for (j = 0; j < (*n_cat - 1); j++)
      tauStore[itau++] = tau[j];
    
    Rprintf("Acceptance ratio: %14g\n", (double)accept[0]/(double) main_loop);
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  FreeMatrix(X, *n_samp+*n_fixed);
  free(gamma0);
  FreeMatrix(gamma, *n_grp);
  FreeMatrix(Psi, *n_random);
  FreeMatrix(A0, *n_fixed);
  FreeMatrix(T0, *n_random);
  FreeMatrix(mtemp, *n_fixed);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp + *n_random);
}


void R2bNegBin(int *Y,          /* count outcome variable */
	       double *dX,      /* model matrix */
	       double *beta,    /* fixed effects coefficients */
	       double *sig2,    /* dispersion parameter */
	       int *n_samp,     /* # of obs */ 
	       int *n_cov,      /* # of covariates */
	       int *n_gen,      /* # of gibbs draws */
	       double *beta0,   /* prior mean for normal */
	       double *dA0,     /* prior precision for normal; can be
				   set to zero to induce improper prior
				   for beta alone
				*/
	       double *a0,      /* prior shape for sig2 */
	       double *b0,      /* prior scale for sig2 */
	       double *varb,    /* proposal variance for beta */
	       double *vars,    /* proposal variance for sig2 */
	       double *betaStore, 
	       double *sig2Store,
	       int *counter
	       ) {

  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int ibeta = 0, isig2 = 0;

  /* matrices */
  double *cont = doubleArray(*n_samp);
  double **X = doubleMatrix(*n_samp, *n_cov);
  double **A0 = doubleMatrix(*n_cov, *n_cov);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];
  
  /* contrast */
  for (i = 0; i < *n_samp; i++)
    cont[i] = 0;

  /* packing the prior */
  itemp = 0; 
  for (k = 0; k < *n_cov; k++)
    for (j = 0; j < *n_cov; j++)
      A0[j][k] = dA0[itemp++];

  /* Gibbs Sampler! */
  counter[0] = 0; counter[1] = 0;
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    Rprintf("%5d done\n", main_loop);
    negbinMetro(Y, X, beta, sig2, *n_samp, *n_cov, beta0, A0,
		*a0, *b0, varb, *vars, cont, 1, counter, 0);

    /* Storing the output */
    for (j = 0; j < *n_cov; j++)
      betaStore[ibeta++] = beta[j];
    sig2Store[isig2++] = *sig2;

    Rprintf("acceptance ratios:%14g%14g\n",
	    (double)counter[0]/(double)main_loop,
	    (double)counter[1]/(double)main_loop);
	    
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(cont);
  FreeMatrix(X, *n_samp);
  FreeMatrix(A0, *n_cov);
} /* end of negative binomial regression */




void R2bnegbinMixedMCMC(int *Y,           /* outcome variable */
			double *dX,       /* model matrix for fixed
					     effects */
			double *dZ,       /* model matrix for random
					     effects */
			int *grp,         /* group indicator: 0, 1, 2,... */
			double *beta,     /* fixed effects coefficients */
			double *dgamma,   /* random effects coefficients */
			double *sig2,     /* dispersion parameter */
			double *dPsi,     /* covariance for random
					     effects */
			int *n_samp,      /* # of obs */ 
			int *n_fixed,     /* # of fixed effects */
			int *n_random,    /* # of random effects */
			int *n_grp,       /* # of groups */
			int *max_samp_grp,/* max # of obs within group */
			double *beta0,    /* prior mean */
			double *dA0,      /* prior precision */
			double *a0,        /* prior shape for sig2 */
			double *b0,        /* prior scale for sig2 */
			int *tau0,        /* prior df for Psi */
			double *dT0,      /* prior scale for Psi */
			double *varb,     /* proposal variance for
					     beta */
			double *varg,     /* proposal variance for
					     gamma */
			double *vars,     /* proposal variance for
					     sig2 */
			int *n_gen,       /* # of gibbs draws */
			/* counters */
			int *counter,     /* counter for beta, sig2 */
			int *icounterg,    /* counter for gamma */
			/* storage of MCMC draws */
			double *betaStore, 
			double *gammaStore,
			double *sig2Store,
			double *PsiStore
			) {
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int **counterg = intMatrix(*n_grp, 2);
  int ibeta = 0, igamma = 0, iPsi =0, isig2 = 0;

  /* matrices */
  int **Ygrp = intMatrix(*n_grp, *max_samp_grp);
  double **X = doubleMatrix(*n_samp, *n_fixed);
  double **gamma = doubleMatrix(*n_grp, *n_random);
  double **Psi = doubleMatrix(*n_random, *n_random);
  double **A0 = doubleMatrix(*n_fixed, *n_fixed);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_random);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    Ygrp[grp[i]][vitemp[grp[i]]] = Y[i];
    for (j = 0; j < *n_random; j++)
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      Psi[j][k] = dPsi[itemp++];
  
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_grp; j++)
      gamma[j][k] = dgamma[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_fixed; k++)
    for (j = 0; j < *n_fixed; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    counterg[j][0] = icounterg[itemp++];

  counter[0] = 0; counter[1] = 0;
  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {
    bnegbinMixedMCMC(Y, Ygrp, X, Zgrp, grp, beta, gamma, sig2, Psi, 
		     *n_samp, *n_fixed, *n_random, *n_grp,
		     *max_samp_grp, beta0, A0, *a0, *b0, *tau0, T0,
		     varb, *vars, varg, counter, counterg, 1);

    /* Storing the output */
    for (j = 0; j < *n_fixed; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_random; j++)
      for (k = j; k < *n_random; k++)
	PsiStore[iPsi++] = Psi[j][k];
    for (j = 0; j < *n_grp; j++)
      for (k = 0; k < *n_random; k++)
	gammaStore[igamma++] = gamma[j][k];
    sig2Store[isig2++] = sig2[0];

    Rprintf("%14g%14g%14g\n", (double)counter[0]/(double)main_loop,
	    (double)counter[1]/(double)main_loop,
	    (double)counterg[0][0]/(double)main_loop);

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    icounterg[itemp++] = counterg[j][0];

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  FreeintMatrix(counterg, *n_grp);
  FreeintMatrix(Ygrp, *n_grp);
  FreeMatrix(X, *n_samp);
  FreeMatrix(gamma, *n_grp);
  FreeMatrix(Psi, *n_random);
  FreeMatrix(A0, *n_fixed);
  FreeMatrix(T0, *n_random);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp);
} /* end of negative binomial mixed effects model */
