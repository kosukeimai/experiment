#include <stddef.h>
#include <stdio.h>      
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "models.h"
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

/*** 
   A Random Walk Metroplis Sampler for Binomial and Multinomial
   Logistic Regression with Independent Normal Prior 
     Proposal distribution is the multivariate normal whose mean is
   the current value and variance is given by the input.
***/
void logitMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		double **X,    /* (N x K) covariate matrix */
		double *beta,  /* (K(J-1)) stacked coefficient vector */
		int n_samp,    /* # of obs */
		int n_dim,     /* # of categories, J-1 */
		int n_cov,     /* # of covariates, K */
		double *beta0, /* (K(J-1)) prior mean vector */
		double **A0,   /* (K(J-1) x K(J-1)) prior precision */
		double *Var,   /* K(J-1) proposal variances */
		int n_gen,     /* # of MCMC draws */
		int *counter   /* # of acceptance for each parameter */
		) {
  
  int i, j, k, main_loop, param;
  double numer, denom;
  double sumall, sumall1, dtemp, dtemp1;
  double *prop = doubleArray(n_dim*n_cov);

  for (j = 0; j < n_cov*n_dim; j++)
    prop[j] = beta[j];

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    for (param = 0; param < n_dim*n_cov; param++) {
      /** Sample from the proposal distribution **/
      prop[param] = beta[param] + norm_rand()*sqrt(Var[param]);
      
      /** Calculating the ratio (log scale) **/
      /* prior */
      numer = dMVN(prop, beta0, A0, n_cov*n_dim, 1);
      denom = dMVN(beta, beta0, A0, n_cov*n_dim, 1);   
      /* likelihood */
      for (i = 0; i < n_samp; i++) {
	sumall = 1.0; sumall1 = 1.0;
	for (j = 0; j < n_dim; j++) {
	  dtemp = 0; dtemp1 = 0;
	  for (k = 0; k < n_cov; k++) {
	    dtemp += X[i][k]*beta[j*n_cov+k];
	    dtemp1 += X[i][k]*prop[j*n_cov+k];
	  }
	  if (Y[i] == (j+1)) {
	    denom += dtemp;
	    numer += dtemp1;
	  } 
	  sumall += exp(dtemp);
	  sumall1 += exp(dtemp1);
	}
	numer -= log(sumall1);
	denom -= log(sumall);
      }
      
      /** Rejection **/
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	counter[param]++;
	beta[param] = prop[param];
      }
    }
  }

  free(prop);
}

/*** 
     Bayesian Normal Regression: see Chap.14 of Gelman et al. (2004) 
       both proper and improper priors (and their combinations)
       allowed for beta and sig2. 
***/
void bNormalReg(double *Y,     /* response variable */
		double **X,    /* model matrix */
		double *beta,  /* coefficients */
		double *sig2,  /* variance */
		int n_samp,    /* sample size */
		int n_cov,     /* # of covariates */
		int pbeta,     /* 0: improper prior 
				  p(beta|X) \propto 1
				  1: proper prior for (beta)
				  p(beta|X) = normal(beta0, A0)
			       */
		double *beta0, /* prior mean for normal */
		double **A0,   /* prior precision for normal; can be
				  set to zero to induce improper prior
				  for beta alone
			       */
		int psig2,     /* 0: improper prior for sig2
				  p(sig2|X) \propto 1/sig2
				  1: proper prior for sig2
				  p(sigma2|X) = InvChi2(nu0, s0)
			       */
		double s0,     /* prior scale for InvChi2 */
		int nu0,       /* prior d.f. for InvChi2 */
		int sig2fixed  /* 1: sig2 fixed, 0: sig2 sampled */ 
		) {
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /* storage parameters and loop counters */
  int i, j, k;  

  /* read the proper prior for beta as additional data points */
  if (pbeta) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }
  for (i = 0; i < n_samp; i++)
    X[i][n_cov] = Y[i];

  /* SS matrix */
  for(j = 0; j <= n_cov; j++)
    for(k = 0; k <= n_cov; k++)
      SS[j][k]=0;
  for(i = 0;i < n_samp; i++)
    for(j = 0;j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++) 
	SS[j][k] += X[i][j]*X[i][k];
  if (pbeta) 
    for(i = n_samp;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
  
  /* SWEEP SS matrix */
  for(j = 0; j < n_cov; j++)
    SWP(SS, j, n_cov+1);

  /* draw sig2 from its marginal dist */
  for(j = 0; j < n_cov; j++)
    mean[j] = SS[j][n_cov];
  if (!sig2fixed)
    if (psig2)
      if (pbeta)
	sig2[0]=(SS[n_cov][n_cov]+nu0*s0)/rchisq((double)n_samp+nu0);
      else
	sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0);
    else
      sig2[0]=SS[n_cov][n_cov]/rchisq((double)n_samp-n_cov);
  
  /* draw beta from its conditional given sig2 */
  for(j = 0; j < n_cov; j++)
    for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2[0];
  rMVN(beta, mean, V, n_cov);

  /* freeing memory */
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/*** 
   A Gibbs Sampler for Binary Probit Regression With and Without
   Marginal Data Augmentation
   
   Marginal Data Augmentation: see p.318 of Imai and van Dyk (2005)
   Journal of Econometrics.
      Prior mean for beta will be set to zero. 
      Improper prior allowed (set A0 to be a matrix of zeros).
***/ 

void bprobitGibbs(int *Y,        /* binary outcome variable */
		  double **X,    /* covariate matrix */
		  double *beta,  /* coefficients */
		  int n_samp,    /* # of obs */ 
		  int n_cov,     /* # of covariates */
		  int prior,     /* Should prior be included in X? */
		  double *beta0, /* prior mean */
		  double **A0,   /* prior precision */
		  int mda,       /* Want to use marginal data augmentation? */ 
		  int n_gen      /* # of gibbs draws */
		  ) {
  
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_cov);            /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double **mtemp = doubleMatrix(n_cov, n_cov);

  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp;
  
  /* marginal data augmentation */
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	if (!mda)
	  X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	X[n_samp+i][j] = mtemp[i][j];
      }
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /* marginal data augmentation */
    if (mda) sig2 = s0/rchisq((double)nu0);
    
    for (i = 0; i < n_samp; i++){
      dtemp = 0;
      for (j = 0; j < n_cov; j++) 
	dtemp += X[i][j]*beta[j]; 
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp-1000,0,dtemp,1,0);
      else 
	W[i] = TruncNorm(0,dtemp+1000,dtemp,1,0);
      X[i][n_cov] = W[i]*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    for(i = n_samp;i < n_samp+n_cov; i++)
      for(j = 0;j <= n_cov; j++)
	for(k = 0; k <= n_cov; k++) 
	  SS[j][k] += X[i][j]*X[i][k];

    /* SWEEP SS matrix */
    for(j = 0; j < n_cov; j++)
      SWP(SS, j, n_cov+1);

    /* draw beta */    
    for(j = 0; j < n_cov; j++)
      mean[j] = SS[j][n_cov];
    if (mda) 
      sig2=(SS[n_cov][n_cov]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mean, V, n_cov);
 
    /* rescaling the parameters */
    if(mda) 
      for (j = 0; j < n_cov; j++) beta[j] /= sqrt(sig2);
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(W);
  free(mean);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(V, n_cov);
  FreeMatrix(mtemp, n_cov);
}


/*** 
   A Standard Gibbs Sampler for Binary Probit Mixed Effects Regression

   MODEL: Y_i = X_i \beta + Z_i \gamma_i + \epsilon_i where i indexes groups.
          \epsilon_i  \sim N(0, I_{n_i})
          \gamma_i \sim N(0, \Psi)
   PRIOR: p(\beta|X, Z) = N(beta_0, A_0^{-1})
          p(\Psi|X,Z) = InvWish(\tau_0, T_0)
   see the docs for bprobitGibbs for the implmentation of marginal
          data augmentation for fixed effects coefficients       
***/ 

void bprobitMixedGibbs(int *Y,          /* binary outcome variable */
		       double **X,      /* model matrix for fixed
					   effects */
		       double **Z,      /* model matrix for random
					   effects */
		       double ***Zgrp,  /* model matrix organized by
					   grous */
		       int *grp,        /* group indicator: 0, 1, 2,... */
		       double *beta,    /* fixed effects coefficients */
		       double **gamma,  /* random effects coefficients */
		       double **Psi,    /* covariance for random
					   effects */
		       int n_samp,      /* # of obs */ 
		       int n_fixed,     /* # of fixed effects */
		       int n_random,    /* # of random effects */
		       int n_grp,       /* # of groups */
		       int *n_samp_grp, /* # of obs within group */
		       int max_samp_grp, /* max # of obs within group */
		       int prior,       /* include prior in X? */
		       double *beta0,   /* prior mean */
		       double **A0,     /* prior precision */
		       int tau0,         /* prior df */
		       double **T0,     /* prior scale */
		       int mda,         /* Want to use marginal data
					   augmentation for fixed effects ? */ 
		       int n_gen        /* # of gibbs draws */
		       ) {
  
  double *gamma0 = doubleArray(n_random);           /* prior mean for gamma */
  double **SS = doubleMatrix(n_fixed+1, n_fixed+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_fixed);              /* means for beta */
  double **V = doubleMatrix(n_fixed, n_fixed);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double **Wgrp = doubleMatrix(n_grp, max_samp_grp);
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  /* storage parameters and loop counters */
  int i, j, k, l, main_loop;  
  int *vitemp = intArray(n_grp);
  double dtemp0, dtemp1;
  double *vdtemp = doubleArray(1);
  vdtemp[0] = 1.0;

  /* marginal data augmentation */
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_fixed, V);
    for(i = 0; i < n_fixed; i++) {
      X[n_samp+i][n_fixed] = 0;
      for(j = 0; j < n_fixed; j++) {
	if (!mda)
	  X[n_samp+i][n_fixed] += V[i][j]*beta0[j];
	X[n_samp+i][j] = V[i][j];
      }
    }
  }

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /** STEP 1: Sample Fixed Effects Given Random Effects **/
    /* marginal data augmentation */
    if (mda) sig2 = s0/rchisq((double)nu0);

    for (i = 0; i < n_samp; i++){
      dtemp0 = 0; dtemp1 = 0;
      for (j = 0; j < n_fixed; j++) 
	dtemp0 += X[i][j]*beta[j]; 
      for (j = 0; j < n_random; j++)
	dtemp1 += Z[i][j]*gamma[grp[i]][j];
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp0+dtemp1-1000,0,dtemp0+dtemp1,1,0);
      else 
	W[i] = TruncNorm(0,dtemp0+dtemp1+1000,dtemp0+dtemp1,1,0);
      X[i][n_fixed] = (W[i]-dtemp1)*sqrt(sig2);
      W[i] *= sqrt(sig2);
    }

    /* SS matrix */
    for(j = 0; j <= n_fixed; j++)
      for(k = 0; k <= n_fixed; k++)
	SS[j][k]=0;
    for(i = 0;i < n_samp; i++)
      for(j = 0;j <= n_fixed; j++)
	for(k = 0; k <= n_fixed; k++) 
	  SS[j][k] += X[i][j]*X[i][k];
    for(i = n_samp;i < n_samp+n_fixed; i++)
      for(j = 0;j <= n_fixed; j++)
	for(k = 0; k <= n_fixed; k++) 
	  SS[j][k] += X[i][j]*X[i][k];

    /* SWEEP SS matrix */
    for(j = 0; j < n_fixed; j++)
      SWP(SS, j, n_fixed+1);

    /* draw beta */    
    for(j = 0; j < n_fixed; j++)
      mean[j] = SS[j][n_fixed];
    if (mda) 
      sig2=(SS[n_fixed][n_fixed]+s0)/rchisq((double)n_samp+nu0);
    for(j = 0; j < n_fixed; j++)
      for(k = 0; k < n_fixed; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mean, V, n_fixed);

    /* rescaling the parameters */
    if(mda) {
      for (j = 0; j < n_fixed; j++) beta[j] /= sqrt(sig2);
      for (i = 0; i < n_samp; i++) W[i] /= sqrt(sig2);
    }

    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      Wgrp[grp[i]][vitemp[grp[i]]] = W[i];
      for (j = 0; j < n_fixed; j++) 
	Wgrp[grp[i]][vitemp[grp[i]]] -= X[i][j]*beta[j]; 
      vitemp[grp[i]]++;
    }
    dinv(Psi, n_random, mtemp);
    for (j = 0; j < n_grp; j++) {
      bNormalReg(Wgrp[j], Zgrp[j], gamma[j], vdtemp, n_samp_grp[j], n_random,
		 1, gamma0, mtemp, 0, 0, 1, 1);
    }

    /** STEP 3: Update Covariance Matrix Given Random Effects **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (j = 0; j < n_grp; j++)
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] += gamma[j][k]*gamma[j][l];
    dinv(mtemp, n_random, mtemp1);
    rWish(mtemp, mtemp1, tau0+n_grp, n_random);
    dinv(mtemp, n_random, Psi);

    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  PdoubleArray(gamma0, n_random);

  /* freeing memory */
  free(W);
  free(mean);
  free(vdtemp);
  free(vitemp);
  free(gamma0);
  FreeMatrix(SS, n_fixed+1);
  FreeMatrix(V, n_fixed);
  FreeMatrix(Wgrp, n_grp);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
}

