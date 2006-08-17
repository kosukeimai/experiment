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
     Bayesian Normal Regression: see Chap.14 of Gelman et al. (2004) 
       both proper and improper priors (and their combinations)
       allowed for beta and sig2. 
***/
void bNormalReg(double **D,    /* data [X Y] */
		double *beta,  /* coefficients */
		double *sig2,  /* variance */
		int n_samp,    /* sample size */
		int n_cov,     /* # of covariates */
		int addprior,  /* Should prior on beta be incorporated
				  into D? */
		int pbeta,     /* Is prior proper for beta? */
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
  if (addprior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      D[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) {
	D[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
	D[n_samp+i][j] = mtemp[i][j];
      }
    }
  } 
  
  /* SS matrix */
  for(j = 0; j <= n_cov; j++)
    for(k = 0; k <= n_cov; k++)
      SS[j][k]=0;
  for(i = 0;i < n_samp + n_cov; i++)
    for(j = 0;j <= n_cov; j++)
      for(k = 0; k <= n_cov; k++) 
	SS[j][k] += D[i][j]*D[i][k];
  
  /* SWEEP SS matrix */
  for(j = 0; j < n_cov; j++)
    SWP(SS, j, n_cov+1);

  /* draw sig2 from its marginal dist */
  for(j = 0; j < n_cov; j++)
    mean[j] = SS[j][n_cov];
  if (!sig2fixed)
    if (psig2)   /* proper prior for sig2 */
      if (pbeta) /* proper prior for beta */
	sig2[0]=(SS[n_cov][n_cov]+nu0*s0)/rchisq((double)n_samp+nu0);
      else       /* improper prior for beta */
	sig2[0]=(n_samp*SS[n_cov][n_cov]/(n_samp-n_cov)+nu0*s0)/rchisq((double)n_samp+nu0);
    else         /* improper prior for sig2 */
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
   A Gibbs Sampler for Ordinal Probit Regression With and Without
   Marginal Data Augmentation
   
   Marginal Data Augmentation for updating coefficients: 
   p.318 of Imai and van Dyk (2005) Journal of Econometrics.
   Prior mean for beta will be set to zero. 
   Improper prior allowed (set A0 to be a matrix of zeros).
   
   Metropolis-Hasting Blocking Step for updating cutpoints:
   Cowles (1996). Statistics and Computing
   
***/ 

void boprobitMCMC(int *Y,        /* ordinal outcome variable: 0, 1,
				    dots, J-1 */
		  double **X,    /* covariate matrix */
		  double *beta,  /* coefficients */
		  double *tau,   /* J cut points: the first
				    cutpoint is set to 0 and the last
				    cutpoint is set to tau_{J-1}+1000 */
		  int n_samp,    /* # of obs */ 
		  int n_cov,     /* # of covariates */
		  int n_cat,     /* # of categories: J */
		  int prior,     /* Should prior be included in X? */
		  double *beta0, /* prior mean */
		  double **A0,   /* prior precision */
		  int mda,       /* use marginal data augmentation? */
		  int mh,        /* use metropolis-hasting step? */
		  double *prop,  /* J-2 proposal variances for MH step */
		  int *accept,   /* counter for acceptance */
		  int n_gen      /* # of gibbs draws */
		  ) {
  
  /* model parameters */
  double **SS = doubleMatrix(n_cov+1, n_cov+1); /* matrix folders for SWEEP */
  double *mean = doubleArray(n_samp);           /* means for each obs */
  double *mbeta = doubleArray(n_cov);           /* means for beta */
  double **V = doubleMatrix(n_cov, n_cov);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double *Wmax = doubleArray(n_cat);  /* max of W in each categry: 0, 1,
					 ..., J-1 */
  double *Wmin = doubleArray(n_cat);  /* min of W in each category: 0, 1, 
					 ..., J-1 */
  
  /* storage parameters and loop counters */
  int i, j, k, main_loop;  
  double dtemp, dtemp1;
  double *dvtemp = doubleArray(n_cat);
  dvtemp[0] = tau[0];
  double **mtemp = doubleMatrix(n_cov, n_cov);
  
  /* marginal data augmentation */
  double sig2 = 1;
  int nu0 = 1;
  double s0 = 1;

  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_cov, mtemp);
    for(i = 0; i < n_cov; i++) {
      X[n_samp+i][n_cov] = 0;
      for(j = 0; j < n_cov; j++) 
	X[n_samp+i][n_cov] += mtemp[i][j]*beta0[j];
    }
  }

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    for (i = 0; i < n_samp; i++){
      mean[i] = 0;
      for (j = 0; j < n_cov; j++) 
	mean[i] += X[i][j]*beta[j]; 
    }
    /* Sampling tau with MH step */
    if (mh) {
      for (j = 1; j < n_cat-1; j++) 
	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1);
      dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000;
      for (j = 1; j < n_cat-1; j++) 
	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) -
			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) -
	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) -
	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0));
      for (i = 0; i < n_samp; i++) {
	if (Y[i] == (n_cat-1))  
	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-mean[i], 0, 1, 0, 1) -
	    pnorm(tau[n_cat-2]-mean[i], 0, 1, 0, 1);
	else if (Y[i] > 0) 
	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-mean[i], 0, 1, 1, 0) -
			      pnorm(dvtemp[Y[i]-1]-mean[i], 0, 1, 1, 0)) -
	    log(pnorm(tau[Y[i]]-mean[i], 0, 1, 1, 0) -
		pnorm(tau[Y[i]-1]-mean[i], 0, 1, 1, 0));
      }
      if (unif_rand() < exp(dtemp)) {
	accept[0]++;
	for (j = 1; j < n_cat; j++)
	  tau[j] = dvtemp[j];
	if (mda) /* marginal data augmentation */ 
	  sig2 = s0/rchisq((double)nu0);
	for (i = 0; i < n_samp; i++){
	  if (Y[i] == 0) 
	    W[i] = TruncNorm(mean[i]-1000,0,mean[i],1,0);
	  else 
	    W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean[i],1,0);
	}
	X[i][n_cov] = W[i]*sqrt(sig2);
      }
    } else {
      /* Sampling the Latent Variable */
      Wmin[0] = tau[0]; Wmax[0] = tau[0]-10;
      for (j = 1; j < n_cat; j++) {
	Wmin[j] = tau[j];
	Wmax[j] = tau[j-1];
      }
      if (mda) /* marginal data augmentation */ 
	sig2 = s0/rchisq((double)nu0);
      for (i = 0; i < n_samp; i++){
	if (Y[i] == 0) 
	  W[i] = TruncNorm(mean[i]-1000,0,mean[i],1,0);
	else 
	  W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean[i],1,0);
	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]);
	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]);
	X[i][n_cov] = W[i]*sqrt(sig2);
      }
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
      mbeta[j] = SS[j][n_cov];
    for(j = 0; j < n_cov; j++)
      for(k = 0; k < n_cov; k++) V[j][k]=-SS[j][k]*sig2;
    rMVN(beta, mbeta, V, n_cov);
    /* rescaling the parameters */
    if (mda)
      for (j = 0; j < n_cov; j++) beta[j] /= sqrt(sig2);
    
    /* sampling taus without MH-step */
    if (!mh) { 
      for (j = 1; j < n_cat-1; j++) 
	tau[j] = runif(fmax2(tau[j-1], Wmax[j]), 
		       fmin2(tau[j+1], Wmin[j+1]));
      tau[n_cat-1] = tau[n_cat-2] + 1000;
    }
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  
  /* freeing memory */
  FreeMatrix(SS, n_cov+1);
  free(mean);
  free(mbeta);
  FreeMatrix(V, n_cov);
  free(W);
  free(Wmax);
  free(Wmin);
  free(dvtemp);
  FreeMatrix(mtemp, n_cov);
}



/*** 
   A Random Walk Metroplis Sampler for Binomial and Multinomial
   Logistic Regression with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
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
  
  int i, j, k, main_loop;
  double numer, denom;
  double *sumall = doubleArray(n_samp); 
  double *sumall1 = doubleArray(n_samp);
  double *prop = doubleArray(n_dim*n_cov);
  double **Xbeta = doubleMatrix(n_samp, n_dim);
  double **Xbeta1 = doubleMatrix(n_samp, n_dim);

  for (j = 0; j < n_cov*n_dim; j++)
    prop[j] = beta[j];
  for (i = 0; i < n_samp; i++) {
    sumall[i] = 1.0; 
    for (j = 0; j < n_dim; j++) {
      Xbeta[i][j] = 0;
      for (k = 0; k < n_cov; k++) 
	Xbeta[i][j] += X[i][k]*beta[j*n_cov+k];
      Xbeta1[i][j] = Xbeta[i][j];
      sumall[i] += exp(Xbeta[i][j]);
    }
    sumall1[i] = sumall[i];
  }

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
    for (j = 0; j < n_dim; j++)
      for (k = 0; k < n_cov; k++) {
	/** Sample from the proposal distribution **/
	prop[j*n_cov+k] = beta[j*n_cov+k] + 
	  norm_rand()*sqrt(Var[j*n_cov+k]);
      
      /** Calculating the ratio (log scale) **/
      /* prior */
      numer = dMVN(prop, beta0, A0, n_cov*n_dim, 1);
      denom = dMVN(beta, beta0, A0, n_cov*n_dim, 1);   
      /* likelihood */
      for (i = 0; i < n_samp; i++) {
	Xbeta1[i][j] = Xbeta[i][j] - X[i][k]*(beta[j*n_cov+k]-prop[j*n_cov+k]);
	if (Y[i] > 0) {
	  denom += Xbeta[i][Y[i]-1];
	  numer += Xbeta1[i][Y[i]-1];
	} 
	sumall1[i] += (exp(Xbeta1[i][j])-exp(Xbeta[i][j]));
	numer -= log(sumall1[i]);
	denom -= log(sumall[i]);
      }
      
      /** Rejection **/
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	counter[j*n_cov+k]++;
	beta[j*n_cov+k] = prop[j*n_cov+k];
	for (i = 0; i < n_samp; i++) {
	  sumall[i] = sumall1[i];
	  Xbeta[i][j] = Xbeta1[i][j];
	}
      }
    }
  }
  
  free(prop);
  free(sumall);
  free(sumall1);
  FreeMatrix(Xbeta, n_samp);
  FreeMatrix(Xbeta1, n_samp);
}


/*** 
   A Standard Gibbs Sampler for Normal Mixed Effects Regression

   MODEL: Y_i = X_i \beta + Z_i \gamma_i + \epsilon_i 
          where 
          \epsilon_i  \sim N(0, \sigma^2 I_{n_i})
          \gamma_i \sim N(0, \Psi^{-1})
          and i indexes groups.
   PRIOR: p(\beta|X, Z) = N(beta_0, A_0^{-1})
          p(\Psi^{-1}|X, Z) = Wish(\tau_0, T_0)
          p(\sigma^2|X, Z) = InvChi2(\nu_0, s_0)

	  standard diffuse improper priro (imp = 1 && A0 = 0)
	  p(\beta, log\sigma, \Psi|X, Z) \propto |\Psi|^{-1/2}
	  The posterior is proper so long as (see van Dyk, 2000, JCGS) 
             n_samp > n_fixed + n_random 
             n_grp > 2 n_random - 1 
***/ 

void bNormalMixedGibbs(double *Y,       /* outcome variable */
		       double **X,      /* model matrix for fixed
					   effects */
		       double ***Zgrp,  /* model matrix for random
					   effects organized by
					   grous */
		       int *grp,        /* group indicator: 0, 1, 2,... */
		       double *beta,    /* fixed effects coefficients */
		       double **gamma,  /* random effects coefficients */
		       double *sig2,    /* variance parameter */
		       double **Psi,    /* precision matrix for random
					   effects */
		       int n_samp,      /* # of obs */ 
		       int n_fixed,     /* # of fixed effects */
		       int n_random,    /* # of random effects */
		       int n_grp,       /* # of groups */
		       int prior,       /* include prior for fixed effects in X? */
		       double *beta0,   /* prior mean */
		       double **A0,     /* prior precision */
		       int imp,         /* use standard improper prior
					   for sig2 and Psi
					   (see above; beta can be proper) */
		       int nu0,         /* prior df for sig2 */
		       double s0,       /* prior scale for sig2 */
		       int tau0,        /* prior df for Psi */
		       double **T0,     /* prior scale for Psi */
		       int n_gen        /* # of gibbs draws */
		       ) {
  
  double *gamma0 = doubleArray(n_random);           /* prior mean for gamma */
  double **V = doubleMatrix(n_fixed, n_fixed);      /* variances for beta */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  /* storage parameters and loop counters */
  int i, j, k, l, main_loop;  
  int *vitemp = intArray(n_grp);
  double dtemp0, dtemp1;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_fixed, V);
    for(i = 0; i < n_fixed; i++) {
      X[n_samp+i][n_fixed] = 0;
      for(j = 0; j < n_fixed; j++) {
	X[n_samp+i][n_fixed] += V[i][j]*beta0[j];
	X[n_samp+i][j] = V[i][j];
      }
    }
  }

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /** STEP 1: Sample Fixed Effects Given Random Effects 
                Also Sample Variance Parameter **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      X[i][n_fixed] = Y[i];
      for (j = 0; j < n_random; j++)
	X[i][n_fixed] -= Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
      vitemp[grp[i]]++;
    }
    if (imp)
      bNormalReg(X, beta, sig2, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1,
		 1, 0);
    else
      bNormalReg(X, beta, sig2, n_samp, n_fixed, 0, 1, beta0, A0, 1, s0,
		 nu0, 0);

    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      Zgrp[grp[i]][vitemp[grp[i]]][n_random] = Y[i];
      for (j = 0; j < n_fixed; j++) 
	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j]; 
      vitemp[grp[i]]++;
    }
    for (j = 0; j < n_grp; j++)
      bNormalReg(Zgrp[j], gamma[j], sig2, vitemp[j], n_random,
		 1, 1, gamma0, Psi, 0, 0, 1, 1);

    /** STEP 3: Update Covariance Matrix Given Random Effects **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	if (imp)
	  mtemp[j][k] = 0;
	else
	  mtemp[j][k] = T0[j][k];
    for (j = 0; j < n_grp; j++)
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] += gamma[j][k]*gamma[j][l];
    dinv(mtemp, n_random, mtemp1);
    if (imp)
      rWish(Psi, mtemp1, n_grp-n_random-1, n_random);
    else
      rWish(Psi, mtemp1, tau0+n_grp, n_random);

    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(gamma0);
  FreeMatrix(V, n_fixed);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
  free(vitemp);
}



/*** 
   A Standard Gibbs Sampler for Binary Probit Mixed Effects Regression

   MODEL: Y_{ij} = 1 if W_{ij} > 0 
                 = 0 otherwise.
          W_i = X_i \beta + Z_i \gamma_i + \epsilon_i 
          where 
          \epsilon_i  \sim N(0, I_{n_i})
          \gamma_i \sim N(0, \Psi^{-1})
          and i indexes groups.
   PRIOR: p(\beta|X, Z) = N(beta_0, A_0^{-1})
          p(\Psi^{-1}|X,Z) = Wish(\tau_0, T_0)
   see the docs for bprobitGibbs for the implmentation of marginal
          data augmentation for fixed effects coefficients       
***/ 

void bprobitMixedGibbs(int *Y,          /* binary outcome variable */
		       double **X,      /* model matrix for fixed
					   effects */
		       double ***Zgrp,  /* model matrix for random
					   effects organized by grous */
		       int *grp,        /* group indicator: 0, 1, 2,... */
		       double *beta,    /* fixed effects coefficients */
		       double **gamma,  /* random effects coefficients */
		       double **Psi,    /* precision matrix for random
					   effects */
		       int n_samp,      /* # of obs */ 
		       int n_fixed,     /* # of fixed effects */
		       int n_random,    /* # of random effects */
		       int n_grp,       /* # of groups */
		       int prior,       /* include prior for fixed effects in X? */
		       double *beta0,   /* prior mean */
		       double **A0,     /* prior precision */
		       int tau0,        /* prior df */
		       double **T0,     /* prior scale */
		       int n_gen        /* # of gibbs draws */
		       ) {
  
  double *gamma0 = doubleArray(n_random);           /* prior mean for gamma */
  double **V = doubleMatrix(n_fixed, n_fixed);      /* variances for beta */
  double *W = doubleArray(n_samp);
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  /* storage parameters and loop counters */
  int i, j, k, l, main_loop;  
  int *vitemp = intArray(n_grp);
  double dtemp0, dtemp1;
  double *vdtemp = doubleArray(1);
  vdtemp[0] = 1.0;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_fixed, V);
    for(i = 0; i < n_fixed; i++) {
      X[n_samp+i][n_fixed] = 0;
      for(j = 0; j < n_fixed; j++) {
	X[n_samp+i][n_fixed] += V[i][j]*beta0[j];
	X[n_samp+i][j] = V[i][j];
      }
    }
  }

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /** STEP 1: Sample Latent Variable **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++){
      dtemp0 = 0; dtemp1 = 0;
      for (j = 0; j < n_fixed; j++) 
	dtemp0 += X[i][j]*beta[j]; 
      for (j = 0; j < n_random; j++)
	dtemp1 += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
      if(Y[i] == 0) 
	W[i] = TruncNorm(dtemp0+dtemp1-1000,0,dtemp0+dtemp1,1,0);
      else 
	W[i] = TruncNorm(0,dtemp0+dtemp1+1000,dtemp0+dtemp1,1,0);
      X[i][n_fixed] = W[i]-dtemp1;
      vitemp[grp[i]]++;
    }
    /** STEP 2: Sample Fixed Effects Given Random Effects **/
    bNormalReg(X, beta, vdtemp, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1,
	       1, 1);

    /** STEP 3: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      Zgrp[grp[i]][vitemp[grp[i]]][n_random] = W[i];
      for (j = 0; j < n_fixed; j++) 
	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j]; 
      vitemp[grp[i]]++;
    }
    for (j = 0; j < n_grp; j++)
      bNormalReg(Zgrp[j], gamma[j], vdtemp, vitemp[j], n_random,
		 1, 1, gamma0, Psi, 0, 0, 1, 1);

    /** STEP 4: Update Covariance Matrix Given Random Effects **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (j = 0; j < n_grp; j++)
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] += gamma[j][k]*gamma[j][l];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, tau0+n_grp, n_random);

    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(W);
  free(vdtemp);
  free(vitemp);
  free(gamma0);
  FreeMatrix(V, n_fixed);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects probit */



/*** 
   A Random Walk Metroplis Sampler for Binomial and Multinomial
   Logistic Mixed Effects Regression with Independent Normal Prior and
   Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 
***/
void logitMixedMetro(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		     double **X,    /* (N x K) covariate matrix for
				       fixed effects */
		     double ***Z,     /* covariates for random effects
					 organized by groups */
		     int *grp,        /* group indicator, 0, 1, ...,
					 G-1 */
		     double *beta,    /* (K(J-1)) stacked coefficient
					 vector for fixed effects */
		     double ***gamma, /*(J-1)G L array of random
					 effects organized by
					 equations and then by
					 groups */
		     double ***Psi,   /* LxL precision matrix for
					 random effecs for each equation */
		     int n_samp,      /* # of obs */
		     int n_dim,       /* # of categories, J-1 */
		     int n_fixed,     /* # of fixed effects, K */
		     int n_random,    /* # of random effects, L */
		     int n_grp,       /* # of groups, G */
		     double *beta0,   /* (K(J-1)) prior mean vector */
		     double **A0,     /* (K(J-1) x K(J-1)) prior
				       precision */
		     int tau0,        /* prior df for Psi */
		     double **T0,     /* prior scale for Psi */
		     double *tune_fixed,   /* K(J-1) proposal variances */
		     double *tune_random,  /* tuning constant for random
					      effects of each random effect */
		     int n_gen,        /* # of MCMC draws */
		     int *acc_fixed,   /* # of acceptance for fixed effects */
		     int *acc_random   /* # of acceptance for random effects */
		     ) {
  
  int i, j, k, l, main_loop;
  int *vitemp = intArray(n_grp);
  double numer, denom;
  double *sumall = doubleArray(n_samp); 
  double *sumall1 = doubleArray(n_samp);
  double *propb = doubleArray(n_dim*n_fixed);
  double *propg = doubleArray(n_random);
  double *gamma0 = doubleArray(n_random);
  double **Xbeta = doubleMatrix(n_samp, n_dim);
  double **Xbeta1 = doubleMatrix(n_samp, n_dim);
  double **Zgamma = doubleMatrix(n_samp, n_dim);
  double **Zgamma1 = doubleMatrix(n_samp, n_dim);
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed*n_dim; j++)
    propb[j] = beta[j];
  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;
  for (j = 0 ; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    sumall[i] = 1.0; 
    for (j = 0; j < n_dim; j++) {
      Xbeta[i][j] = 0; Zgamma[i][j] = 0;
      for (k = 0; k < n_fixed; k++) 
	Xbeta[i][j] += X[i][k]*beta[j*n_fixed+k];
      Xbeta1[i][j] = Xbeta[i][j];
      for (k = 0; k < n_random; k++)
	Zgamma[i][j] += Z[j][vitemp[grp[i]]][k]*gamma[j][grp[i]][k];
      sumall[i] += exp(Xbeta[i][j] + Zgamma[i][j]);
      Zgamma1[i][j] = Zgamma[i][j];
    }
    sumall1[i] = sumall[i];
    vitemp[grp[i]]++;
  }

  for (main_loop = 0; main_loop < n_gen; main_loop++) {
     /** STEP 1: Update Fixed Effects Given Random Effects **/
    for (j = 0; j < n_dim; j++)
      for (k = 0; k < n_fixed; k++) {
	/** Sample from the proposal distribution **/
	propb[j*n_fixed+k] = beta[j*n_fixed+k] + 
	  norm_rand()*sqrt(tune_fixed[j*n_fixed+k]);
	/** Calculating the ratio (log scale) **/
	/* prior */
	numer = dMVN(propb, beta0, A0, n_fixed*n_dim, 1);
	denom = dMVN(beta, beta0, A0, n_fixed*n_dim, 1);   
	/* likelihood */
	for (i = 0; i < n_samp; i++) {
	  Xbeta1[i][j] = Xbeta[i][j] - X[i][k]*(beta[j*n_fixed+k]-propb[j*n_fixed+k]);
	  if (Y[i] > 0) {
	    denom += (Xbeta[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	    numer += (Xbeta1[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	  } 
	  sumall1[i] += (exp(Xbeta1[i][j] + Zgamma[i][Y[i]-1]) - 
			 exp(Xbeta[i][j] + Zgamma[i][Y[i]-1]));
	  numer -= log(sumall1[i]);
	  denom -= log(sumall[i]);
	}
	/** Rejection **/
	if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	  acc_fixed[j*n_fixed+k]++;
	  beta[j*n_fixed+k] = propb[j*n_fixed+k];
	  for (i = 0; i < n_samp; i++) {
	    sumall[i] = sumall1[i];
	    Xbeta[i][j] = Xbeta1[i][j];
	  }
	}
      }
 
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_dim; j++) {
      dinv(Psi[j], n_random, mtemp);
      for (i = 0; i < n_random; i++)
	for (k = 0; k < n_random; k++)
	  mtemp[i][k] *= tune_random[j];
      for (k = 0; k < n_grp; k++) {
	rMVN(propg, gamma[j][k], mtemp, n_random);
	/** Calculating the ratio (log scale) **/
	/* prior */
	numer = dMVN(propg, gamma0, Psi[j], n_random, 1);
	denom = dMVN(gamma[j][k], gamma0, Psi[j], n_random, 1); 
 	/* likelihood */
	for (l = 0; l < n_grp; l++)
	  vitemp[l] = 0;
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == k)
	    for (l = 0; l < n_random; l++)
	      Zgamma1[i][j] -= Z[k][vitemp[k]][l]*(gamma[j][k][l]-propg[l]);
	  vitemp[grp[i]]++;
	  if (Y[i] > 0) {
	    denom += (Xbeta[i][Y[i]-1] + Zgamma[i][Y[i]-1]);
	    numer += (Xbeta[i][Y[i]-1] + Zgamma1[i][Y[i]-1]);
	  } 
	  sumall1[i] += (exp(Xbeta[i][j] + Zgamma1[i][j]) -
			 exp(Xbeta[i][j] + Zgamma[i][j]));
	  numer -= log(sumall1[i]);
	  denom -= log(sumall[i]);
	}
	/* Rejection */
	if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	  acc_random[j*n_grp+k]++;
	  for (l = 0; l < n_random; l++)
	    gamma[j][k][l] = propb[l];
	  for (i = 0; i < n_samp; i++) {
	    sumall[i] = sumall1[i];
	    Zgamma[i][j] = Zgamma1[i][j];
	  }      
	}
      }
    }
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_dim; j++) {
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] = T0[k][l];
      for (i = 0; i < n_grp; i++)
	for (k = 0; k < n_random; k++)
	  for (l = 0; l < n_random; l++)
	    mtemp[k][l] += gamma[j][i][k]*gamma[j][i][l];
      dinv(mtemp, n_random, mtemp1);
      rWish(Psi[j], mtemp1, tau0+n_grp, n_random);
    }
  }

  /* freeing memory */
  free(vitemp);
  free(sumall);
  free(sumall1);
  free(propb);
  free(propg);
  free(gamma0);
  FreeMatrix(Xbeta, n_samp);
  FreeMatrix(Xbeta1, n_samp);
  FreeMatrix(Zgamma, n_samp);
  FreeMatrix(Zgamma1, n_samp);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */




/*** 
   A Standard Gibbs Sampler for Ordinal Probit Mixed Effects Regression
***/ 

void boprobitMixedMCMC(int *Y,          /* binary outcome variable */
		       double **X,      /* model matrix for fixed
					   effects */
		       double ***Zgrp,  /* model matrix for random
					   effects organized by grous */
		       int *grp,        /* group indicator: 0, 1, 2,... */
		       double *beta,    /* fixed effects coefficients */
		       double **gamma,  /* random effects coefficients */
		       double *tau,     /* cutpoints */
		       double **Psi,    /* precision matrix for random
					   effects */
		       int n_samp,      /* # of obs */ 
		       int n_cat,       /* number of categories */
		       int n_fixed,     /* # of fixed effects */
		       int n_random,    /* # of random effects */
		       int n_grp,       /* # of groups */
		       int prior,       /* include prior for fixed effects in X? */
		       double *beta0,   /* prior mean */
		       double **A0,     /* prior precision */
		       int tau0,        /* prior df */
		       double **T0,     /* prior scale */
		       int mh,          /* metropolis-hastings step
					   for cutpoints? */
		       double *prop,    /* proposal variance for MH
					   step */
		       int *accept,     /* counter for acceptance */
		       int n_gen        /* # of gibbs draws */
		       ) {
  
  double *gamma0 = doubleArray(n_random);         /* prior mean for gamma */
  double *Xbeta = doubleArray(n_samp);            /* X beta */
  double *Zgamma = doubleArray(n_samp);
  double **V = doubleMatrix(n_fixed, n_fixed);    /* variances for beta */
  double *W = doubleArray(n_samp);
  double *Wmax = doubleArray(n_cat);  /* max of W in each categry: 0, 1,
					 ..., J-1 */
  double *Wmin = doubleArray(n_cat);  /* min of W in each category: 0, 1, 
					 ..., J-1 */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);
  double *mean0 = doubleArray(n_samp);
  double *mean1 = doubleArray(n_samp);

  /* storage parameters and loop counters */
  int i, j, k, l, main_loop;  
  int *vitemp = intArray(n_grp);
  double dtemp;
  double *vdtemp = doubleArray(1);
  double *dvtemp = doubleArray(n_cat);
  vdtemp[0] = 1.0;
  
  /* read the prior as additional data points */
  if (prior) {
    dcholdc(A0, n_fixed, V);
    for(i = 0; i < n_fixed; i++) {
      X[n_samp+i][n_fixed] = 0;
      for(j = 0; j < n_fixed; j++) {
	X[n_samp+i][n_fixed] += V[i][j]*beta0[j];
	X[n_samp+i][j] = V[i][j];
      }
    }
  }

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= n_gen; main_loop++){
    /** STEP 1: Sample Latent Variable **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++){
      Xbeta[i] = 0; Zgamma[i] = 0;
      for (j = 0; j < n_fixed; j++) 
	Xbeta[i] += X[i][j]*beta[j]; 
      for (j = 0; j < n_random; j++)
	Zgamma[i] += Zgrp[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
      vitemp[grp[i]]++;
    }
    /* Sampling tau with MH step */
    if (mh) {
      for (j = 1; j < n_cat-1; j++) 
	dvtemp[j] = TruncNorm(dvtemp[j-1], tau[j+1], tau[j], prop[j-1], 1);
      dtemp = 0; dvtemp[n_cat-1] = dvtemp[n_cat-2] + 1000;
      for (j = 1; j < n_cat-1; j++) 
	dtemp = dtemp + log(pnorm(tau[j+1]-tau[j], 0, sqrt(prop[j-1]), 1, 0) -
			    pnorm(dvtemp[j-1]-tau[j], 0, sqrt(prop[j-1]), 1, 0)) -
	  log(pnorm(dvtemp[j+1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0) -
	      pnorm(tau[j-1]-dvtemp[j], 0, sqrt(prop[j-1]), 1, 0));
      for (i = 0; i < n_samp; i++) {
	if (Y[i] == (n_cat-1))  
	  dtemp = dtemp + pnorm(dvtemp[n_cat-2]-mean0[i]-mean1[i], 0, 1, 0, 1) -
	    pnorm(tau[n_cat-2]-mean0[i]-mean1[i], 0, 1, 0, 1);
	else if (Y[i] > 0) 
	  dtemp = dtemp + log(pnorm(dvtemp[Y[i]]-mean0[i]-mean1[i], 0, 1, 1, 0) -
			      pnorm(dvtemp[Y[i]-1]-mean0[i]-mean1[i], 0, 1, 1, 0)) -
	    log(pnorm(tau[Y[i]]-mean0[i]-mean1[i], 0, 1, 1, 0) -
		pnorm(tau[Y[i]-1]-mean0[i]-mean1[i], 0, 1, 1, 0));
      }
      if (unif_rand() < exp(dtemp)) {
	accept[0]++;
	for (j = 1; j < n_cat; j++)
	  tau[j] = dvtemp[j];
	for (i = 0; i < n_samp; i++){
	  if (Y[i] == 0) 
	    W[i] = TruncNorm(mean0[i]+mean1[i]-1000,0,mean0[i]+mean1[i],1,0);
	  else 
	    W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean0[i]+mean1[i],1,0);
	X[i][n_fixed] = W[i]-mean1[i];
	}
      }
    } else {
      /* Sampling the Latent Variable */
      Wmin[0] = tau[0]; Wmax[0] = tau[0]-10;
      for (j = 1; j < n_cat; j++) {
	Wmin[j] = tau[j];
	Wmax[j] = tau[j-1];
      }
      for (i = 0; i < n_samp; i++){
	if (Y[i] == 0) 
	  W[i] = TruncNorm(mean0[i]+mean1[i]-1000,0,mean0[i]+mean1[i],1,0);
	else 
	  W[i] = TruncNorm(tau[Y[i]-1],tau[Y[i]],mean0[i]+mean1[i],1,0);
	Wmax[Y[i]] = fmax2(Wmax[Y[i]], W[i]);
	Wmin[Y[i]] = fmin2(Wmin[Y[i]], W[i]);
	X[i][n_fixed] = W[i];
      }
    }

    /** STEP 2: Sample Fixed Effects Given Random Effects **/
    bNormalReg(X, beta, vdtemp, n_samp, n_fixed, 0, 1, beta0, A0, 0, 1,
	       1, 1);

    /** STEP 3: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++)
      vitemp[j] = 0;
    for (i = 0; i < n_samp; i++) {
      Zgrp[grp[i]][vitemp[grp[i]]][n_random] = W[i];
      for (j = 0; j < n_fixed; j++) 
	Zgrp[grp[i]][vitemp[grp[i]]][n_random] -= X[i][j]*beta[j]; 
      vitemp[grp[i]]++;
    }
    for (j = 0; j < n_grp; j++)
      bNormalReg(Zgrp[j], gamma[j], vdtemp, vitemp[j], n_random,
		 1, 1, gamma0, Psi, 0, 0, 1, 1);

    /** STEP 4: Update Covariance Matrix Given Random Effects **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (j = 0; j < n_grp; j++)
      for (k = 0; k < n_random; k++)
	for (l = 0; l < n_random; l++)
	  mtemp[k][l] += gamma[j][k]*gamma[j][l];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, tau0+n_grp, n_random);

    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  /* freeing memory */
  free(gamma0);
  free(mean0);
  free(mean1);
  free(Xbeta);
  free(Zgamma);
  FreeMatrix(V, n_fixed);
  free(W);
  free(Wmax);
  free(Wmin);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
  free(vdtemp);
  free(dvtemp);
  free(vitemp);
} /* end of mixed effects ordinal probit */


