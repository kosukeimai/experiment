/* normal regression */
void bNormalReg(double **D, double *beta, double *sig2, 
		int n_samp, int n_cov, int addprior, int pbeta, 
		double *beta0, double **A0, int psig2, double s0, 
		int nu0, int sig2fixed);
  
/* binomial probit regression */
void bprobitGibbs(int *Y, double **X, double *beta, int n_samp, 
		  int n_cov, int prior, double *beta0, double **A0, 
		  int mda, int n_gen);

/* ordinal probit regression */
void boprobitMCMC(int *Y, double **X, double *beta, 
		  double *tau, int n_samp, int n_cov, int n_cat, 
		  int prior, double *beta0, double **A0, int mda, 
		  int mh, double *prop, int *accept, int n_gen);

/* binomial and mulitnomial logistic regression */
void logitMetro(int *Y, double **X, double *beta, int n_samp,      
		int n_dim, int n_cov, double *beta0, double **A0,     
		double *Var, int n_gen, int *counter);

/* Normal mixed effects regression */
void bNormalMixedGibbs(double *Y, double **X, double ***Zgrp,
		       int *grp, double *beta, double **gamma, double *sig2,    
		       double **Psi, int n_samp, int n_fixed, int n_random,
		       int n_grp, int prior, double *beta0, double **A0, 
		       int imp, int nu0, double s0, int tau0, double **T0, 
		       int n_gen0); 

/* binomial mixed effects probit regression */
void bprobitMixedGibbs(int *Y, double **X, double ***Zgrp, 
		       int *grp, double *beta, double **gamma, 
		       double **Psi, int n_samp, int n_fixed, 
		       int n_random, int n_grp, 
		       int prior, double *beta0, double **A0, 
		       int tau0, double **T0, int n_gen);

/* (binomial/multinomial) logistic mixed effects regression */
void logitMixedMetro(int *Y, double **X, double ***Z, int *grp,
		     double *beta, double ***gamma, double ***Psi,
		     int n_samp, int n_dim, int n_fixed,
		     int n_random, int n_grp, double *beta0,
		     double **A0, int tau0, double **T0,
		     double *tune_fixed, double *tune_random,
		     int n_gen, int *acc_fixed, int *acc_random);

/* ordinal probit mixed effects regression */
void boprobitMixedMCMC(int *Y, double **X, double ***Zgrp, int *grp,
		       double *beta, double **gamma, double *tau,
		       double **Psi, int n_samp, int n_cat,
		       int n_fixed, int n_random, int n_grp,
		       int prior, double *beta0, double **A0, int tau0,
		       double **T0, int mh, double *prop, int *accept,
		       int n_gen);

/* negative binomial regression */
void negbinMetro(int *Y, double **X, double *beta, double *sig2,
		 int n_samp, int n_cov, double *beta0, double **A0, 
		 double a0, double b0, double *varb, double vars,
		 double *cont, int n_gen, int *counter, int sig2fixed);

/* mixed effects negative binomial regression */
void bnegbinMixedMCMC(int *Y, int **Ygrp, double **X, double ***Zgrp,
		      int *grp, double *beta, double **gamma,
		      double *sig2, double **Psi, int n_samp,
		      int n_fixed, int n_random, int n_grp,
		      int max_samp_grp, double *beta0,
		      double **A0, double a0, double b0,
		      int tau0, double **T0, double *varb, double vars,
		      double *varg, int *counter, int **counterg,
		      int n_gen);
