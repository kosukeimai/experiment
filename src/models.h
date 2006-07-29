void logitMetro(int *Y, double **X, double *beta, int n_samp,      
		int n_dim, int n_cov, double *beta0, double **A0,     
		double *Var, int n_gen, int *counter);

void bNormalReg(double **D, double *beta, double *sig2, 
		int n_samp, int n_cov, int addprior, int pbeta, 
		double *beta0, double **A0, int psig2, double s0, 
		int nu0, int sig2fixed);
  
void bprobitGibbs(int *Y, double **X, double *beta, int n_samp, 
		  int n_cov, int prior, double *beta0, double **A0, 
		  int mda, int n_gen);

void bprobitMixedGibbs(int *Y, double **X,  double **Z, 
		       double ***Zgrp, int *grp, double *beta, 
		       double **gamma, double **Psi, int n_samp, 
		       int n_fixed, int n_random, int n_grp, 
		       int *n_samp_grp, int prior, double *beta0, 
		       double **A0, int tau0, double **T0, int mda, 
		       int n_gen);

void bNormalMixedGibbs(double *Y, double **X, double **Z, double ***Zgrp,
		       int *grp, double *beta, double **gamma, double *sig2,    
		       double **Psi, int n_samp, int n_fixed, int n_random,
		       int n_grp, int *n_samp_grp, int max_samp_grp,
		       int prior, double *beta0, double **A0, int imp,
		       int nu0, double s0, int tau0, double **T0, int n_gen0); 
