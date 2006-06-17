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
  /* mean vector for the outcome model */
  double *meano = doubleArray(n_samp);
  /* mean vector for the compliance model */
  double *meanc = doubleArray(n_samp);

  /*** model parameters ***/
  /* probability of Y = 1 for a complier */
  double *pC = doubleArray(n_samp); 
  double *pN = doubleArray(n_samp);
  /* probability of being a complier and never-taker */
  double *qC = doubleArray(n_samp);
  double *qN = doubleArray(n_samp);

  /* probability of being a always-taker */
  double *pA = doubleArray(n_samp);
  double *meana = doubleArray(n_samp);
  /* subset of the data */
  double **Xtemp = doubleMatrix(n_samp+n_covC, n_covC+1);
  int *Atemp = intArray(n_samp);
  
  double **A0C = doubleMatrix(n_covC, n_covC);
  double **A0O = doubleMatrix(n_covO, n_covO);

  /* quantities of interest: ITT, CACE  */
  double ITT, CACE;
  double Y1barC, Y0barC, YbarN, YbarA;
  int n_comp;          /* number of compliers */
  int n_never;
  int p_comp, p_never; /* prob. of being a particular type */

  /*** storage parameters and loop counters **/
  int progress = 1;
  int keep = 1;
  int i, j, k, l, main_loop;  
  int itemp, itempP = ftrunc((double) n_gen/10), itempA, itempC, itempO, itempQ;
  double dtemp;
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
      A0O[j][k] = dA0O[itemp++];

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
  itempA = 0; itempC = 0; itempO = 0; itempQ = 0;   
  for(main_loop = 1; main_loop <= n_gen; main_loop++){

    /** COMPLIANCE MODEL **/    
    bprobitGibbs(C, Xc, betaC, n_samp, n_covC, 0, beta0, A0C, *mda, 1);

    /** ALWAYS-TAKERS MODEL **/
    if (AT) {
      /* subset the data */
      itemp = 0;
      for (i = 0; i < n_samp; i++)
	if (C[i] == 0) {
	  Atemp[itemp] = A[i];
	  for (j = 0; j < n_covC; j++)
	    Xtemp[itemp][j] = Xc[i][j];
	  itemp++;
	}
      for (i = n_samp; i < n_samp + n_covC; i++) {
	for (j = 0; j < n_covC; j++)
	  Xtemp[itemp][j] = Xc[i][j];
	itemp++;
      }
      bprobitGibbs(Atemp, Xtemp, betaA, itemp-n_covC, n_covC, 0, beta0, A0C, *mda, 1);
    }      

    /* Sample complier status for control group */
    if (AT) { /* some always-takers */
      for(i = 0; i < n_samp; i++) {
	meanc[i] = 0;
	for(j = 0; j < n_covC; j++) 
	  meanc[i] += Xc[i][j]*betaC[j];
	qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	meana[i] = 0;
	for(j = 0; j < n_covC; j++) 
	  meana[i] += Xc[i][j]*betaA[j];
	qN[i] = (1-qC[i])*pnorm(meana[i], 0, 1, 0, 0);
	if (Z[i] == 0 & D[i] == 0){
	  if (unif_rand() < (qC[i]*pC[i]/(qC[i]*pC[i]+qN[i]*pN[i]))) { 
	    C[i] = 1; 
	    Xo[i][1] = 1;
	  }
	  else {
	    C[i] = 0; 
	    Xo[i][1] = 0;
	  }  
	}
	if (Z[i] == 1 & D[i] == 1){
	  if (unif_rand() < (qC[i]*pC[i]/(qC[i]*pC[i]+(1-qC[i]-qN[i])*pN[i]))) { 
	    C[i] = 1;
	    A[i] = 0;
	    Xo[i][0] = 1;
	    Xo[i][2] = 0;
	  }
	  else {
	    C[i] = 0;
	    A[i] = 1;
	    Xo[i][0] = 0;
	    Xo[i][2] = 1;
	  }  
	}
      }
    } else { /* no always-takers */
      for(i = 0; i < n_samp; i++)
	if(Z[i] == 0){
	  meanc[i] = 0;
	  for(j = 0; j < n_covC; j++) 
	    meanc[i] += Xc[i][j]*betaC[j];
	  qC[i] = pnorm(meanc[i], 0, 1, 1, 0);
	  if (unif_rand() < (qC[i]*pC[i]/(qC[i]*pC[i]+(1-qC[i])*pN[i]))) { 
	    C[i] = 1; 
	    Xo[i][1] = 1; 
	  }
	  else {
	    C[i] = 0; 
	    Xo[i][1] = 0;
	  }
	}
    }


    /** OUTCOME MODEL **/
    bprobitGibbs(Y, Xo, gamma, n_samp, n_covO, 0, gamma0, A0O, *mda, 1);

    /** Compute probabilities of Y = 1 **/ 
    if (AT) { /* always-takers */
      for (i = 0; i < n_samp; i++) {
	meano[i] = 0;
	for(j = 3; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (Z[i] == 0 & D[i] == 0) {
	  pC[i] = Y[i]*pnorm(meano[i]+gamma[1], 0, 1, 1, 0) +
	    (1-Y[i])*pnorm(meano[i]+gamma[1], 0, 1, 0, 0);
	  pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	    (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	} 
	if (Z[i] == 1 & D[i] == 1) {
	  pC[i] = Y[i]*pnorm(meano[i]+gamma[0], 0, 1, 1, 0) +
	    (1-Y[i])*pnorm(meano[i]+gamma[0], 0, 1, 0, 0);
	  pA[i] = Y[i]*pnorm(meano[i]+gamma[2], 0, 1, 1, 0) +
	    (1-Y[i])*pnorm(meano[i]+gamma[2], 0, 1, 0, 0);
	}
      }
    } else { /* no always-takers */
      for(i = 0; i < n_samp; i++){
	meano[i] = 0;
	for(j = 2; j < n_covO; j++)
	  meano[i] += Xo[i][j]*gamma[j];
	if (Z[i] == 0) {
	  pC[i] = Y[i]*pnorm(meano[i]+gamma[1], 0, 1, 1, 0) + 
	    (1-Y[i])*pnorm(meano[i]+gamma[1], 0, 1, 0, 0);
	  pN[i] = Y[i]*pnorm(meano[i], 0, 1, 1, 0) +
	    (1-Y[i])*pnorm(meano[i], 0, 1, 0, 0);
	}
      } 
    }
    
    /** storing the results **/
    if (main_loop > *iBurnin) {
      if (keep == *iKeep) {
	/** Computing Quantities of Interest **/
	ITT = 0; n_comp = 0; n_never = 0;
	p_comp = 0; p_never = 0; 
	Y1barC = 0; Y0barC = 0; YbarN = 0; YbarA = 0;
	for(i = 0; i < n_samp; i++){
	  p_comp += qC[i];
	  p_never += qN[i];
	  if(C[i] == 1) { /* ITT effects */
	    n_comp++;
	    Y1barC += pnorm(meano[i]+gamma[0], 0, 1, 1, 0);
	    Y0barC += pnorm(meano[i]+gamma[1], 0, 1, 1, 0); 
	    ITT += (pnorm(meano[i]+gamma[0], 0, 1, 1, 0) - 
		     pnorm(meano[i]+gamma[1], 0, 1, 1, 0));
	  } else {
	    if (AT)
	      if (A[i])
		YbarA += pnorm(meano[i]+gamma[2], 0, 1, 1, 0);
	      else {
		n_never++;
		YbarN += pnorm(meano[i], 0, 1, 1, 0);
	      }
	    else {
	      n_never++;
	      YbarN += pnorm(meano[i], 0, 1, 1, 0);
	    }
	  }
	}
	ITT /= n_comp;     /* ITT effect */
	p_comp /= n_samp;  /* ITT effect on D; Prob. of being a complier */
	CACE = ITT/p_comp; /* CACE */
	p_never /= n_samp; /* Prob. of being a never-taker */
	Y1barC /= n_comp; Y0barC /= n_comp; /* E[Y_i(j)|C_i=1] for j=0,1 */
	YbarN /= n_never;
	if (AT)
	  YbarA /= (n_samp-n_comp-n_never);

	QoI[itempQ++] = ITT;   
	QoI[itempQ++] = CACE;   
	QoI[itempQ++] = p_comp; 	  
	QoI[itempQ++] = p_never;
	QoI[itempQ++] = Y1barC;
	QoI[itempQ++] = Y0barC;
	QoI[itempQ++] = YbarN;
	if (AT)
	  QoI[itempQ++] = YbarA;

	if (param) {
	  for (j = 0; j < n_covC; j++)
	    coefC[itempC++] = betaC[j];
	  if (AT)
	    for (j = 0; j < n_covC; j++)
	      coefA[itempA++] = betaA[j];
	  for (j = 0; j < n_covO; j++)
	    coefO[itempO++] = gamma[j];
	}

	keep = 1;
      }
      else
	keep++;
    }

    /** Imputing missing Y **/
    for(i = 0; i < n_samp; i++){
      if (R[i] == 1) { 
	if (AT) 
	  for (j = 0; j < 3; j++)
	    meano[i] += Xo[i][j]*gamma[i];
	else
	  for (j = 0; j < 2; j++)
	    meano[i] += Xo[i][j]*gamma[i];
	if (unif_rand() < pnorm(meano[i], 0, 1, 1, 0))
	  Y[i] = 1;
	else
	  Y[i] = 0;
      }
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
  FreeMatrix(Xc, n_samp+n_covC);
  FreeMatrix(Xo, n_samp+n_covO);
  free(meano);
  free(meanc);
  free(pC);
  free(pN);
  free(pA);
  free(meana);
  FreeMatrix(Xtemp, n_samp+n_covC);
  free(Atemp);
  FreeMatrix(A0C, n_covC);
  FreeMatrix(A0O, n_covO);
  FreeMatrix(mtempC, n_covC);
  FreeMatrix(mtempO, n_covO);

} /* main */

