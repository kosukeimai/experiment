/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai, Jordan R. Vance, and 
  David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

double TruncNorm(double lb, double ub, double mu, double var, int invcdf);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);


