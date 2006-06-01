/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai, Jordan R. Vance, and 
  David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
