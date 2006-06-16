/******************************************************************
  This file is a part of eco: R Package for Fitting Bayesian Models
  of Ecological Inference for 2x2 Tables
  by Kosuke Imai and Ying Lu
  Copyright: GPL version 2 or later.
*******************************************************************/

void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dinv2D(double *X, int size, double *X_inv,char* emsg);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
double ddet2D(double **X, int size, int give_log);
void dcholdc2D(double *X, int size, double *L);
void matrixMul(double **A, double **B, int r1, int c1, int r2, int c2, double **C); 
