/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai, Jordan R. Vance, and 
  David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stdlib.h>
#include <assert.h>

int *intArray(int num);
int **intMatrix(int row, int col);

double *doubleArray(int num);
double **doubleMatrix(int row, int col);
double ***doubleMatrix3D(int x, int y, int z);

long *longArray(int num);

void FreeMatrix(double **Matrix, int row);
void FreeintMatrix(int **Matrix, int row);
void Free3DMatrix(double ***Matrix, int index, int row);
