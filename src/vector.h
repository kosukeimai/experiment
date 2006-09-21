#include <stdlib.h>
#include <assert.h>

int *intArray(int num);
void PintArray(int *ivector, int length);
int **intMatrix(int row, int col);
void PintMatrix(int **imatrix, int row, int col);

double *doubleArray(int num);
void PdoubleArray(double *dvector, int length);
double **doubleMatrix(int row, int col);
void PdoubleMatrix(double **dmatrix, int row, int col);

double ***doubleMatrix3D(int x, int y, int z);
void PdoubleMatrix3D(double ***dmatrix3D, int x, int y, int z);

long *longArray(int num);

void FreeMatrix(double **Matrix, int row);
void FreeintMatrix(int **Matrix, int row);
void Free3DMatrix(double ***Matrix, int index, int row);
