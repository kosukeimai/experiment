double TruncNorm(double lb, double ub, double mu, double var, int invcdf);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
double dMVN(double *Y, double *MEAN, double **SIG_INV, int dim, int give_log);
void rWish(double **Sample, double **S, int df, int size);
double dnegbin(int Y, double mu, double theta, int give_log);
double rnegbin(double mu, double theta);
