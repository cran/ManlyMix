#ifndef MANLY_H
#define MANLY_H
void Manly_transn(int n, double la, double *x, double *y);
void Manly_trans(int p, double *la, double *x, double *y);
void Manly_transX(int n, int p, double *la, double **X, double **Y);
int vec_(int a, double *Res, double *Y);
double vAvt(double *v, int p, double **A);
void Manly_dens(int n, int p, double **X, double *la, double *Mu, double **S, double *dens);
void Manly_mix(int n, int p, int K, double **X, double *tau, double **Mu, double ***S, double **la, double *mix);
double Manly_logl(int n, int p, int K, double **X, double *tau, double **Mu, double ***S, double **la);
double Manly_bic(double ll, int n, int M);
void E_step(int n, int K, int p, double **X, double *tau, double **Mu, double ***S, double **la, double **gamma);
double Q(int n, int p, double *la_nonzero, int *index, double **X, double *gamma_k);
void cpyv(double **A, int col, int nrows, double *V);
void Anull3(double ***X, int ax, int bx, int cx);
double M_step(int n, int p, int K, double *misc_double, double **X, double **gamma, double **la, double *tau, double **Mu, double ***S);
void Manly_EM(int n, int p, int K, double **X, int *id, int max_iter, double *misc_double, double **la, double *tau, double **Mu, double ***S, double **gamma, double *ll, int *conv);
void Manly_EM2(int n, int p, int K, double **X, int max_iter, double *misc_double, double *tau, double **Mu, double ***S, double **la, double **gamma, int *id, double *ll, int *conv);
void array1to2(int a, int b, double *y, double **x);
void array2to1(int a, int b, double *y, double **x);
void array1to3(int a, int b, int c, double *y, double ***x);
void array3to1(int a, int b, int c, double *y, double ***x);
double simplex(double (*func)(int, int, double *, int *, double **, double *), int n1, int p, int *index, double **X, double *gamma_k, double *start, double EPSILON, double scale);

void EigValDec(int size, double *W, double **A, double (*determinant));
void Anull(double **X, int ax, int bx);
void anull(double *x, int p);
void anulli(int *x, int p);
void XAXt(double **X, int p, double **A, double **Res);
void cpy1(double ***a, int k, int nrows, int ncols, double **b);
void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);


void E_stepk(int n, int K, int p, double **X, double **Mu, double *sigma2, double **la, double **distance);
double Qk(int n, int p, double *la, double **X);
double M_stepk(int n, int p, int K, double *misc_double, double **X, int *id, double **la, double **Mu, double *sigma2);
void Manly_CEM(int n, int p, int K, double **X, int *id, int max_iter, double *misc_double, double **la, double **Mu, double *sigma2, int *conv);
void Manly_CEM2(int n, int p, int K, double **X, int max_iter, double *misc_double, double **la, double **Mu, double *sigma2, int *id, int *conv);
double simplexk(double (*func)(int, int, double *, double **), int n1, int p, double *start, double **X, double EPSILON, double scale);
void extract(int n, int p, double **X, int *index, double **Y);


/* WCC: "libEVD.c" and "libEVD_LAPACK.c" */
#ifndef __HAVE_R_
	void cephes_symmeigens_down(int p, double *eval, double **A, double (*determinant));
#else
	void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
	void EigValDec(int size, double *W, double **A, double (*determinant));
#endif

#endif /* MANLY_H */
