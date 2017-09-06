
#include "array.h"
#include "Manly.h"





void runProAgree(int (*n), int (*K1), int (*K2), int *id1, int *id2, double (*maxPro), int * combination){

	double maxPro1;

	int n1, K11, K21;

	n1 = (*n);
	K11 = (*K1);
	K21 = (*K2);		

	maxPro1 = (*maxPro);
	
	proAgree(n1, K11, K21, id1, id2, &maxPro1, combination);

	(*maxPro) = maxPro1;

}




void run_Manly_transX(double *x1, int *misc_int, double *la1, double *y1){


	double **X, **Y;

	int p, n;


	p = misc_int[0];
	n = misc_int[1];


	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(Y, n, p);
	
	array1to2(n, p, x1, X);
	
        Manly_transX(n, p, la1, X, Y);	

	array2to1(n, p, y1, Y);
		
	FREE_MATRIX(X);
	FREE_MATRIX(Y);


}




// Manly mixture
void run_Manly_mix(double *x1, int *misc_int, double *tau, double *Mu1, double *S1, double *la1, double *mix){


	double **X, **la, **Mu, ***S;

	int p, n, K;


	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];	


	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	
	array1to2(n, p, x1, X);
	array1to2(K, p, la1, la);
	array1to2(K, p, Mu1, Mu);
	array1to3(K, p, p, S1, S);

	
	Manly_mix(n, p, K, X, tau, Mu, S, la, mix);
	
	array2to1(K, p, la1, la);
	array2to1(K, p, Mu1, Mu);
	array3to1(K, p, p, S1, S);
		
	FREE_MATRIX(X);
	FREE_MATRIX(la);
	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);

}


void run_Manly(double *x1, int *id, int *misc_int, double *misc_double, double *la1, double *tau, double *Mu1, double *S1, double *gamma1, double *ll, int *conv){
		
	double **X, **la, **Mu, ***S, **gamma;

	int p, n, K, max_iter;


	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];	
	max_iter = misc_int[3];	

	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	MAKE_MATRIX(gamma, n, K);
	
	array1to2(n, p, x1, X);
	array1to2(K, p, la1, la);
	array1to2(K, p, Mu1, Mu);
	array1to3(K, p, p, S1, S);
	array1to2(n, K, gamma1, gamma);

	
	Manly_EM(n, p, K, X, id, max_iter, misc_double, la, tau, Mu, S, gamma, ll, conv);
	
	array2to1(K, p, la1, la);
	array2to1(K, p, Mu1, Mu);
	array3to1(K, p, p, S1, S);
	array2to1(n, K, gamma1, gamma);
		
	FREE_MATRIX(X);
	FREE_MATRIX(la);
	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);
	FREE_MATRIX(gamma);
		
}






void run_Manly2(double *x1, int *misc_int, double *misc_double, double *la1, double *tau, double *Mu1, double *S1, double *gamma1, int *id, double *ll, int *conv){
		
	double **X, **la, **Mu, ***S, **gamma;

	int p, n, K, max_iter;


	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];	
	max_iter = misc_int[3];	

	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	MAKE_MATRIX(gamma, n, K);
	
	array1to2(n, p, x1, X);
	array1to2(K, p, la1, la);
	array1to2(K, p, Mu1, Mu);
	array1to3(K, p, p, S1, S);
	array1to2(n, K, gamma1, gamma);

	
	Manly_EM2(n, p, K, X, max_iter, misc_double, tau, Mu, S, la, gamma, id, ll, conv);
	
	array2to1(K, p, la1, la);
	array2to1(K, p, Mu1, Mu);
	array3to1(K, p, p, S1, S);
	array2to1(n, K, gamma1, gamma);
		
	FREE_MATRIX(X);
	FREE_MATRIX(la);
	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);
	FREE_MATRIX(gamma);
		
}




void run_Manlyk(double *x1, int *id, int *misc_int, double *misc_double, double *la1, double *Mu1, double *S1, int *conv){
		
	double **X, **la, **Mu;

	int p, n, K, max_iter;


	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];	
	max_iter = misc_int[3];	

	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(Mu, K, p);
	
	array1to2(n, p, x1, X);
	array1to2(K, p, la1, la);
	array1to2(K, p, Mu1, Mu);
	
	Manly_CEM(n, p, K, X, id, max_iter, misc_double, la, Mu, S1, conv);
	
	array2to1(K, p, la1, la);
	array2to1(K, p, Mu1, Mu);

		
	FREE_MATRIX(X);
	FREE_MATRIX(la);
	FREE_MATRIX(Mu);
		
}






void run_Manlyk2(double *x1, int *misc_int, double *misc_double, double *la1, double *Mu1, double *S1, int *id, int *conv){
		
	double **X, **la, **Mu;

	int p, n, K, max_iter;


	p = misc_int[0];
	n = misc_int[1];
	K = misc_int[2];	
	max_iter = misc_int[3];	

	MAKE_MATRIX(X, n, p);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(Mu, K, p);

	
	array1to2(n, p, x1, X);
	array1to2(K, p, la1, la);
	array1to2(K, p, Mu1, Mu);
	
	Manly_CEM2(n, p, K, X, max_iter, misc_double, la, Mu, S1, id, conv);
	
	array2to1(K, p, la1, la);
	array2to1(K, p, Mu1, Mu);


	FREE_MATRIX(X);
	FREE_MATRIX(la);
	FREE_MATRIX(Mu);

		
}


