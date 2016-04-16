 

#include "array.h"
#include "math.h"
 
#define Inf 1e+140

#include "Manly.h"

#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif





//Manly E-step

void E_stepk(int n, int K, int p, double **X, double **Mu, double *sigma2, double **la, double **distance){
	int i, j, k;
	double c, Res;
	double **Y;

	MAKE_MATRIX(Y, n, p);
	
	
	for(k=0; k<K; k++){

		Manly_transX(n, p, la[k], X, Y);
		 

		for(i=0; i<n; i++){

	 		c = 0;
			for(j=0; j<p; j++){
				c = c + la[k][j] * X[i][j];
		
			}

			
			Res = 0;
			vec_(p, Y[i], Mu[k]);
			for(j=0; j<p; j++){
				Res = Res + pow(Y[i][j], 2);
		
			}
			
			distance[i][k] = p / 2.0 * log(sigma2[k]) + 1.0 / 2.0 / sigma2[k] * Res - c;
			


		}		

	}

	FREE_MATRIX(Y);

}





// Q-function
double Qk(int n, int p, double *la, double **X){
  
	int i, j;
	double sigma2;
	double res, *Mu, **Y;
	double *jac;


	MAKE_VECTOR(jac, n);	
	MAKE_VECTOR(Mu, p);
	MAKE_MATRIX(Y, n, p);


	Manly_transX(n, p, la, X, Y);


	anull(Mu, p);

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			Mu[j] += Y[i][j];
		}
		Mu[j] /= n;
	}

	sigma2 = 0;
	for(i=0; i<n; i++){

		vec_(p, Y[i], Mu);

		for(j=0; j<p; j++){
			sigma2 = sigma2 + pow(Y[i][j], 2);
		
		}

	}
	sigma2 = 1.0 * sigma2 / p / n;

	anull(jac, n);

	
	res = -p * n / 2.0 * log(sigma2); 


	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			jac[i] += X[i][j] * la[j]; 
		}
		res = res + jac[i];
	}

	FREE_VECTOR(jac);
	FREE_VECTOR(Mu);
	FREE_MATRIX(Y);
	return(-res);

}


// M_step

double M_stepk(int n, int p, int K, double *misc_double, double **X, int *id, double **la, double **Mu, double *sigma2){
  
	int i, j, k;
	int *index, *sum_id;
	double 	eps, *Q_value, Q_value0=0, min;
	double **Y, **Z;

	MAKE_VECTOR(sum_id, K);
	MAKE_VECTOR(index, n);
	MAKE_VECTOR(Q_value, K);


	Anull(Mu, K, p);
	anull(sigma2, K);
 	eps = misc_double[0];
	
	anulli(sum_id, K);
	
	for(k=0; k<K; k++){

		anulli(index, n);
		for(i=0; i<n; i++){

			if(id[i] == k+1) {
				sum_id[k] += 1;
				index[i] = 1;

			
			}		
			
		}

		
		MAKE_MATRIX(Z, sum_id[k], p);
		MAKE_MATRIX(Y, sum_id[k], p);

		extract(n, p, X, index, Z);
		
		min = simplexk(Qk, sum_id[k], p, la[k], Z, eps, 0.1);


		Q_value[k] = min;


		Manly_transX(sum_id[k], p, la[k], Z, Y);

		FREE_MATRIX(Z);
		
		for(j=0; j<p; j++){
			for(i=0; i<sum_id[k]; i++){
				Mu[k][j] += Y[i][j];
			}
			Mu[k][j] /= sum_id[k];
		}

		
		for(i=0; i<sum_id[k]; i++){

			vec_(p, Y[i], Mu[k]);

			for(j=0; j<p; j++){
				sigma2[k] = sigma2[k] + pow(Y[i][j], 2);
		
			}

		}
		sigma2[k] = 1.0 * sigma2[k] / p / sum_id[k];

		
		FREE_MATRIX(Y);
		Q_value0 += Q_value[k];


	}


	//printf("%lf\n", Q_value0);
	FREE_VECTOR(sum_id);
	FREE_VECTOR(index);
	FREE_VECTOR(Q_value);
	
	return Q_value0;	

}



void Manly_CEM(int n, int p, int K, double **X, int *id, int max_iter, double *misc_double, double **la, double **Mu, double *sigma2, int *conv){

	int i, k, iter;
	double min, Q_value, Q_value_new, eps, **distance;


	MAKE_MATRIX(distance, n, K);

 	eps = misc_double[0];

	iter = 0;
	Q_value = -INFINITY;


	do{
		Q_value_new = Q_value;
 
		iter += 1;
		
		Q_value = M_stepk(n, p, K, misc_double, X, id, la, Mu, sigma2);


		E_stepk(n, K, p, X, Mu, sigma2, la, distance);

		for(i=0; i<n; i++){
			min = INFINITY;
			for(k=0; k<K; k++){
				if(distance[i][k] < min){
					id[i] = k+1;
					min = distance[i][k];
				}

			}

		}

			
	}


	
	while ((iter < max_iter) && (fabs(Q_value_new - Q_value) / fabs(Q_value) > eps));



	conv[0] = iter;
	if(fabs(Q_value_new - Q_value) / fabs(Q_value) <= eps){
		conv[1] = 0;
	} else{
		conv[1] = 1;
	}
	
	FREE_MATRIX(distance);

}




void Manly_CEM2(int n, int p, int K, double **X, int max_iter, double *misc_double, double **la, double **Mu, double *sigma2, int *id, int *conv){

	int i, k, iter;
	double min, Q_value, Q_value_new, eps, **distance;


	MAKE_MATRIX(distance, n, K);

 	eps = misc_double[0];

	iter = 0;
	Q_value = -INFINITY;


	do{
		Q_value_new = Q_value;
 

		iter += 1;
		

		E_stepk(n, K, p, X, Mu, sigma2, la, distance);
		//printf("%lf %lf\n", distance[0][0], distance[0][1]);
		//printf("%lf %lf\n", distance[1][0], distance[1][1]);

		for(i=0; i<n; i++){
			min = INFINITY;
			for(k=0; k<K; k++){
				if(distance[i][k] < min){
					id[i] = k+1;
					min = distance[i][k];
				}

			}

		}


		Q_value = M_stepk(n, p, K, misc_double, X, id, la, Mu, sigma2);

		//printf("%lf %lf\n", Mu[0][0], Mu[0][1]);
		//printf("%lf %lf\n", Mu[1][0], Mu[1][1]);


			
	}


	
	while ((iter < max_iter) && (fabs(Q_value_new - Q_value) / fabs(Q_value) > eps));

	//printf("%d\n",iter);
	//printf("%lf %lf\n",Mu[0][0], Mu[0][1]);

	//printf("%lf %lf\n",la[0][0], la[0][1]);
	

	conv[0] = iter;
	if(fabs(Q_value_new - Q_value) / fabs(Q_value) <= eps){
		conv[1] = 0;
	} else{
		conv[1] = 1;
	}
	

	FREE_MATRIX(distance);

}



