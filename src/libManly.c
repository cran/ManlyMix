

#include "array.h"
#include "math.h"
#include "Manly.h"
 
#define Inf 1e+140

#include "Manly.h"

#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


void proAgree(int n, int K1, int K2, int *id1, int *id2, double *maxPro, int *combination){

	int sch, i, j, v, w, finish, flag, ind, v1;
	int size, h, idsmall, trclass;
	double curr;

    	int *cn;
	double **pat;

	if (K1 < K2){
		size = K1;
		idsmall = 1;
	} else {
		size = K2;
		idsmall = 2;
	}
		
	curr = 0.0;
	
	sch = 0;
	i = 0;
	j = -1;
	flag = 0;
	finish = 0;
	ind = 0;

	MAKE_MATRIX(pat, size, size);
	for (v=0; v<size; v++){
		for (w=0; w<size; w++){
			pat[v][w] = 0;
		}
	}

	MAKE_VECTOR(cn,size);
	for (v=0; v<size; v++){
		cn[v] = 0;
	}
  
	while (finish == 0){
    
		if (j != (size-1)){
			j = j+1;
		} else {
			if (flag == 1){
				j = 0;
				i = i + 1;
				flag = 0;
			}
		}
    
		if (pat[i][j] == 0){
			for (v=0; v<size; v++){
				pat[i][v] = 1;
				pat[v][j] = 1;
			}
      
			sch = sch + 1;
			cn[sch-1] = j;
			flag = 1;
		}

		if ((sch == size) & (flag == 1)){
      
//			for (v=0; v<size; v++){
//				printf(" %i", cn[v]);
//			}
//			printf(" \n");

//			####################################

			trclass = 0;

			if (idsmall == 1){
				for (h=0; h<n; h++){				
					if (cn[id1[h]] == id2[h]) trclass++;
				}
			} else {
				for (h=0; h<n; h++){				
					if (cn[id2[h]] == id1[h]) trclass++;
				}
			}

			curr = (double)trclass / n;

			if ((*maxPro) < curr) {
	
				(*maxPro) = curr;
				anulli(combination, size);
				for (v1=0; v1<size; v1++){
					combination[v1] = cn[v1];
				}
				
			}
			
//			printf("    %lf\n", curr);

//			####################################

			ind++;
			flag = 0;
			sch = sch - 1;
			i = i - 1;
			j = cn[sch-1];
			sch = sch - 1;
      
			for (v=0; v<size; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 0;
				}
			}

			for (v=0; v<sch; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 1;
					pat[w][cn[v]] = 1;
				}
			}    
      
		}



		if ((j == (size-1)) & (flag == 0)){
			i = i - 1;
			
			sch = sch - 1;

			if (sch >= 0){

				j = cn[sch];

				for (v=0; v<size; v++){
					for (w=0; w<size; w++){
						pat[v][w] = 0;
					}
				}

				if (sch > 0){
					for (v=0; v<sch; v++){
						for (w=0; w<size; w++){
							pat[v][w] = 1;
							pat[w][cn[v]] = 1;
						}
					}    
				}

			}

			if (i >= 0){
				pat[i][j] = 1;
			}

		}

		if (sch == -1){
			finish = 1;
		}

	}

	FREE_MATRIX(pat);
	FREE_VECTOR(cn);

}



// Manly transformation for a column
void Manly_transn(int n, double la, double *x, double *y){
  int i;
  if(fabs(la)<1e-12){
    for(i=0; i<n; i++){
       y[i] = x[i];
    }
  }
  
  else{
    for(i=0; i<n; i++) {
      y[i] = (exp(x[i] * la)-1) / la;
    }
  }

}

// Manly transformation for an observation
void Manly_trans(int p, double *la, double *x, double *y){
  int j;

  for(j=0; j<p; j++) {
    if(fabs(la[j])<1e-12){
      y[j] = x[j];
    } else{
      y[j] = (exp(x[j] * la[j])-1) / la[j];
    }
  }

}






// Manly transformation for n vectors
void Manly_transX(int n, int p, double *la, double **X, double **Y){

  int i, j;

  for (j=0; j<p; j++){
    if (fabs(la[j])<1e-12){
      for(i=0; i<n; i++) {
	Y[i][j] = X[i][j];
      }
    } else{
      for(i=0; i<n; i++) {
	Y[i][j] = (exp(X[i][j] * la[j]) - 1) / la[j];
      }
    }
      
  }
  
}




// Manly density calculation
void Manly_dens(int n, int p, double **X, double *la, double *Mu, double **S, double *dens){

	int i, j;
	double c, detS, Res;
	double *Eig;
	double **Y, **L, **Sinv;

	MAKE_MATRIX(Y, n, p);
	MAKE_MATRIX(Sinv, p, p);
	MAKE_MATRIX(L, p, p);
	MAKE_VECTOR(Eig, p);
 
 	Manly_transX(n, p, la, X, Y);
	
	#ifdef __HAVE_R_
		EigValDec(p, Eig, S, &detS);
	#else
		cephes_symmeigens_down(p, Eig, S, &detS);
	#endif

	Anull(L, p, p);

	for (j=0; j<p; j++){
		L[j][j] = 1 / Eig[j];
	}

	XAXt(S, p, L, Sinv);


	for(i=0; i<n; i++){
		
		vec_(p, Y[i], Mu);
		Res = vAvt(Y[i], p, Sinv);

		dens[i] = exp(-Res / 2.0) / sqrt(pow((2*PI), p) * detS);
		
		c = 0;
		for(j=0; j<p; j++){
			c = c + la[j] * X[i][j];
		
		}

		dens[i] *= exp(c);
	
	}


	FREE_MATRIX(Y);
	FREE_MATRIX(Sinv);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);

}








// Manly mixture
void Manly_mix(int n, int p, int K, double **X, double *tau, double **Mu, double ***S, double **la, double *mix){
  	int i, k;	
	double **gamma, **Ga, *dens;	

	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(Ga, p, p);
	MAKE_VECTOR(dens, n);


	for(k=0; k<K; k++){

		cpy1(S, k, p, p, Ga);
	
		Manly_dens(n, p, X, la[k], Mu[k], Ga, dens);
		
		for(i=0; i<n; i++){

			gamma[i][k] = tau[k] * dens[i];
			
		}	

	}
  
	anull(mix, n);
	for(i=0; i<n; i++){
		for(k=0; k<K; k++){
			mix[i] += gamma[i][k];
		}

	}

	FREE_VECTOR(dens);
	FREE_MATRIX(gamma);
 	FREE_MATRIX(Ga); 
 }







// Manly log likelihood
double Manly_logl(int n, int p, int K, double **X, double *tau, double **Mu, double ***S, double **la){

	int i;
	double *mix, ll;
	
	ll = 0;
	MAKE_VECTOR(mix, n);

	Manly_mix(n, p, K, X, tau, Mu, S, la, mix);
	
	for(i=0; i<n; i++){
		ll = ll + log(mix[i]);
	}
	
	FREE_VECTOR(mix);

	return(ll);

}




// BIC
double Manly_bic(double ll, int n, int M){
	double bic;

	bic = -2 * ll + M * log(n);

	return(bic);  
}




// E-step
void E_step(int n, int K, int p, double **X, double *tau, double **Mu, double ***S, double **la, double **gamma){
	int i, k;
	
	double *dens, *sum_gamma;
	double **Ga;
	MAKE_VECTOR(dens, n);
	MAKE_VECTOR(sum_gamma, n);
	MAKE_MATRIX(Ga, p, p);

	anull(sum_gamma, n);
	
	for(k=0; k<K; k++){
		
		cpy1(S, k, p, p, Ga);

		Manly_dens(n, p, X, la[k], Mu[k], Ga, dens);

		for(i=0; i<n; i++){
			gamma[i][k] = tau[k] * dens[i];
		}		


	}
	
	for(i=0; i<n; i++){

	  for(k=0; k<K; k++){
			sum_gamma[i] += gamma[i][k];
		}

		for(k=0; k<K; k++){
			gamma[i][k] /= sum_gamma[i];
		}

	}		
	
	FREE_VECTOR(sum_gamma);
	FREE_MATRIX(Ga);
	FREE_VECTOR(dens);

}


// Q-function
double Q(int n, int p, double *la_nonzero, int *index, double **X, double *gamma_k){
  
	int i, j, l, count;
	double sum =0;
	double res, *la, *Mu, **Y, **S;
	double detS;
	double *Eig, *jac;


	MAKE_VECTOR(jac, n);
	MAKE_VECTOR(Eig, p);	
	MAKE_VECTOR(Mu, p);
	MAKE_VECTOR(la, p);
	MAKE_MATRIX(Y, n, p);
	MAKE_MATRIX(S, p, p);

	count = 0;

	for(j=0; j<p; j++){
		if(index[j]==1){
			la[j] = la_nonzero[count];
			count += 1;
		}
		else{
			la[j] = 0; 

		}
	}	

	Manly_transX(n, p, la, X, Y);

	anull(Mu, p);
	Anull(S, p, p);
	anull(Eig, p);
	
	for(i=0; i<n; i++){
		sum += gamma_k[i];
	}

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			Mu[j] += Y[i][j] * gamma_k[i];
		}
		Mu[j] /= sum;
	}

	for(j=0; j<p; j++){
		for(l=0; l<p; l++){
			for(i=0; i<n; i++){
				S[j][l] += (Y[i][j] - Mu[j]) * (Y[i][l] - Mu[l]) * gamma_k[i];
			}
			S[j][l] /= sum;
		}
	}	

	#ifdef __HAVE_R_
		EigValDec(p, Eig, S, &detS);
	#else
		cephes_symmeigens_down(p, Eig, S, &detS);
	#endif

	anull(jac, n);

	res = -sum / 2.0 * log(detS);

	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			jac[i] += X[i][j] * la[j]; 
		}
		res = res + gamma_k[i] * jac[i];
	}

	FREE_VECTOR(jac);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Mu);
	FREE_VECTOR(la);
	FREE_MATRIX(Y);
	FREE_MATRIX(S);

	return(-res);

}





// M-step
double M_step(int n, int p, int K, double *misc_double, double **X, double **gamma, double **la, double *tau, double **Mu, double ***S){
  
	int i, j, k, l, count;
	int sum_index, *index;
	double 	eps, *sum_gamma, *gamma_k, *Q_value, Q_value0, min;
	double **Y;

	MAKE_VECTOR(sum_gamma, K);
	MAKE_VECTOR(index, p);
	MAKE_VECTOR(Q_value, K);
	MAKE_VECTOR(gamma_k, n);
	MAKE_MATRIX(Y, n, p);

	anull(sum_gamma, K);
	Anull(Mu, K, p);
	Anull3(S, K, p, p);

 	eps = misc_double[0];
	
	for(k=0; k<K; k++){
		for(i=0; i<n; i++){
		
			sum_gamma[k] += gamma[i][k];

		}

		tau[k] = sum_gamma[k] / n; 
	}

	Q_value0 = 0;
	
	for(k=0; k<K; k++){
		sum_index = 0;
	
		cpyv(gamma, k, n, gamma_k);		
		
		for(j=0; j<p; j++){
			index[j] = (la[k][j] != 0.0);
			sum_index += index[j];
		}

		if(sum_index > 0){

			double *la_nonzero;

			MAKE_VECTOR(la_nonzero, sum_index);
			count = 0;
			for(j=0; j<p; j++){
				if(index[j] == 1){
					la_nonzero[count] = la[k][j];
					count += 1;
				}
			}
			
			min = simplex(Q, n, p, index, X, gamma_k, la_nonzero, eps, 0.1);

			count = 0;
			for(j=0; j<p; j++){
				if(index[j] == 1){
					la[k][j] = la_nonzero[count];
					count += 1;
				}	
				
				else{
					la[k][j] = 0.0;
				}
			}			
			
			FREE_VECTOR(la_nonzero);

			Q_value[k] = min;
		
		} else {
			double *la_nonzero;

			MAKE_VECTOR(la_nonzero, p);

			anull(la_nonzero, p);

			Q_value[k] = Q(n, p, la_nonzero, index, X, gamma_k);

			FREE_VECTOR(la_nonzero);
		}

 		Manly_transX(n, p, la[k], X, Y);

		for(j=0; j<p; j++){
		
			for(i=0; i<n; i++){
		
				Mu[k][j] += Y[i][j] * gamma_k[i];

			}

			Mu[k][j] /= sum_gamma[k];
		
		}


		for(j=0; j<p; j++){
		
			for(l=0; l<p; l++){

				for(i=0; i<n; i++){
		
					S[k][j][l] += (Y[i][j] - Mu[k][j]) * (Y[i][l] - Mu[k][l]) * gamma_k[i];
	
				}

				S[k][j][l] /= sum_gamma[k];
			}

		}

		Q_value0 += Q_value[k];

	}
	
	
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(index);
	FREE_VECTOR(gamma_k);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(Y);
	
	return Q_value0;	

}


// EM algorithm (runs from id vector)
void Manly_EM(int n, int p, int K, double **X, int *id, int max_iter, double *misc_double, double **la, double *tau, double **Mu, double ***S, double **gamma, double *ll, int *conv){

	int i, k, iter;
	double max, Q_value, Q_value_new, eps;

 	eps = misc_double[0];

	iter = 0;
	Q_value = -INFINITY;

	for(i=0; i<n; i++){
		for(k=0; k<K; k++){
			if(id[i] == k+1){
				gamma[i][k] = 1.0;
			}
			else{
				gamma[i][k] = 0.0;
			}
		}
	}


	do{
		Q_value_new = Q_value; 
		
		iter += 1;
		
		Q_value = M_step(n, p, K, misc_double, X, gamma, la, tau, Mu, S);

		E_step(n, K, p, X, tau, Mu, S, la, gamma);
			
	}

	
	while ((iter < max_iter) && (fabs(Q_value_new - Q_value) / fabs(Q_value) > eps));
	
	ll[0] = Manly_logl(n, p, K, X, tau, Mu, S, la);

	conv[0] = iter;
	if(fabs(Q_value_new - Q_value) / fabs(Q_value) <= eps){
		conv[1] = 0;
	} else{
		conv[1] = 1;
	}
	
	anulli(id, n);
	
	for(i=0; i<n; i++){
		max = -INFINITY;
		for(k=0; k<K; k++){
			if(gamma[i][k] > max){
				id[i] = k+1;
				max = gamma[i][k];
			}

		}

	}


}



// EM algorithm (runs from parameter values)
void Manly_EM2(int n, int p, int K, double **X, int max_iter, double *misc_double, double *tau, double **Mu, double ***S, double **la, double **gamma, int *id, double *ll, int *conv){
  
	int i, k, iter;
	double max, Q_value, Q_value_new, eps;

 	eps = misc_double[0];


	iter = 0;
	Q_value = -INFINITY;

	do{
		Q_value_new = Q_value; 
		
		iter += 1;

		E_step(n, K, p, X, tau, Mu, S, la, gamma);

		Q_value = M_step(n, p, K, misc_double, X, gamma, la, tau, Mu, S);


	}

	
	while ((iter < max_iter) && (fabs(Q_value_new - Q_value) / fabs(Q_value) > eps));
	

	ll[0] = Manly_logl(n, p, K, X, tau, Mu, S, la);

	conv[0] = iter;
	if(fabs(Q_value_new - Q_value) / fabs(Q_value) <= eps){
		conv[1] = 0;
	} else {
		conv[1] = 1;
	}

	anulli(id, n);

	for(i=0; i<n; i++){
		max = -INFINITY;
		for(k=0; k<K; k++){
			if(gamma[i][k] > max){
				id[i] = k+1;
				max = gamma[i][k];
			}
		}
	}


}

