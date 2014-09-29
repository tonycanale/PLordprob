//  Pairwise log-likelihood for multivariate ordered probit model
//  C code by Euloge Clovis Kenne Pagui and Antonio Canale

#include <R.h>   
#include <Rmath.h>  #include <R_ext/RS.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rdefines.h> 
#include <R_ext/Linpack.h>

#define MAT2(mat, i, j, nI) (mat[i+j*nI])

// To compute the bivariate normal cdf we use the code by A. Genz with minor modification by A. Azzalini in his package "mnormt"

void F77_NAME(sadmvn)( int* , double* , double*, int* , double* , int* , double* , double* , double* , double * , int*) ;

void mult_pmnorm( int *nvar , double *lower , double *upper , int *infin , double *corr , int *maxpts , double *abseps , double *releps , double *esterror , double *result , int *fail )
{ 
  F77_CALL(sadmvn)( nvar, lower, upper, infin, corr, maxpts, abseps, releps, esterror, result, fail ) ;
} 




void cat_pair_llik_real2( double *res, double *data, double *corr, double *mu, double *xi, double *sig2, 
		      double *tresh, int *n, int *q, int *mod, int *two)
{
	
  // The function computes the pairwise log-likelihood based on the multivariate ordered probit model
 	
  // res: calculation's result
  // data: observed data
  // corr: array of correlation coefficients
  // mu: mean structure
  // xi: item-specific means
  // sig2: variance if the data
  // tresh: vector of thresholds
  // n: sample size
  //	q: number of variables
  //  mod: number of levels of the categorical responses	
  //  two: dimension of the integral	
  register int i , r,s;
  int infin[2], maxpts=2000**q, fail=100;
  double corr_s, res2=0.0, lower[2], upper[2], abseps=1.0e-6, releps=0.0, esterror=10.0;
  //int *TWO, ini_TWO[1];
  //TWO = &ini_TWO[0];
  //TWO = 2;

  for( i = 0; i < *n; i++ ){
    for(r = 0; r < (*q-1); r++){
      for(s = (r+1); s < *q; s++){
	if(data[i+r**n]==1){
	  lower[0] = 0.0;
	  upper[0] = (tresh[0]-mu[i]-xi[r])/sqrt(*sig2);
	  infin[0] = 0;
	}
	if(data[i+r**n]==*mod){
	  lower[0] = (tresh[(*mod-2)]-mu[i]-xi[r])/sqrt(*sig2);
	  upper[0] = 0.0;
	  infin[0] = 1;
	}
	if(data[i+r**n]>1 & data[i+r**n] < *mod){
	  lower[0] = (tresh[(int)data[i+r**n]-2]-mu[i]-xi[r])/sqrt(*sig2);
	  upper[0] = (tresh[(int)data[i+r**n]-1]-mu[i]-xi[r])/sqrt(*sig2);
	  infin[0] = 2;
	}
	if(data[i+s**n]==1){
	  lower[1] = 0.0;
	  upper[1] = (tresh[0]-mu[i]-xi[s])/sqrt(*sig2);
	  infin[1] = 0;
	}
	if(data[i+s**n]==*mod){
	  lower[1] = (tresh[(*mod-2)]-mu[i]-xi[s])/sqrt(*sig2);
	  upper[1] = 0.0;
	  infin[1] = 1;
	}
	if(data[i+s**n]>1 & data[i+s**n] < *mod){
	  lower[1] = (tresh[(int)data[i+s**n]-2]-mu[i]-xi[s])/sqrt(*sig2);
	  upper[1] = (tresh[(int)data[i+s**n]-1]-mu[i]-xi[s])/sqrt(*sig2);
	  infin[1] = 2;
	}
	corr_s = MAT2(corr, r, s, *q);
	mult_pmnorm( two, lower, upper, infin, &corr_s, &maxpts, &abseps, &releps, &esterror, &res2, &fail );
	if(res2<=0){
	  res2 = 0.0;
	}
	*res += log(res2);
	res2=0.0;
      }
    }
  }
}
