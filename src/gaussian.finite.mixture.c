/****************************************************************************************
*
*
* My attempt to produce and MCMC algorith for a finite gaussian mixture using 
* the .Call function
*
*
****************************************************************************************/
#include "Rutil.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>

#include <math.h>
#include <stdio.h>
#include <stddef.h>

// inputs
// y - data vector
// n - length of y
// N - number of mixture components
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
// a - prior shape of the variance from mixture components
// b - prior rate of the variance from mixture components
// alpha - prior shape parameters for mixture component weights (dirichlet)
// niter - number of MCMC iterates
// nburn - number of MCMC iterates to be discarded
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for component means
// sigma2 - matrix containing MCMC draws for component variances
// pi - matrix containing MCMC draws for component weights
// z - matrix containing MCMC draws for component labels

static void gibbs(double* y, int* n, int *N,
            double* m, double* v, double* a, double* b, double *alpha,
            int* niter, int *nburn, int *nthin,
            double* mu, double* sigma2, double* w,
            int* z){

  int i, j, k, ii, nk;
  int nout = (*niter - *nburn)/(*nthin);
  double mstar, vstar, astar, bstar;
  double alphastar[*N];
  
  double _mu[*N], _sigma2[*N], _w[*N];
  int _z[*n];
  
  /* Initialize variables */
  for(k = 0; k < *N; k++){ 
    _mu[k] = rnorm(0,1);
    _sigma2[k] = rgamma(1,1);
  }

  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*N, 0.25);
  }
                     
  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it
  
  // scratch vectors to update parameters
  double sumdens, cprob, uu, sumy=0.0, ssq=0.0;
  double scr1[*N], dens[*N];

 ii = 0;

  for(i=0; i<*niter; i++){


    // 1) update w
    for(k=0; k<*N; k++){
      alphastar[k] = *alpha;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          alphastar[k] = alphastar[k] + 1;
        }
      }
    }
    ran_dirich(alphastar, *N, scr1, _w);
  
  
    // 2) update zi - component labels
    for(j=0; j<*n; j++){
      sumdens = 0.0;
      for(k=0; k<*N; k++){
        dens[k] = dnorm(y[j], _mu[k], sqrt(_sigma2[k]), 0)*_w[k]; 
        sumdens = sumdens + dens[k];
      }
//      RprintVecAsMat("dens", dens, 1, *N);
      for(k=0; k<*N; k++){
        scr1[k] = dens[k]/sumdens;
      }
//      RprintVecAsMat("scr1", scr1, 1, *N);

	  uu = runif(0.0,1.0);
	  cprob= 0.0;
      for(k = 0; k < *N; k++){
	    cprob = cprob + scr1[k];
		if (uu < cprob){
		  _z[j] = k+1;
		  break;
		}
	  }
    }
    
    // update mu
    for(k=0; k<*N; k++){
      sumy = 0.0;
      nk = 0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          sumy = sumy + y[j];  
          nk = nk+1;
        }
      }

      vstar = 1/((nk/_sigma2[k]) + 1/(*v));
      mstar = vstar*((1.0/_sigma2[k])*sumy + (*m)/(*v));
      
      _mu[k] = rnorm(mstar, sqrt(vstar));
    }


    // update sigma2
    for(k=0; k<*N; k++){
      ssq=0.0, nk=0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          ssq = ssq + (y[j]-_mu[k])*(y[j]-_mu[k]);
          nk = nk + 1;
        }
      }
      astar = 0.5*(nk) + *a;
      bstar = 0.5*ssq + *b;
      _sigma2[k] = 1/rgamma(astar, 1/bstar);
    }
  
    // keep iterates
    if((i > (*nburn-1)) & ((i+1) % *nthin ==0)){
      for(k=0; k<*N; k++){
        mu[ii + nout*k] = _mu[k];
        sigma2[ii + nout*k] = _sigma2[k];
        w[ii + nout*k] = _w[k];
      }
      for(j=0; j<*n; j++){
        z[ii + nout*j] = _z[j];
      }
      ii = ii + 1;
    }
  }
  
}


SEXP GIBBS(SEXP y, SEXP n, SEXP N, SEXP m, SEXP v, SEXP a, SEXP b, SEXP alpha,
	     SEXP niter, SEXP nburn, SEXP nthin) {
  int nprot = 0;
  int _n = asInteger(n);
  int _N = asInteger(N);
  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);
  double _m = asReal(m);
  double _v = asReal(v);
  double _a = asReal(a);
  double _b = asReal(b);  
  double _alpha = asReal(alpha);  

  double nout = (_niter-_nburn)/_nthin;
  
  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _N)); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++; 

  double *MUout, *SIGMA2out, *Wout;
  int *Zout;
  MUout = REAL(MU);
  SIGMA2out = REAL(SIGMA2);
  Wout = REAL(W);
  Zout = INTEGER(Z);

  GetRNGstate();

  gibbs(REAL(y), &_n, &_N, &_m, &_v, &_a, &_b, &_alpha, &_niter, &_nburn, &_nthin, 
        MUout, SIGMA2out, Wout, Zout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 4)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, W);
  SET_VECTOR_ELT(ans, 3, Z);


  SEXP nm = allocVector(STRSXP, 4);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
  SET_STRING_ELT(nm, 2, mkChar("w"));
  SET_STRING_ELT(nm, 3, mkChar("z"));

  UNPROTECT(nprot);
  return(ans);
}



