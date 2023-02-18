/****************************************************************************************
*
*
* My attempt to produce and MCMC algorith for a DP gaussian mixture using 
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
// N - upper bound on the number of mixture components
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
// a - prior shape of the variance from mixture components
// b - prior rate of the variance from mixture components
// alpha - DP precision parameter
// niter - number of MCMC iterates
// nburn - number of MCMC iterates to be discarded
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for component means
// sigma2 - matrix containing MCMC draws for component variances
// pi - matrix containing MCMC draws for component weights
// z - matrix containing MCMC draws for component labels

static void dp_gibbs(double *y, int *n, int *N,
            double *m, double *v, double *a, double *b, double *alpha,
            int *ndens_y, double *y_grid, 
            int *subsample, int *ss_size,
            int *niter, int *nburn, int *nthin,
            double *mu, double *sigma2, double *w,
            int *z, double *density_y){


  int i, j, k, ii, nk;
  int nout = (*niter - *nburn)/(*nthin);
  double mstar, vstar, astar, bstar, tmp, cc, cmc;
  
  double _mu[*N], _sigma2[*N], _w[*N], _v[*N];
  int _z[*n];
  
  /* Initialize variables */
  for(k = 0; k < *N; k++){ 
    _mu[k] = rnorm(0,1);
    _sigma2[k] = rgamma(1,1);
  }

  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*N, 0.25);
    
  }
           
  // initializing everything to do a subsample;
  int arry[*n];
  int n_sub, tmpI, randomIndex; 
  for(j=0; j<*n; j++){
    arry[j] = j;
  }
  n_sub = *ss_size;
  
//  Rprintf("subsample = %d\n", *subsample);
  double ytmp[*n]; 
  if(*subsample == 0){
    for(j=0; j<*n; j++){
      ytmp[j] = y[j];
      arry[j] = j;
    }
    n_sub = *n;    
  }         	 
  
  Rprintf("n_sub = %d\n", n_sub);         
//  RprintIVecAsMat("arry", arry, 1, n_sub);         
//  RprintVecAsMat("ytmp", ytmp, 1, n_sub);         

  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it
  
  // scratch vectors to update parameters
  double sumdens, cprob, uu, sumy=0.0, ssq=0.0;
  double scr1[*N], dens[*N], dens_y[*ndens_y];

  ii = 0;

  for(i=0; i<*niter; i++){
    // if subsample==TRUE, then I need to draw a subsample from 
    // the data.  I try to do this here.
    if(*subsample==1){

      for (j = 0; j < *n; j++) {    // shuffle array
        tmpI = arry[j];
        randomIndex = rand() % (*n);

        arry[j] = arry[randomIndex];
        arry[randomIndex] = tmpI;
        if((i > (*nburn-1)) & ((i+1) % *nthin ==0)){
          z[ii + nout*j] = 0;
        }
      }   
      
      for(j = 0; j < n_sub; j++){
        ytmp[j] = y[arry[j]];
      }
    }
    
//    RprintIVecAsMat("arry", arry, 1, n_sub);
//    RprintVecAsMat("ytmp", ytmp, 1, n_sub);

    // 2) update zi - component labels
    for(j=0; j<n_sub; j++){
      sumdens = 0.0;
      for(k=0; k<*N; k++){
        dens[k] = dnorm(ytmp[j], _mu[k], sqrt(_sigma2[k]), 0)*_w[k]; 
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

//    RprintIVecAsMat("_z", _z, 1, *n);
    
    // 1) update stick-breaking weights
    tmp=1.0;
    for(k=1; k<*N; k++){
      cc = 0.0;
      cmc = 0.0;
      for(j=0; j<n_sub; j++){
        if(_z[j] == (k+1)){
          cc = cc + 1.0;
        }
        if(_z[j] > (k+1)){
          cmc = cmc + 1.0;
        }
      }
      astar = 1.0 + cc;
      bstar = *alpha + cmc;
      
      _v[k-1] = rbeta(astar, bstar);
    
      if(k==(*N-1)) _v[k-1] = 1.0;
      if(k==1) _w[k-1] = _v[k-1];
      if(k>1){
        tmp = tmp*(1-_v[k-2]);
        _w[k-1] = _v[k-1]*tmp;
      }
    }

  
  
    
    // update mu
    for(k=0; k<*N; k++){
      sumy = 0.0;
      nk = 0;
      for(j=0; j<n_sub; j++){
        if(_z[j] == (k+1)){
          sumy = sumy + ytmp[j];  
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
      for(j=0; j<n_sub; j++){
        if(_z[j] == (k+1)){
          ssq = ssq + (ytmp[j]-_mu[k])*(ytmp[j]-_mu[k]);
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
        w[ii + nout*k] = _w[k];
        mu[ii + nout*k] = _mu[k];
        sigma2[ii + nout*k] = _sigma2[k];
      }
      for(j=0; j<n_sub; j++){
//        mu[ii + nout*j] = _mu[_z[j]-1];
//        sigma2[ii + nout*j] = _sigma2[_z[j]-1];
//        Rprintf("arry[j] = %d\n", arry[j]);
//        Rprintf("_z[j]=%d\n", _z[j]);
        z[ii + nout*(arry[j])] = _z[j];
      }
      if(*ndens_y > 0){
        for(j=0; j<*ndens_y; j++){
          dens_y[j] = 0.0;
          for(k=0; k<*N; k++){
            dens_y[j] = dens_y[j] + _w[k]*dnorm(y_grid[j], _mu[k], sqrt(_sigma2[k]), 0);
          }
          density_y[ii + nout*j] = dens_y[j];  
        }  
      }
      
      ii = ii + 1;
    }
  }
}





SEXP DP_GIBBS(SEXP y, SEXP n, SEXP N, SEXP m, SEXP v, SEXP a, SEXP b, SEXP alpha,
              SEXP ndens_y, SEXP y_grid,
	          SEXP subsample, SEXP ss_size, SEXP niter, SEXP nburn, SEXP nthin) {
	          
  int nprot = 0;
  int _n = asInteger(n);
  int _N = asInteger(N);
  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);
  int _ndens_y = asInteger(ndens_y);
  int _subsample = asInteger(subsample);
  int _ss_size = asInteger(ss_size);
  double _m = asReal(m);
  double _v = asReal(v);
  double _a = asReal(a);
  double _b = asReal(b);  
  double _alpha = asReal(alpha);  

  double nout = (_niter-_nburn)/_nthin;
  
  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  y_grid = PROTECT(coerceVector(y_grid, REALSXP)); nprot++;
  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _N)); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++; 
  SEXP DENSITY_Y = PROTECT(allocMatrix(REALSXP, nout, (_ndens_y))); nprot++; 

  double *MUout, *SIGMA2out, *Wout, *DENSITY_Yout;
  int *Zout;
  MUout = REAL(MU);
  SIGMA2out = REAL(SIGMA2);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  DENSITY_Yout = REAL(DENSITY_Y);

  GetRNGstate();

  dp_gibbs(REAL(y), &_n, &_N, &_m, &_v, &_a, &_b, &_alpha, 
           &_ndens_y, REAL(y_grid), 
           &_subsample, &_ss_size, 
           &_niter, &_nburn, &_nthin, 
           MUout, SIGMA2out, Wout, Zout, DENSITY_Yout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 5)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, W);
  SET_VECTOR_ELT(ans, 3, Z);
  SET_VECTOR_ELT(ans, 4, DENSITY_Y);


  SEXP nm = allocVector(STRSXP, 5);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
  SET_STRING_ELT(nm, 2, mkChar("w"));
  SET_STRING_ELT(nm, 3, mkChar("z"));
  SET_STRING_ELT(nm, 4, mkChar("density"));

  UNPROTECT(nprot);
  return(ans);
}



