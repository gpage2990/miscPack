/****************************************************************************************
*
*
* My attempt to produce and MCMC algorith for a finite gaussian mixture using
* the .Call function
*
*
****************************************************************************************/
#include "Rutil.h"
#include "matrix.h"

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
// Xmat - n x p covariate/design matrix
// n - length of y
// N - number of mixture components
// p - number of columns in X matrix (this may or may not include an intercept)
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

static void gibbs_reg(double* y,
                      double *Xmat,
                      int* n, int *N, int *p,
                      double* m, double* v, double* a, double* b, double *alpha,
                      int* niter, int *nburn, int *nthin,
                      double* beta, double* sigma2, double* w,
                      int* z, double* fitted_line){

  int i, j, k, ii, d, dd, csobs;
  int nout = (*niter - *nburn)/(*nthin);
  double astar, bstar, xb, ld;
  double alphastar[*N];



  double _beta[(*N)*(*p)], _sigma2[*N], _w[*N];
  int _z[*n];

  /* Initialize variables */
  for(k = 0; k < *N; k++){
    for(d = 0; d < *p; d++){
      _beta[d*(*p) + k] = rnorm(0,1);
    }
    _sigma2[k] = rgamma(1,1);
  }

  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*N, 0.25);
  }

  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it

  // scratch vectors to update parameters
  double sumdens, cprob, uu, ssq=0.0;
  double scr1[*N], dens[*N], scr2[*n], scr3[*n], y_tmp[*n], X_tmp[(*n)*(*p)];
  double Xb[*n], outrmvnorm[*p], sumXy[*p], Xy[*p], tX[(*n)*(*p)], XtX[(*p)*(*p)];
  double Sstar[(*p)*(*p)], Mstar[*p];

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
        xb = 0.0;
        for(d=0; d<*p; d++){
          xb = xb + Xmat[j*(*p)+d]*_beta[k*(*p) + d];
        }
        dens[k] = dnorm(y[j], xb, sqrt(_sigma2[k]), 0)*_w[k];
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


    //RprintIVecAsMat("z", _z, 1, *n);



    // update beta and sigma2 in same loop
    for(k=0; k<*N; k++){
    
    // 3a) update beta
	  csobs=0;
	  for(j=0; j<*n; j++){
	    if(_z[j] == k+1){
	      y_tmp[csobs] = y[j];
	      for(d=0; d<*p; d++){
	        X_tmp[csobs*(*p) + d] = Xmat[j*(*p) + d];
	      }
	      csobs = csobs + 1;
        }
	  }
      
      mat_transpose(X_tmp, tX, csobs, *p);
      
      matrix_product(tX, y_tmp, Xy, *p, 1, csobs);

      matrix_product(tX, X_tmp, XtX, *p, *p, csobs);
      
      
	  for(d=0; d<(*p); d++){
	    sumXy[d] = (1.0/_sigma2[k])*Xy[d];

	    for(dd = 0; dd < (*p); dd++){	      
		  Sstar[d*(*p)+dd] = (1.0/_sigma2[k])*XtX[d*(*p)+dd];
		  if(d == dd){
		    Sstar[d*(*p)+dd] = (1.0/_sigma2[k])*XtX[d*(*p)+dd] + (1.0/(*v));
		  }
		}
	  	  
	  }

	  cholesky(Sstar, (*p), &ld);
	  inverse_from_cholesky(Sstar, scr2, scr3, (*p));

	  matrix_product(Sstar, sumXy, Mstar, (*p), 1, (*p));
	  

	  cholesky(Sstar, (*p) , &ld);


	  ran_mvnorm(Mstar, Sstar, (*p), scr2, outrmvnorm);

	  for(d=0; d < (*p); d++){
	    _beta[k*(*p) + d] = outrmvnorm[d];
	  }

      matrix_product(X_tmp, outrmvnorm, Xb, csobs, 1, (*p));
      
//      RprintVecAsMat("_beta", _beta, *N, *p);
//	    RprintVecAsMat("Xb", Xb, 1, csobs);


      // 3b) update sigma2
      ssq=0.0;
      for(j=0; j<csobs; j++){
      	ssq = ssq + (y_tmp[j] - Xb[j])*(y_tmp[j] - Xb[j]);
      }

      astar = 0.5*(csobs) + *a;
      bstar = 0.5*ssq + *b;

      //bstar is rate and rgamma requires scale hence inverse

     _sigma2[k] = 1/rgamma(astar, 1/bstar);
//      Rprintf("sigma2 = %f\n", _sigma2[j]);

    }


    // keep iterates
    if((i > (*nburn-1)) & ((i+1) % *nthin ==0)){
      dd = 0;
      for(k=0; k<*N; k++){
        sigma2[ii + nout*k] = _sigma2[k];
        w[ii + nout*k] = _w[k];
        for(d=0; d<*p; d++){
          beta[ii + nout*dd] = _beta[k*(*p)+d];
          dd = dd + 1;
        }
      }
      dd = 0;
      for(j=0; j<*n; j++){
        z[ii + nout*j] = _z[j];

        xb = 0.0;
	    for(d = 0; d < *p; d++){
	      xb = xb + _beta[(_z[j]-1)*(*p)+d]*Xmat[j*(*p)+d];	        
	    }
	    fitted_line[ii + nout*dd] = xb;
	    dd = dd+1;
      }

      ii = ii + 1;
    }
  
  }

}


SEXP GIBBS_REG(SEXP y, SEXP Xmat, SEXP n, SEXP N, SEXP p,
               SEXP m, SEXP v, SEXP a, SEXP b, SEXP alpha,
	           SEXP niter, SEXP nburn, SEXP nthin) {
  int nprot = 0;
  int _n = asInteger(n);
  int _N = asInteger(N);
  int _p = asInteger(p);
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
  Xmat = PROTECT(coerceVector(Xmat, REALSXP)); nprot++;
  SEXP BETA = PROTECT(allocMatrix(REALSXP, nout, (_N)*(_p))); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++;
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++;
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++;
  SEXP FITTED_LINE = PROTECT(allocMatrix(REALSXP, nout, (_n))); nprot++;

  double *BETAout, *SIGMA2out, *Wout, *FITTED_LINEout;
  int *Zout;
  BETAout = REAL(BETA);
  SIGMA2out = REAL(SIGMA2);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  FITTED_LINEout = REAL(FITTED_LINE);

  GetRNGstate();

  gibbs_reg(REAL(y), REAL(Xmat), &_n, &_N, &_p, &_m, &_v, &_a, &_b, &_alpha, &_niter, &_nburn, &_nthin,
            BETAout, SIGMA2out, Wout, Zout, FITTED_LINEout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 5)); nprot++;
  SET_VECTOR_ELT(ans, 0, BETA);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, W);
  SET_VECTOR_ELT(ans, 3, Z);
  SET_VECTOR_ELT(ans, 4, FITTED_LINE);


  SEXP nm = allocVector(STRSXP, 5);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("beta"));
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
  SET_STRING_ELT(nm, 2, mkChar("w"));
  SET_STRING_ELT(nm, 3, mkChar("z"));
  SET_STRING_ELT(nm, 4, mkChar("fitted_line"));

  UNPROTECT(nprot);
  return(ans);
}



