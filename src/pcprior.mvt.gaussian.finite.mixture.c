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
// y - nxp data matrix
// n - number of observations
// p - dimension of response vector
// K - number of mixture components
// alpha_prior_type - identifies the type of mixture model that is used
//   1 - sparse - symmetric dirichlet
//   2 - centered - asymmetric dirichlet
// alpha_prior_dist - identifies if gamma prior or pc prior with be used for alpha
//   1 - pc prior
//   2 - gamma prior
// U - value that centers asymmetric dirichlet (not needed for sparse model)
// alpha1_val - value for alpha1 if it is fixed
// alpha2_val - value for alpha2 if it is fixed
// update_alpha1 - logical indicating if alpha1 should be updated in MCMC algorithm
// update_alpha2 - logical indicating if alpha2 should be updated in MCMC algorithm
// y_grid - array of y values that are used to estimate density
// ndens_y - number of values in y-grid
// alpha_grid - values of alpha that were used to evaluate the density of alpha
// alpha_density - density values of a corresponding to alpha-grid values
// ndens_alpha - number of values in alpha-grid
// M0 - prior mean vector of the mean from mixture components
// Sigma0 - prior covariance matrix of the mean from mixture components
// A - prior shape of the variance from mixture components
// A0 - prior rate of the variance from mixture components
// a_gam - shape and scale rate parameter when using gamma prior for alpha
// niter - number of MCMC iterates
// nburn - number of MCMC iterates to be discarded
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for component means
// sigma2 - matrix containing MCMC draws for component variances
// w - matrix containing MCMC draws for component weights
// z - matrix containing MCMC draws for component labels
// alpha - array containing MCMC draws for alpha
// mu0 - array containing MCMC draws for mu0
// sigma20 - array containing MCMC draws for sigma20
// density_y - matrix containing MCMC draw for the density estimate of y
//

static void mvt_ifmm_gibbs(double *y,
                           int *n, int *p,
                           int *K,
                           int *alpha_prior_type,
                           int *alpha_prior_dist,
                           int *mu_sigma_prior,
                           int *U,
                           int *basemodel,
                           double *alpha1_val,
                           double *alpha2_val,
                           int *update_alpha1,
                           int *update_alpha2,
                           double *y_grid, int *ndens_y, int *hierarchy,
                           double *alpha_grid, double *alpha_density, int *ndens_alpha,
                           double *alpha1_grid, double *alpha1_density, int *ndens_alpha1,
                           double *alpha2_grid, double *alpha2_density, int *ndens_alpha2,
                           double *M0, double *Sigma0,
                           int *nu0, double *K0,
                           double *A, double *a0, double *b0,
                           double *a_gam, double *b_gam,
                           double *a1_gam, double *b1_gam,
                           double *a2_gam, double *b2_gam,
                           int* niter, int *nburn, int *nthin,
                           double* mu, double* Sigma, double* w, int* z,
                           double *alpha, double *alpha1, double *alpha2,
                           double *mu0, double *sigma20, double *density_y){


  // alpha_prior_type = 1, then sparse.  If alpha_prior_type = 2, then centered
  // alpha_prior_dist = 1, then pc.      If alpha_prior_dist = 2, then gamma.

  Rprintf("alpha_prior_type = %d\n", *alpha_prior_type);
  Rprintf("alpha_prior_dist = %d\n", *alpha_prior_dist);
  Rprintf("update_alpha1 = %d\n", *update_alpha1);
  Rprintf("update_alpha2 = %d\n", *update_alpha2);
  RprintVecAsMat("M0", M0, 1, *p);
  RprintVecAsMat("Sigma0", Sigma0, *p, *p);
  Rprintf("nu0 = %d\n", *nu0);
  RprintVecAsMat("K0", K0, *p, *p);
//  Rprintf("n = %d\n", *n);
//  Rprintf("N = %d\n", *K);

  if(*alpha_prior_type==2) Rprintf("U=%d\n", *U);
  
  int i, j, k, d, ii, dd, nk;
  int nout = (*niter - *nburn)/(*nthin);

  double _mu[(*K)*(*p)], _Sigma[(*K)*(*p)*(*p)], _w[*K];
  int _z[*n];
  
  double Sigma_tmp[(*K)*(*p)*(*p)], mu_tmp[(*K)*(*p)], w_tmp[*K];
  int z_tmp[*n];
  
  double  _alpha, _alpha1, _alpha2;
  double _mu0[*p], _Sigma0[(*p)*(*p)];

  int nk_vec2[*K],  nk_tmp[*n], alpha_loc[*K], alpha_loc_tmp[*K];

  int *nk_vecOrd = (int *) R_alloc(*K, sizeof(int));
  SEXP nk_vec = PROTECT(Rf_allocVector(REALSXP, *K));


  // Initialize variables
  _alpha = 0.5;
  _alpha1 = *alpha1_val;
  _alpha2 = *alpha2_val;

  for(k = 0; k < *K; k++){
    nk_vec2[k] = 0;

    for(d=0; d< *p; d++){
      _mu[k*(*p) + d] = rnorm(0,1);
      mu_tmp[k*(*p) + d] = _mu[k*(*p) + d];
      _mu0[d] = M0[d];

      for(dd=0; dd < *p; dd++){
        _Sigma[k*(*p)*(*p) + (d)*(*p)+dd] = 0.0;
        if(d == dd) _Sigma[k*(*p)*(*p) + (d)*(*p)+dd] = runif(0, 1);

        Sigma_tmp[k*(*p)*(*p) + (d)*(*p)+dd] = _Sigma[k*(*p)*(*p) + (d)*(*p)+dd];
        _Sigma0[d*(*p)+dd]=Sigma0[d*(*p)+dd];
      }
    }


    _w[k] = 1.0/(*K);
    w_tmp[k] = _w[k];

    if(*alpha_prior_type == 2){//asymmetric dirichlet
      if(k<(*U)){
        alpha_loc[k] = 1;
        alpha_loc_tmp[k] = 1;
      }
      if(k>=(*U)){
        alpha_loc[k] = 2;
        alpha_loc_tmp[k] = 2;
      }
    }
  }



  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*K, 0.25);
    z_tmp[j] = _z[j];
//    if(j < 100){_z[j] = 1; z_tmp[j] = 1;}
//    if(j >= 100 & j < 200){_z[j] = 2; z_tmp[j] = 2;}
//    if(j >= 200){_z[j] = 3; z_tmp[j] = 3;}
  }

//  for(k=0; k<*K; k++){
//    REAL(nk_vec)[k] = 0.0;
//    nk_vec2[k]=0;
//    nk_tmp[k]=0;
//  }
//  REAL(nk_vec)[0]=100.0; REAL(nk_vec)[1]=100.0; REAL(nk_vec)[2]=100.0;
//  nk_vec2[0]=100; nk_vec2[1]=100; nk_vec2[2]=100;
//  nk_tmp[0]=100; nk_tmp[1]=100; nk_tmp[2]=100;
  
//  RprintIVecAsMat("alpha_loc_tmp", alpha_loc_tmp, 1, *K);
//  RprintIVecAsMat("alpha_loc", alpha_loc, 1, *K);

  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it

  // scratch vectors to update parameters
  int nustar;
  double sumdens, cprob,  ssq=0.0, summu=0.0, kap2;
  double mstar, vstar, astar, bstar, max_dens, ns, os, ld;
  double p0, an, ao, log_ao_density=0.0, log_an_density, llo,lln,llr,uu;
  double scr1[*K*(*p)*(*p)], scr2[(*p)*(*p)], scr3[(*p)*(*p)], dens[*ndens_y];
  double dens_vals[*K], ao_vec[*K], an_vec[*K], alphastar[*K];;
  double y_tmp[(*p)],  mean_tmp[(*p)], S_tmp[(*p)*(*p)],  Sstar[(*p)*(*p)];
  double iSy[*p], iS0M0[*p], sumy[*p], Mstar[*p];
  double outrmvnorm[*p], outrwish[(*p)*(*p)];
  double *zerov = R_VectorInit(*p, 0.0);

//  double A = 100.0, A0=15.0;

//  RprintVecAsMat("alpha_grid",alpha_grid, 1, *ndens_alpha);
//  RprintVecAsMat("alpha_density",alpha_density, 1, *ndens_alpha);

  // Find the log f(alpha) value
  for(d=0; d<(*ndens_alpha); d++){
    if(log(_alpha) < log(alpha_grid[d])){
      if(d==0){
        log_ao_density = alpha_density[d];
      } else {
        p0 = (_alpha - alpha_grid[d-1])/(alpha_grid[d] - alpha_grid[d-1]);
        log_ao_density = (1-p0)*alpha_density[d-1] + p0*alpha_density[d];
      }
      break;
    }
  }

  Rprintf("log_ao_density = %f\n", log_ao_density);

  // Metropolis tunning parameter
  double csiga = 0.5;

//  Rprintf("Prior values: m0 = %.2f, s20 = %.2f, A = %.2f, A0 = %.2f, \n  %11s a0 = %.2f, b0 = %.2f, nu0 = %.2f, k0 = %.2f\n\n",
//             *m0, *s20, *A, *A0, *a0, *b0, *nu0, *k0);


  if(*mu_sigma_prior == 1){
//    Rprintf("Prior values: m0 = %0.2f, s20 = %0.2f, A = %0.2f\n", *m0, *s20, *A);
  }

  if(*mu_sigma_prior == 2){
//    Rprintf("Prior values: m0 = %0.2f, s20 = %0.2f, a0 = %0.2f, b0 = %0.2f\n", *m0, *s20, *a0, *b0);
  }

  if(*mu_sigma_prior == 3){
//    Rprintf("Prior values: m0 = %0.2f, k0 = %0.2f, nu0 = %0.2f, s20 = %0.2f\n", *m0, *k0, *nu0, *s20);
  }

  ii = 0;

  for(i=0; i<*niter; i++){

    Rprintf("i =========================================== %d\n", i);


    for(k=0; k<*K; k++){
      REAL(nk_vec)[k] = 0.0;
      nk_vec2[k]=0;
      nk_tmp[k]=0;
    }


    // 1) update zi - component labels
    for(j=0; j<*n; j++){

      for(d=0; d<*p; d++){
        y_tmp[d] = y[j*(*p) + d];
        mean_tmp[d] = _mu[0*(*p) + d];
        for(dd = 0; dd < *p; dd++){

          S_tmp[d*(*p)+dd] = _Sigma[0*(*p)*(*p) + (d)*(*p)+dd];
        }
      }

      cholesky(S_tmp, (*p), &ld);
	  inverse_from_cholesky(S_tmp, scr1, scr2, *p);

      max_dens = dmvnorm(y_tmp, mean_tmp, S_tmp, *p, ld, scr1, 1) + log(_w[0]);

      for(k=0; k<*K; k++){
        for(d=0; d<*p; d++){
          mean_tmp[d] = _mu[k*(*p)+d];
          for(dd = 0; dd < *p; dd++){
            S_tmp[d*(*p)+dd] = _Sigma[k*(*p)*(*p) + (d)*(*p)+dd];
          }
        }

        
        cholesky(S_tmp, (*p), &ld);
	    inverse_from_cholesky(S_tmp, scr2, scr3, *p);
	    
        dens_vals[k] = dmvnorm(y_tmp, mean_tmp, S_tmp, *p, ld, scr1, 1) + log(_w[k]);

        if(dens_vals[k] > max_dens) max_dens = dens_vals[k];
      }

      sumdens = 0.0;
      for(k=0; k<*K; k++){
        dens_vals[k] = exp(dens_vals[k] - max_dens);
        sumdens = sumdens + dens_vals[k];
      }
      for(k=0; k<*K; k++){
        scr1[k] = dens_vals[k]/sumdens;
      }

	  uu = runif(0.0,1.0);
	  cprob= 0.0;
      for(k = 0; k < *K; k++){
	    cprob = cprob + scr1[k];
		if (uu < cprob){
		  _z[j] = k+1;
		  z_tmp[j] = k+1;
		  REAL(nk_vec)[k] = REAL(nk_vec)[k] + 1.0;
		  nk_vec2[k] = nk_vec2[k] + 1;
		  nk_tmp[k] = nk_tmp[k] + 1;
		  break;
		}
	  }
    }
    
    
//    RprintIVecAsMat("nk_vec2", nk_vec2, 1, *K);
//    RprintIVecAsMat("nk_tmp", nk_tmp, 1, *K);
//    RprintVecAsMat("nk_vec", REAL(nk_vec), 1, *K);

    // Should I relabel here and make the first component always that which has the most
    // number of members?
    RprintIVecAsMat("z = ", _z, 1, *n);
    RprintIVecAsMat("z_tmp", z_tmp, 1, *n);

    if(*alpha_prior_type==2){ // centered prior need to reorder based on cluster size
      R_orderVector(nk_vecOrd, *K, Rf_lang1(nk_vec), TRUE, TRUE);
      for(k=0; k<*K; k++){

        for(j=0; j<*n; j++){
          if(z_tmp[j] == (nk_vecOrd[k]+1)) _z[j] = k+1;
        }
        _w[k] = w_tmp[nk_vecOrd[k]];
        alpha_loc[k] = alpha_loc_tmp[nk_vecOrd[k]];
        nk_vec2[k] = nk_tmp[nk_vecOrd[k]];
        for(d=0; d<(*p); d++){
          _mu[k*(*p) + d] = mu_tmp[nk_vecOrd[k]*(*p) + d];
          for(dd=0; dd<(*p); dd++){
            _Sigma[k*(*p)*(*p) + (d)*(*p)+dd] = Sigma_tmp[nk_vecOrd[k]*(*p)*(*p) + (d)*(*p)+dd];
          }
        }
      }
//      rsort_with_index()
    }
//    RprintIVecAsMat("nk_vec2", nk_vec2, 1, *K);


    

    RprintIVecAsMat("z_tmp", z_tmp, 1, *n);
    RprintIVecAsMat("_z", _z, 1, *n);

    // 2) update w
    for(k=0; k<*K; k++){
      if(*alpha_prior_type == 1){//symmetric dirichlet
        alphastar[k] = _alpha;
      }
      if(*alpha_prior_type == 2){//asymmetric dirichlet
        if(alpha_loc[k] == 1){
          alphastar[k] = _alpha1;
        }
        if(alpha_loc[k] == 2){
          alphastar[k] = _alpha2;
        }
      }

      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          alphastar[k] = alphastar[k] + 1;
        }
      }
    }

    ran_dirich(alphastar, *K, scr1, _w);

    for(k=0; k<*K; k++){
      if(_w[k] < 1e-323) _w[k] = 1e-323;
      w_tmp[k] = _w[k];
      alpha_loc_tmp[k] = alpha_loc[k];
    }


    RprintVecAsMat("_w", _w, 1, *K);
//    Rprintf("alpha = %f\n", _alpha);



    // 3) update Mu
    for(k=0; k<*K; k++){
    
      for(d=0; d<(*p); d++) sumy[d] = 0.0;

      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          for(d=0; d<(*p); d++) sumy[d] = sumy[d] + y[j*(*p) + d];
        }
      }

      for(d=0; d<*p; d++){
        mean_tmp[d] = _mu[k*(*p)+d];
        for(dd = 0; dd < *p; dd++){
          S_tmp[d*(*p)+dd] = _Sigma[k*(*p)*(*p) + (d)*(*p)+dd];
        }
      }

      cholesky(S_tmp, (*p), &ld);
	  inverse_from_cholesky(S_tmp, scr1, scr2, *p);

      matrix_product(S_tmp, sumy, iSy, *p, 1, *p);

      cholesky(_Sigma0, (*p), &ld);
	  inverse_from_cholesky(_Sigma0, scr1, scr2, *p);

      matrix_product(_Sigma0, M0, iS0M0, *p, 1, *p);


	  for(d = 0; d < (*p); d++){

	    scr1[d] = iSy[d] + iS0M0[d];

	    for(dd = 0; dd < (*p); dd++){
	      // Note that S_tmp and _Sigma0 are inverses right now.
		  Sstar[d*(*p)+dd] = nk_vec2[k]*S_tmp[d*(*p)+dd] + _Sigma0[d*(*p)+dd];
		}
	  }

//      RprintVecAsMat("iSy + iS0M0", scr1, 1, *p);

	  cholesky(Sstar, (*p), &ld);
	  inverse_from_cholesky(Sstar, scr2, scr3, (*p));
//      RprintVecAsMat("Sstar", Sstar, *p, *p);


//      RprintVecAsMat("iSy", iSy, 1, *p);

	  matrix_product(Sstar, scr1, Mstar, (*p), 1, (*p));

//      RprintVecAsMat("Mstar", Mstar, 1, *p);

	  cholesky(Sstar, (*p) , &ld);

	  ran_mvnorm(Mstar, Sstar, (*p), scr2, outrmvnorm);

      // This is the independent prior for muk,sigma2k so
      // muk ~ N(mu0, sigma20) and sigma2k ~ UN or sigma2k ~ IG
/*      if(*mu_sigma_prior != 3 ){
        vstar = 1.0/((nk/_sigma2[k]) + 1.0/(_sigma20));
        mstar = vstar*((1.0/_sigma2[k])*sumy + (_mu0)/(_sigma20));
      }

      // This is for NIG prior on (muk|sigma^2k)(sigma^2k)
      if(*mu_sigma_prior == 3){
        vstar = _sigma2[k]/(nk + *k0);
        mstar = (sumy + (*k0)*(_mu0))/(nk + *k0);
      }
*/
      for(d=0; d<(*p); d++){
        _mu[k*(*p) + d] = outrmvnorm[d];
        mu_tmp[k*(*p) + d] = outrmvnorm[d];;
      }
    }

    RprintVecAsMat("muk = ", _mu, *K, *p);





    //
    // 4) update Sigma
    //
    for(k=0; k<*K; k++){
      
      Rprintf("k = %d\n", k);

      // This is y|z ~ MVN(mu_z, Sigma_z); mu_z ~ MVN; Sigma_z ~ IW
      if(*mu_sigma_prior == 1){
    
        for(d=0; d < (*p)*(*p); d++) Sstar[d] = K0[d];
      
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            for(d=0; d<*p; d++){
              for(dd=0; dd<*p; dd++){
                Sstar[d*(*p)+dd] = Sstar[d*(*p)+dd]  +
                                            (y[j*(*p)+d] -  _mu[k*(*p) + d])*
                                            (y[j*(*p)+dd] - _mu[k*(*p) + dd]);

              }
            }
          }
        }

        RprintVecAsMat("K0", K0, *p, *p);
        RprintVecAsMat("Sstar", Sstar, *p, *p);

        nustar = nk_vec2[k] + *nu0;
//	    Rprintf("nk_vec2[k] = %d\n", nk_vec2[k]);
//	    Rprintf("nu0 = %d\n", *nu0);
//      Rprintf("nustar = %d\n", nustar);

        cholesky(Sstar, *p, &ld);

	    // I need to invert Astar to get an inverse wishart from the ran_wish function;
	    inverse_from_cholesky(Sstar, scr1, scr2, *p);
	    cholesky(Sstar, *p, &ld); // Sstar is now the choleksy decomposition of Sstar^{-1};

	    ran_wish(nustar, Sstar, *p, scr1, scr2, zerov, outrwish);
	
	    RprintVecAsMat("outrwish", outrwish, (*p), (*p));

	    cholesky(outrwish, *p, &ld);
	    inverse_from_cholesky(outrwish, scr1, scr2, *p);

	    RprintVecAsMat("outrwish", outrwish, (*p), (*p));


        for(d = 0; d < (*p)*(*p); d++){
	      _Sigma[k*(*p)*(*p) + d] = outrwish[d];
	      Sigma_tmp[k*(*p)*(*p) + d] = outrwish[d];
	    }
      }
     
      // This is y|z ~ MVN(mu_z, kappa^2_zI); mu_z ~ MVN; kappa_z ~ UN(0,A)
      if(*mu_sigma_prior == 2){
        
        kap2 = _Sigma[k*(*p)*(*p) + (0)*(*p)+0];
        os = sqrt(kap2);
        ns = rnorm(os, 0.5);
        if(ns > 0){
          llo=0.0, lln=0.0;
          ssq=0.0; nk=0;
          for(j=0; j<*n; j++){
            if(_z[j] == (k+1)){
              for(d=0; d<*p; d++){
                ssq = ssq + (y[j*(*p)+d] -  _mu[k*(*p) + d])*(y[j*(*p)+d] -  _mu[k*(*p) + d]);
              }
              nk = nk + 1; 
            }
          }
          llo = -0.5*(*p)*nk*log(os) - (1/os)*0.5*ssq + dunif(os, 0, *A, 1);
          lln = -0.5*(*p)*nk*log(ns) - (1/ns)*0.5*ssq + dunif(ns, 0, *A, 1);;

          llr = lln - llo;
          uu = runif(0,1);
          if(llr > log(uu)){
            for(d=0; d<*p; d++){
              for(dd=0; dd<*p; dd++){
                if(d==dd){
                  _Sigma[k*(*p)*(*p) + (d)*(*p)+dd] = ns*ns;
                  Sigma_tmp[k*(*p)*(*p) + (d)*(*p)+dd] = ns*ns;
                } 
              }
            }
          } 
        }
      }

      // This is y|z ~ MVN(mu_z, kappa^2_zI); mu_z ~ MVN; kappa_z ~ IG(a0,b0)
      if(*mu_sigma_prior == 3){
        ssq=0.0, nk=0;
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            for(d=0; d<*p; d++){
              ssq = ssq + (y[j*(*p)+d] -  _mu[k*(*p) + d])*(y[j*(*p)+d] -  _mu[k*(*p) + d]);
            }
            nk = nk + 1;
          }
        }
        astar = 0.5*(nk)*(*p) + *a0;
        bstar = 0.5*ssq + *b0;

        kap2 = 1/rgamma(astar, 1/bstar);
        for(d=0; d<*p; d++){
          for(dd=0; dd<*p; dd++){
            if(d==dd){
              _Sigma[k*(*p)*(*p) + (d)*(*p)+dd] = kap2;
              Sigma_tmp[k*(*p)*(*p) + (d)*(*p)+dd] = kap2;
            } 
          }
        }
      }

/*
      // This is the dependent prior for (muk,sigma2k) ~ NIW(m0,k0,v0,s20)
      if(*mu_sigma_prior == 4){
        ssq=0.0, nk=0;
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            ssq = ssq + (y[j]-_mu[k])*(y[j]-_mu[k]);
            nk = nk + 1;
          }
        }
        astar = 0.5*(nk + 1 + *nu0);
        bstar = 0.5*ssq + 0.5*(*k0)*(mu[k]- _mu0)*(mu[k]- _mu0) + 2*((*s20)/(*nu0));

        _Sigma[k] = 1/rgamma(astar, 1/bstar);
      }

*/

    }
//   RprintVecAsMat("_Sigma", _Sigma, *K, (*p)*(*p));

//    RprintVecAsMat("sigma2 = ", _sigma2, 1, *K);




//

/*

    if(*hierarchy==1){

      // 5) update mu0
      summu = 0.0;
      for(k=0; k<*K; k++){
        summu = summu + _mu[k];
      }
      vstar = 1.0/((*K)/_sigma20 + 1.0/(*s20));
      mstar = vstar*((1.0/_sigma20)*summu + (*m0)/(*s20));

      _mu0 = rnorm(mstar, sqrt(vstar));


      // 6) update sigma20
      os = sqrt(_sigma20);
      ns = rnorm(os, 0.5);
      if(ns > 0){
        llo = 0.0, lln = 0.0;
        for(k=0; k<*K; k++){
          llo = llo + dnorm(_mu[k], _mu0, os, 1);
          lln = lln + dnorm(_mu[k], _mu0, ns, 1);
        }
        llo = llo + dunif(os, 0, *A0, 1);
        lln = lln + dunif(ns, 0, *A0, 1);

        llr = lln - llo;
        uu = runif(0,1);

        if(llr > log(uu)) _sigma20 = ns*ns;

      }

      // ssq = 0.0;
      // for(k=0; k<*K; k++){
      //  ssq = ssq + (_mu[k] - *m0)*(_mu[k] - *m0);
      // }
      // astar = 0.5*(*K) + *a;
      // bstar = 0.5*ssq + *b;
      // _sigma20 = 1.0/rgamma(astar, 1.0/bstar);
    }

*/

    // 7) update alpha

    // This code is for the sparse model for K+ with symmetric dirichlet
    if(*alpha_prior_type==1){ // sparse model

      // pc prior rather than gamma
      if(*alpha_prior_dist==1){ // pc prior
        if(*basemodel==0){ // model is mixture with 1 component

          ao = _alpha;
          an = rnorm(ao, csiga);

          if((an > alpha_grid[0]) & (an < alpha_grid[*ndens_alpha-1])){
            // Need to find the density value of the induce prior for ao and an;
            for(d=0; d<(*ndens_alpha); d++){
              if(an < alpha_grid[d]){
                if(d==0){
                  log_an_density = alpha_density[d];
                } else {
                  p0 = (an - alpha_grid[d-1])/(alpha_grid[d] - alpha_grid[d-1]);
                  log_an_density = (1-p0)*alpha_density[d-1] + p0*alpha_density[d];
                }
                break;
              }
            }
            for(k=0; k<*K; k++){
              ao_vec[k] = ao;
              an_vec[k] = an;
            }


            llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density);
            lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density);
//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);

            llr = lln - llo;
//          Rprintf("llr = %f\n", llr);

            uu = runif(0,1);

//          Rprintf("log(uu) = %f\n", log(uu));
            if(llr > log(uu)){
              _alpha = an;
              log_ao_density = log_an_density;
            }
//          Rprintf("_alpha = %f\n", _alpha);
//          Rprintf("log_ao_density = %f\n", log_ao_density);
//          Rprintf("log_an_density = %f\n", log_an_density);
          }
        }


        // I have to use the density of the log(a) to update this parameter.
        // I am unsure if the jacobian comes into play or not here.

        if(*basemodel==1){// base model is mixture with K components

//        ao = log(_alpha);
          ao = _alpha;
          an = rnorm(ao, csiga);

//        if((exp(an) > alpha_density[0]) & (exp(an) < alpha_grid[*ndens_alpha-1])){
          if((an > alpha_grid[0]) & (an < alpha_grid[*ndens_alpha-1])){

            // Need to find the density value of the induce prior for ao and an;
            for(d=0; d<(*ndens_alpha); d++){
//              if(an < log(alpha_grid[d])){
              if(an < alpha_grid[d]){
                if(d==0){
                  log_an_density = alpha_density[d];
                } else {
//                  p0 = (exp(an) - alpha_grid[d-1])/(alpha_grid[d] - alpha_grid[d-1]);
                  p0 = (an - alpha_grid[d-1])/(alpha_grid[d] - alpha_grid[d-1]);
                  log_an_density = (1-p0)*alpha_density[d-1] + p0*alpha_density[d];
                }
                break;
              }
            }

            for(k=0; k<*K; k++){
              ao_vec[k] = ao;
              an_vec[k] = an;
            }


            llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density);
            lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density);
  //        Rprintf("llo = %f\n", llo);
  //        Rprintf("lln = %f\n", lln);

            llr = lln - llo;

            uu = runif(0,1);

            if(llr > log(uu)){
              _alpha = an;
              log_ao_density = log_an_density;
            }
          }
        }

      }
      // Use gamma prior rather than pc prior, still symmetric dirichlet
      if(*alpha_prior_dist==2){ // gamma prior
        ao = _alpha;
        an = rnorm(ao, csiga);
        if(an > 0){
          for(k=0; k<*K; k++){
            ao_vec[k] = ao;
            an_vec[k] = an;
          }

          llo = ddirich(_w, ao_vec, *K, 1) + dgamma(ao, *a_gam, *b_gam, 1);
          lln = ddirich(_w, an_vec, *K, 1) + dgamma(an, *a_gam, *b_gam, 1);

          llr = lln - llo;
          uu = runif(0,1);

          if(llr > log(uu)){
            _alpha = an;
          }
        }
      }
    }



    // This code is for the centered prior on K+, now using an asymmetric dirichlet
    if(*alpha_prior_type==2){ // centered prior on K+

      // update a1
      if(*update_alpha1 == 1){
        if(*alpha_prior_dist==1){// pc prior
          ao = _alpha1;
          an = rnorm(ao, csiga);

          if((an > alpha1_grid[0]) & (an < alpha1_grid[*ndens_alpha1-1])){

            // Need to find the density value of the induce prior for ao and an;
            for(d=0; d<(*ndens_alpha1); d++){
              if(an < alpha1_grid[d]){
                if(d==0){
                  log_an_density = alpha1_density[d];
                } else {
                  p0 = (an - alpha1_grid[d-1])/(alpha1_grid[d] - alpha1_grid[d-1]);
                  log_an_density = (1-p0)*alpha1_density[d-1] + p0*alpha1_density[d];
                }
                break;
              }
            }

            for(k=0; k<*K; k++){
              if(alpha_loc[k]==1){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(alpha_loc[k]==2){
                ao_vec[k] = _alpha2;
                an_vec[k] = _alpha2;
              }
            }

            llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density);
            lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density);
//            Rprintf("llo = %f\n", llo);
//            Rprintf("lln = %f\n", lln);



            llr = lln - llo;

            uu = runif(0,1);

            if(llr > log(uu)){
              _alpha1 = an;
              log_ao_density = log_an_density;
            }
          }

        }

        if(*alpha_prior_dist == 2){ // gamma prior
          ao = _alpha1;
          an = rnorm(ao, csiga);
          if(an > 0){
            for(k=0; k<*K; k++){
              if(alpha_loc[k]==1){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(alpha_loc[k]==2){
                ao_vec[k] = _alpha2;
                an_vec[k] = _alpha2;
              }
            }
            llo = ddirich(_w, ao_vec, *K, 1) + dgamma(ao, *a1_gam, *b1_gam, 1);
            lln = ddirich(_w, an_vec, *K, 1) + dgamma(an, *a1_gam, *b1_gam, 1);

            llr = lln - llo;
            uu = runif(0,1);

            if(llr > log(uu)){
              _alpha1 = an;
            }
          }
        }
      }
      // update a2
      if(*update_alpha2 == 1){

        if(*alpha_prior_dist==1){// pc prior

          ao = _alpha2;
          an = rnorm(ao, csiga);

          if((an > alpha2_grid[0]) & (an < alpha2_grid[*ndens_alpha2-1])){

            // Need to find the density value of the induce prior for ao and an;
            for(d=0; d<(*ndens_alpha2); d++){
              if(an < alpha2_grid[d]){
                if(d==0){
                  log_an_density = alpha2_density[d];
                } else {
                  p0 = (an - alpha2_grid[d-1])/(alpha2_grid[d] - alpha2_grid[d-1]);
                  log_an_density = (1-p0)*alpha2_density[d-1] + p0*alpha2_density[d];
                }
                break;
              }
            }

            for(k=0; k<*K; k++){
              if(alpha_loc[k]==1){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(alpha_loc[k]==2){
                ao_vec[k] = _alpha2;
                an_vec[k] = _alpha2;
              }
            }

            llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density);
            lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density);
  //        Rprintf("llo = %f\n", llo);
  //        Rprintf("lln = %f\n", lln);



            llr = lln - llo;

            uu = runif(0,1);

            if(llr > log(uu)){
              _alpha2 = an;
              log_ao_density = log_an_density;
            }
          }


        }

        if(*alpha_prior_dist==2){// gamma prior
          ao = _alpha2;
          an = rnorm(ao, (1e-3*csiga));
          if(an > 0){
            for(k=0; k<*K; k++){
              if(alpha_loc[k]==1){
                ao_vec[k] = _alpha1;
                an_vec[k] = _alpha1;
              }
              if(alpha_loc[k]==2){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
            }
//            RprintVecAsMat("an_vec", an_vec, 1, *K);
//            llo = ddirich(_w, ao_vec, *K, 1) + dunif(ao, 1e-10, 1e-3, 1);
//            lln = ddirich(_w, an_vec, *K, 1) + dunif(an, 1e-10, 1e-3, 1);
            llo = ddirich(_w, ao_vec, *K, 1) + dgamma(ao, *a2_gam, *b2_gam, 1);
            lln = ddirich(_w, an_vec, *K, 1) + dgamma(an, *a2_gam, *b2_gam, 1);
//            llo = ddirich(_w, ao_vec, *K, 1) + dgamma(ao, 1.0, 1.0/((*a_gam)*(*K-*U)), 1); //scale
//            lln = ddirich(_w, an_vec, *K, 1) + dgamma(an, 1.0, 1.0/((*a_gam)*(*K-*U)), 1);

            llr = lln - llo;
            uu = runif(0,1);

            if(llr > log(uu)){
              _alpha2 = an;
            }
          }
        }
      }
    }
//    Rprintf("alpha1 = %f\n", alpha1);


    // keep iterates
    if((i > (*nburn-1)) & ((i) % *nthin ==0)){
//      Rprintf("ii = %d\n", ii);
      for(k=0; k<*K; k++){
        for(d=0; d<*p; d++){
          mu[ii + nout*(k*(*p)+d)] = _mu[k*(*p)+d];
          for(dd=0; dd<*p; dd++){
            Sigma[ii + nout*(k*(*p)*(*p) + (d)*(*p)+dd)] = _Sigma[k*(*p)*(*p) + (d)*(*p)+dd];
          }
        }
        w[ii + nout*k] = _w[k];
      }
      for(j=0; j<*n; j++){
        z[ii + nout*j] = _z[j];
      }
//      sigma20[ii] = _sigma20;
      alpha[ii] = _alpha;
      alpha1[ii] = _alpha1;
      alpha2[ii] = _alpha2;


//      if(*ndens_y > 0){
//        for(j=0; j<*ndens_y; j++){
//          dens[j] = 0.0;
//          for(k=0; k<*K; k++){
//            dens[j] = dens[j] + _w[k]*dnorm(y_grid[j], _mu[k], sqrt(_sigma2[k]), 0);
//          }
//          density_y[ii + nout*j] = dens[j];
//        }
//      }


      ii = ii + 1;
    }


  }

  UNPROTECT(1); // seems like you forgot this
}




SEXP MVT_IFMM_GIBBS(SEXP y, SEXP n, SEXP p, SEXP K,
                SEXP alpha_prior_type, SEXP alpha_prior_dist,
                SEXP mu_sigma_prior, SEXP U,
                SEXP basemodel,
                SEXP alpha1_val, SEXP alpha2_val,
                SEXP update_alpha1, SEXP update_alpha2,
                SEXP y_grid, SEXP ndens_y, SEXP hierarchy,
                SEXP alpha_grid, SEXP alpha_density, SEXP ndens_alpha,
                SEXP alpha1_grid, SEXP alpha1_density, SEXP ndens_alpha1,
                SEXP alpha2_grid, SEXP alpha2_density, SEXP ndens_alpha2,
                SEXP M0, SEXP Sigma0,
                SEXP nu0, SEXP K0,
                SEXP A, SEXP a0, SEXP b0,
                SEXP a_gam, SEXP b_gam, SEXP a1_gam, SEXP b1_gam,
                SEXP a2_gam, SEXP b2_gam,
 	            SEXP niter, SEXP nburn, SEXP nthin){


  int nprot = 0;

  int _n = asInteger(n);
  int _p = asInteger(p);
  int _K = asInteger(K);
  int _alpha_prior_type = asInteger(alpha_prior_type);
  int _alpha_prior_dist = asInteger(alpha_prior_dist);
  int _mu_sigma_prior = asInteger(mu_sigma_prior);
  int _U = asInteger(U);
  int _basemodel = asInteger(basemodel);
  int _update_alpha1 = asInteger(update_alpha1);
  int _update_alpha2 = asInteger(update_alpha2);
  int _ndens_y = asInteger(ndens_y);
  int _hierarchy = asInteger(hierarchy);
  int _ndens_alpha = asInteger(ndens_alpha);
  int _ndens_alpha1 = asInteger(ndens_alpha1);
  int _ndens_alpha2 = asInteger(ndens_alpha2);
  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);
  double _alpha1_val = asReal(alpha1_val);
  double _alpha2_val = asReal(alpha2_val);
  int _nu0 = asInteger(nu0);
  double _A = asReal(A);
  double _a0 = asReal(a0);
  double _b0 = asReal(b0);  
  double _a_gam = asReal(a_gam);
  double _b_gam = asReal(b_gam);
  double _a1_gam = asReal(a1_gam);
  double _b1_gam = asReal(b1_gam);
  double _a2_gam = asReal(a2_gam);
  double _b2_gam = asReal(b2_gam);

  double nout = (_niter-_nburn)/_nthin;


  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  y_grid = PROTECT(coerceVector(y_grid, REALSXP)); nprot++;
  alpha_grid = PROTECT(coerceVector(alpha_grid, REALSXP)); nprot++;
  alpha1_grid = PROTECT(coerceVector(alpha1_grid, REALSXP)); nprot++;
  alpha2_grid = PROTECT(coerceVector(alpha2_grid, REALSXP)); nprot++;
  alpha_density = PROTECT(coerceVector(alpha_density, REALSXP)); nprot++;
  alpha1_density = PROTECT(coerceVector(alpha1_density, REALSXP)); nprot++;
  alpha2_density = PROTECT(coerceVector(alpha2_density, REALSXP)); nprot++;
  M0 = PROTECT(coerceVector(M0, REALSXP)); nprot++;
  Sigma0 = PROTECT(coerceVector(Sigma0, REALSXP)); nprot++;
  K0 = PROTECT(coerceVector(K0, REALSXP)); nprot++;
  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, (_K)*(_p))); nprot++;
  SEXP SIGMA = PROTECT(allocMatrix(REALSXP, nout, (_K)*(_p)*(_p))); nprot++;
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++;
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++;
  SEXP ALPHA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA1 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP MU0 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP SIGMA20 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP DENSITY_Y = PROTECT(allocMatrix(REALSXP, nout, (_ndens_y))); nprot++;


  double *MUout, *SIGMAout, *Wout, *MU0out, *SIGMA20out;
  double  *ALPHAout,  *ALPHA1out,  *ALPHA2out, *DENSITY_Yout;
  int *Zout;
  MUout = REAL(MU);
  SIGMAout = REAL(SIGMA);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  MU0out = REAL(MU0);
  SIGMA20out = REAL(SIGMA20);
  ALPHAout = REAL(ALPHA);
  ALPHA1out = REAL(ALPHA1);
  ALPHA2out = REAL(ALPHA2);
  DENSITY_Yout = REAL(DENSITY_Y);

  GetRNGstate();

  mvt_ifmm_gibbs(REAL(y), &_n, &_p, &_K, &_alpha_prior_type, &_alpha_prior_dist,
              &_mu_sigma_prior, &_U,
              &_basemodel,
              &_alpha1_val, &_alpha2_val,
              &_update_alpha1, &_update_alpha2,
              REAL(y_grid), &_ndens_y, &_hierarchy,
              REAL(alpha_grid), REAL(alpha_density), &_ndens_alpha,
              REAL(alpha1_grid), REAL(alpha1_density), &_ndens_alpha1,
              REAL(alpha2_grid), REAL(alpha2_density), &_ndens_alpha2,
              REAL(M0), REAL(Sigma0), &_nu0, REAL(K0), &_A, &_a0, &_b0,
              &_a_gam, &_b_gam, &_a1_gam, &_b1_gam, &_a2_gam, &_b2_gam,
              &_niter, &_nburn, &_nthin,
              MUout, SIGMAout, Wout, Zout, ALPHAout, ALPHA1out, ALPHA2out,
              MU0out, SIGMA20out, DENSITY_Yout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 10)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA);
  SET_VECTOR_ELT(ans, 2, W);
  SET_VECTOR_ELT(ans, 3, Z);
  SET_VECTOR_ELT(ans, 4, ALPHA);
  SET_VECTOR_ELT(ans, 5, ALPHA1);
  SET_VECTOR_ELT(ans, 6, ALPHA2);
  SET_VECTOR_ELT(ans, 7, MU0);
  SET_VECTOR_ELT(ans, 8, SIGMA20);
  SET_VECTOR_ELT(ans, 9, DENSITY_Y);


  SEXP nm = allocVector(STRSXP, 10);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("Sigma"));
  SET_STRING_ELT(nm, 2, mkChar("w"));
  SET_STRING_ELT(nm, 3, mkChar("z"));
  SET_STRING_ELT(nm, 4, mkChar("alpha"));
  SET_STRING_ELT(nm, 5, mkChar("alpha1"));
  SET_STRING_ELT(nm, 6, mkChar("alpha2"));
  SET_STRING_ELT(nm, 7, mkChar("mu0"));
  SET_STRING_ELT(nm, 8, mkChar("sigma20"));
  SET_STRING_ELT(nm, 9, mkChar("density"));

  UNPROTECT(nprot);
  return(ans);

}



