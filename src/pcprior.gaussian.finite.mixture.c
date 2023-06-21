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
// K - number of mixture components
// alpha_prior_type - identifies the type of mixture model that is used
//   1 - sparse - symmetric dirichlet
//   2 - centered - asymmetric dirichlet
// alpha_prior_dist - identifies if gamma prior or pc prior with be used for alpha
//   1 - pc prior
//   2 - gamma prior
// U - value that centers asymmetric dirichlet (not needed for sparse model)
// y_grid - array of y values that are used to estimate density
// ndens_y - number of values in y-grid
// alpha_grid - values of alpha that were used to evaluate the density of alpha
// alpha_density - density values of a corresponding to alpha-grid values
// ndens_alpha - number of values in alpha-grid
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
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

static void ifmm_gibbs(double* y, int* n, int *K,
            int *alpha_prior_type, int *alpha_prior_dist, 
            int *mu_sigma_prior, int *U,
            int *basemodel, 
            double *alpha1_val, double *alpha2_val,
            int *update_alpha1, int *update_alpha2,
            double *y_grid, int *ndens_y, int *hierarchy,
            double *alpha_grid, double *alpha_density, int *ndens_alpha,
            double *alpha1_grid, double *alpha1_density, int *ndens_alpha1,
            double *alpha2_grid, double *alpha2_density, int *ndens_alpha2,
            double *m0, double *s20, double *A, double *A0,
            double *a0, double *b0, double *nu0, double *k0,
            double *a_gam, double *b_gam,
            double *a1_gam, double *b1_gam,
            double *a2_gam, double *b2_gam,
            int* niter, int *nburn, int *nthin,
            double* mu, double* sigma2, double* w, int* z,
            double *alpha, double *alpha1, double *alpha2,
            double *mu0, double *sigma20, double *density_y){


  // alpha_prior_type = 1, then sparse.  If alpha_prior_type = 2, then centered
  // alpha_prior_dist = 1, then pc.      If alpha_prior_dist = 2, then gamma.

  Rprintf("alpha_prior_type = %d\n", *alpha_prior_type);
  Rprintf("alpha_prior_dist = %d\n", *alpha_prior_dist);
  Rprintf("update_alpha1 = %d\n", *update_alpha1);
  Rprintf("update_alpha2 = %d\n", *update_alpha2);
//  Rprintf("n = %d\n", *n);
//  Rprintf("N = %d\n", *K);
  if(*alpha_prior_type==2) Rprintf("U=%d\n", *U);
  int i, j, k, d, ii, nk;
  int nout = (*niter - *nburn)/(*nthin);
  double alphastar[*K];

  double _mu[*K], _sigma2[*K], _w[*K], _alpha, _alpha1, _alpha2, _mu0, _sigma20;
  double mu_tmp[*K], sigma2_tmp[*K], w_tmp[*K];
  int _z[*n], nk_vec2[*K], z_tmp[*n], alpha_loc[*K], alpha_loc_tmp[*K];
  
  int *nk_vecOrd = (int *) R_alloc(*K, sizeof(int));
  SEXP nk_vec = PROTECT(Rf_allocVector(REALSXP, *K));


  // Initialize variables
  _alpha = 0.5;
  _alpha1 = *alpha1_val;
  _alpha2 = *alpha2_val;
  _mu0 = *m0;
  _sigma20 = *s20;
  for(k = 0; k < *K; k++){
    nk_vec2[k] = 0;

    _mu[k] = rnorm(0,1);
    mu_tmp[k] = _mu[k];

    _sigma2[k] = runif(0, *A);
    sigma2_tmp[k] = _sigma2[k];

    _w[k] = 1/(*K);
    w_tmp[k] = _w[k];
  }

  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*K, 0.25);
    z_tmp[j] = _z[j];
  }

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


  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it

  // scratch vectors to update parameters
  double sumdens, cprob, sumy=0.0, ssq=0.0, summu=0.0;
  double mstar, vstar, astar, bstar, max_dens, ns, os;
  double p0, an, ao, log_ao_density, log_an_density, llo,lln,llr,uu;
  double scr1[*K], dens[*ndens_y], ao_vec[*K], an_vec[*K];

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

  // Metropolis tunning parameter
  double csiga = 0.5;

//  Rprintf("Prior values: m0 = %.2f, s20 = %.2f, A = %.2f, A0 = %.2f, \n  %11s a0 = %.2f, b0 = %.2f, nu0 = %.2f, k0 = %.2f\n\n",
//             *m0, *s20, *A, *A0, *a0, *b0, *nu0, *k0);


  if(*mu_sigma_prior == 1){
    Rprintf("Prior values: m0 = %0.2f, s20 = %0.2f, A = %0.2f\n", *m0, *s20, *A);
  }

  if(*mu_sigma_prior == 2){
    Rprintf("Prior values: m0 = %0.2f, s20 = %0.2f, a0 = %0.2f, b0 = %0.2f\n", *m0, *s20, *a0, *b0);
  }

  if(*mu_sigma_prior == 3){
    Rprintf("Prior values: m0 = %0.2f, k0 = %0.2f, nu0 = %0.2f, s20 = %0.2f\n", *m0, *k0, *nu0, *s20);
  }
  
  ii = 0;

  for(i=0; i<*niter; i++){

//    Rprintf("i =========================================== %d\n", i);

    for(k=0; k<*K; k++){
      REAL(nk_vec)[k] = 0.0;
      nk_vec2[k]=0;
    }


    // 1) update zi - component labels
    for(j=0; j<*n; j++){
//      Rprintf("j = %f\n", j);
      max_dens = dnorm(y[j], _mu[0], sqrt(_sigma2[0]), 1) + log(_w[0]);
      for(k=0; k<*K; k++){
        dens[k] = dnorm(y[j], _mu[k], sqrt(_sigma2[k]), 1) + log(_w[k]);
        if(dens[k] > max_dens) max_dens = dens[k];
      }
//      RprintVecAsMat("dens", dens, 1, *K);
//      Rprintf("max_dens = %f\n", max_dens);
      sumdens = 0.0;
      for(k=0; k<*K; k++){
        dens[k] = exp(dens[k] - max_dens);
        sumdens = sumdens + dens[k];
      }
//      RprintVecAsMat("dens", dens, 1, *K);
      for(k=0; k<*K; k++){
        scr1[k] = dens[k]/sumdens;
      }
//      RprintVecAsMat("scr1", scr1, 1, *K);

	  uu = runif(0.0,1.0);
	  cprob= 0.0;
      for(k = 0; k < *K; k++){
	    cprob = cprob + scr1[k];
		if (uu < cprob){
		  _z[j] = k+1;
		  z_tmp[j] = k+1;
		  REAL(nk_vec)[k] = REAL(nk_vec)[k] + 1.0;
		  nk_vec2[k] = nk_vec2[k] + 1;
		  break;
		}
	  }
    }

    // Should I relabel here and make the first component always that which has the most
    // number of members?
//    RprintIVecAsMat("z = ", _z, 1, *n);
//    RprintIVecAsMat("z_tmp", z_tmp, 1, *n);
//    RprintVecAsMat("mu", _mu, 1, *K);
//    RprintVecAsMat("sigma2", _sigma2, 1, *K);
    if(*alpha_prior_type==2){ // centered prior need to reorder based on cluster size
      R_orderVector(nk_vecOrd, *K, Rf_lang1(nk_vec), TRUE, TRUE);
      for(k=0; k<*K; k++){

        for(j=0; j<*n; j++){
          if(z_tmp[j] == (nk_vecOrd[k]+1)) _z[j] = k+1;
        }
        _sigma2[k] = sigma2_tmp[nk_vecOrd[k]];
        _mu[k] = mu_tmp[nk_vecOrd[k]];
        _w[k] = w_tmp[nk_vecOrd[k]];
        alpha_loc[k] = alpha_loc_tmp[nk_vecOrd[k]];
        
      }
//      rsort_with_index()
    }




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
//      Rprintf("alphastar = %f\n", alphastar[k]);
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          alphastar[k] = alphastar[k] + 1;
        }
      }
    }

    ran_dirich(alphastar, *K, scr1, _w);
//    RprintVecAsMat("_w", _w, 1, *K);
    for(k=0; k<*K; k++){
      if(_w[k] < 1e-323) _w[k] = 1e-323;
      w_tmp[k] = _w[k];
      alpha_loc_tmp[k] = alpha_loc[k];
    }
//    RprintVecAsMat("_w", _w, 1, *K);



 //    Rprintf("alpha = %f\n", _alpha);



    // 3) update mu
    for(k=0; k<*K; k++){
    
      sumy = 0.0;
      nk = 0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          sumy = sumy + y[j];
          nk = nk+1;
        }
      }
      
      // This is the independent prior for muk,sigma2k so 
      // muk ~ N(mu0, sigma20) and sigma2k ~ UN or sigma2k ~ IG
      if(*mu_sigma_prior != 3 ){
        vstar = 1.0/((nk/_sigma2[k]) + 1.0/(_sigma20));
        mstar = vstar*((1.0/_sigma2[k])*sumy + (_mu0)/(_sigma20));
      }

      // This is for NIG prior on (muk|sigma^2k)(sigma^2k)
      if(*mu_sigma_prior == 3){   
        vstar = _sigma2[k]/(nk + *k0);
        mstar = (sumy + (*k0)*(_mu0))/(nk + *k0);
      }
      
      _mu[k] = rnorm(mstar, sqrt(vstar));
      mu_tmp[k] = _mu[k];
    }
   
//    RprintVecAsMat("muk = ", _mu, 1, *K);
   
   

    // 4) update sigma2
    for(k=0; k<*K; k++){
    
      // This is the independent prior for muk,sigma2k and sigmak ~ UN(0,A)
      if(*mu_sigma_prior == 1){
    
        os = sqrt(_sigma2[k]);
        ns = rnorm(os, 0.5);
        if(ns > 0){
          llo=0.0, lln=0.0;
          for(j=0; j<*n; j++){
            if(_z[j] == (k+1)){
              llo = llo + dnorm(y[j],_mu[k], os, 1);
              lln = lln + dnorm(y[j],_mu[k], ns, 1);
            }
          }
          llo = llo + dunif(os, 0, *A, 1);
          lln = lln + dunif(ns, 0, *A, 1);

          llr = lln - llo;
          uu = runif(0,1);
          if(llr > log(uu)) _sigma2[k] = ns*ns;

        }
      }
        
      // This is the independent prior for muk,sigma2k and sigmak ~ IG(a0,b0)
      if(*mu_sigma_prior == 2){
        ssq=0.0, nk=0;
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            ssq = ssq + (y[j]-_mu[k])*(y[j]-_mu[k]);
            nk = nk + 1;
          }
        }
        astar = 0.5*(nk) + *a0;
        bstar = 0.5*ssq + *b0;
        
        _sigma2[k] = 1/rgamma(astar, 1/bstar);
      }


      // This is the dependent prior for (muk,sigma2k) ~ NIG(m0,k0,v0,s20)
      if(*mu_sigma_prior == 3){
        ssq=0.0, nk=0;
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            ssq = ssq + (y[j]-_mu[k])*(y[j]-_mu[k]);
            nk = nk + 1;
          }
        }
        astar = 0.5*(nk + 1 + *nu0);
        bstar = 0.5*ssq + 0.5*(*k0)*(mu[k]- _mu0)*(mu[k]- _mu0) + 2*((*s20)/(*nu0));
        
        _sigma2[k] = 1/rgamma(astar, 1/bstar);
      }
      
      sigma2_tmp[k] = _sigma2[k];
    }
    
//    RprintVecAsMat("sigma2 = ", _sigma2, 1, *K);




//



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
//            ao_vec[k] = exp(ao);
//            an_vec[k] = exp(an);
              ao_vec[k] = ao;
              an_vec[k] = an;
            }

  //        Rprintf("exp(ao) = %f\n", exp(ao));
  //        Rprintf("exp(an) = %f\n", exp(an));
  //        RprintVecAsMat("w", _w, 1, *K);
  //        Rprintf("ddirich(_w, ao_vec, *K, 1) = %f\n", ddirich(_w, ao_vec, *K, 1));
  //        Rprintf("ddirich(_w, an_vec, *K, 1) = %f\n", ddirich(_w, an_vec, *K, 1));
          // add ao and an for the jacobian.
//          llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density) + ao;
//          lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density) + an;

            llo = ddirich(_w, ao_vec, *K, 1) + log(log_ao_density);
            lln = ddirich(_w, an_vec, *K, 1) + log(log_an_density);
  //        Rprintf("llo = %f\n", llo);
  //        Rprintf("lln = %f\n", lln);

            llr = lln - llo;

            uu = runif(0,1);

            if(llr > log(uu)){
//              _alpha = exp(an);
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
//          llo = ddirich(_w, ao_vec, *K, 1) + dgamma(ao, *a_gam, 1.0/((*a_gam)*(*K)), 1);
//          lln = ddirich(_w, an_vec, *K, 1) + dgamma(an, *a_gam, 1.0/((*a_gam)*(*K)), 1);
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
//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);



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

        if(*alpha_prior_dist==2){// uniform prior
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
    // keep iterates
    if((i > (*nburn-1)) & ((i) % *nthin ==0)){
//      Rprintf("ii = %d\n", ii);
      for(k=0; k<*K; k++){
        mu[ii + nout*k] = _mu[k];
        sigma2[ii + nout*k] = _sigma2[k];
        w[ii + nout*k] = _w[k];
      }
      for(j=0; j<*n; j++){
        z[ii + nout*j] = _z[j];
      }
      mu0[ii] = _mu0;
      sigma20[ii] = _sigma20;
      alpha[ii] = _alpha;
      alpha1[ii] = _alpha1;
      alpha2[ii] = _alpha2;


      if(*ndens_y > 0){
        for(j=0; j<*ndens_y; j++){
          dens[j] = 0.0;
          for(k=0; k<*K; k++){
            dens[j] = dens[j] + _w[k]*dnorm(y_grid[j], _mu[k], sqrt(_sigma2[k]), 0);
          }
          density_y[ii + nout*j] = dens[j];
        }
      }


      ii = ii + 1;
    }


  }
  
  UNPROTECT(1); // seems like you forgot this
}






SEXP IFMM_GIBBS(SEXP y, SEXP n, SEXP K, SEXP alpha_prior_type, SEXP alpha_prior_dist,
                SEXP mu_sigma_prior, SEXP U,
                SEXP basemodel, 
                SEXP alpha1_val, SEXP alpha2_val,
                SEXP update_alpha1, SEXP update_alpha2,
                SEXP y_grid, SEXP ndens_y, SEXP hierarchy,
                SEXP alpha_grid, SEXP alpha_density, SEXP ndens_alpha,
                SEXP alpha1_grid, SEXP alpha1_density, SEXP ndens_alpha1,
                SEXP alpha2_grid, SEXP alpha2_density, SEXP ndens_alpha2,
                SEXP m0, SEXP s20, SEXP A, SEXP A0,
                SEXP a0, SEXP b0, SEXP nu0, SEXP k0,
                SEXP a_gam, SEXP b_gam, SEXP a1_gam, SEXP b1_gam, 
                SEXP a2_gam, SEXP b2_gam,
 	            SEXP niter, SEXP nburn, SEXP nthin){


  int nprot = 0;

  int _n = asInteger(n);
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
  double _m0 = asReal(m0);
  double _s20 = asReal(s20);
  double _A = asReal(A);
  double _A0 = asReal(A0);
  double _a0 = asReal(a0);
  double _b0 = asReal(b0);
  double _nu0 = asReal(nu0);
  double _k0 = asReal(k0);
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
  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _K)); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++;
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++;
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++;
  SEXP ALPHA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA1 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP MU0 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP SIGMA20 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP DENSITY_Y = PROTECT(allocMatrix(REALSXP, nout, (_ndens_y))); nprot++;


  double *MUout, *SIGMA2out, *Wout, *MU0out, *SIGMA20out;
  double  *ALPHAout,  *ALPHA1out,  *ALPHA2out, *DENSITY_Yout;
  int *Zout;
  MUout = REAL(MU);
  SIGMA2out = REAL(SIGMA2);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  MU0out = REAL(MU0);
  SIGMA20out = REAL(SIGMA20);
  ALPHAout = REAL(ALPHA);
  ALPHA1out = REAL(ALPHA1);
  ALPHA2out = REAL(ALPHA2);
  DENSITY_Yout = REAL(DENSITY_Y);

  GetRNGstate();

  ifmm_gibbs(REAL(y), &_n, &_K, &_alpha_prior_type, &_alpha_prior_dist, 
              &_mu_sigma_prior, &_U,
              &_basemodel, 
              &_alpha1_val, &_alpha2_val,
              &_update_alpha1, &_update_alpha2,
              REAL(y_grid), &_ndens_y, &_hierarchy,
              REAL(alpha_grid), REAL(alpha_density), &_ndens_alpha,
              REAL(alpha1_grid), REAL(alpha1_density), &_ndens_alpha1,
              REAL(alpha2_grid), REAL(alpha2_density), &_ndens_alpha2,
              &_m0, &_s20, &_A, &_A0, &_a0, &_b0, &_nu0, &_k0, 
              &_a_gam, &_b_gam, &_a1_gam, &_b1_gam, &_a2_gam, &_b2_gam,
              &_niter, &_nburn, &_nthin,
              MUout, SIGMA2out, Wout, Zout, ALPHAout, ALPHA1out, ALPHA2out,
              MU0out, SIGMA20out, DENSITY_Yout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 10)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
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
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
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



