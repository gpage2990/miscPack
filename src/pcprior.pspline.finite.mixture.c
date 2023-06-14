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
// y - data vector (long vector of all subjects curves)
// Xmat - matrix of bspline basis functions (this is nobs*nb if balanced)
// n - number of subjects (or curves)
// nobs - vector containing number of observed points for each subject
// S - smoothing matrix for the splines
// nb - number of basis functions (dimension of S)
// K - number of mixture components
// alpha_prior_type - identifies if sparse or centered prior is employed
// alpha_prior_dist - identifies if gamma or pc prior is used (only on alpha1 currently)
// U - the "centered" value for 
// basemodel - identifies the base model used in the PC prior specification (currently set at sparse)
// alpha0 - specifies dirichlet distribution that is the base model
// alpha_grid - values of alpha that were used to evaluate the density of alpha
// alpha_density - density values of a corresponding to alpha-grid values
// ndens_alpha - number of values in alpha-grid
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
// a - prior shape of the variance from mixture components
// b - prior rate of the variance from mixture components
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
// 

static void pspline_ifmm_gibbs(double *y, 
                               double *Xmat, 
                               int *n, 
                               int *nobs,
                               double *S, 
                               int *nb, 
                               int *K,
                               int *alpha_prior_type, 
                               int *alpha_prior_dist, 
                               int *U,
                               int *basemodel, 
                               double *alpha1_val, 
                               double *alpha2_val,
                               int *update_alpha1, 
                               int *update_alpha2,
                               double *alpha_grid, 
                               double *alpha_density, 
                               int *ndens_alpha,
                               double *alpha1_grid, 
                               double *alpha1_density, 
                               int *ndens_alpha1,
                               double *alpha2_grid, 
                               double *alpha2_density, 
                               int *ndens_alpha2,                               
                               double *A, 
                               double *e_tau, 
                               double *e_omega, 
                               double *a_gam, 
                               double *b_gam,
                               double *a1_gam, 
                               double *b1_gam,
                               double *a2_gam, 
                               double *b2_gam,
                               int *niter, 
                               int *nburn, 
                               int *nthin,
                               double *beta, 
                               double *sigma2, 
                               double *theta, 
                               double *lambda2,
                               double *w, 
                               int *z, 
                               double *alpha,
                               double *alpha1,
                               double *alpha2, 
                               double *tau2, 
                               double *line_fit){


  Rprintf("basemodel = %d\n", *basemodel);
  Rprintf("n = %d\n", *n);
  Rprintf("K = %d\n", *K);
  Rprintf("nobs = %d\n", *nobs);
  
  if(*alpha_prior_type==2) Rprintf("U=%d\n", *U);
  int i, j, k, d, s, ss, t, ii, jj, nk, csobs;
  int nout = (*niter - *nburn)/(*nthin);
  double alphastar[*K], astar, bstar;
  
  double _beta[(*n)*(*nb)], _sigma2[*n];
  double _theta[(*K)*(*nb)],  _w[*K], _tau2[*K], _lambda2[*K]; 
  double _alpha, _alpha1, _alpha2;
  double theta_tmp[(*K)*(*nb)], lambda2_tmp[*K], w_tmp[*K];
  int _z[*n], nk_vec2[*K], z_tmp[*n];
 
  int *nk_vecOrd = (int *) R_alloc(*K, sizeof(int));
  SEXP nk_vec = PROTECT(Rf_allocVector(REALSXP, *K));
 
  
  // Initialize variables 
  _alpha = 0.5;
  _alpha1 = *alpha1_val;
  _alpha2 = *alpha2_val;

  for(k = 0; k < *K; k++){ 
  
    _tau2[k] = rgamma(1,1);
    
    _w[k] = 1/(*K);
    w_tmp[k] = _w[k];
    
    _lambda2[k] = 0.5;
    lambda2_tmp[k] = _lambda2[k];
    
    for(s = 0; s < *nb; s++){
      _theta[k*(*nb)+s] = rnorm(0,1);
      theta_tmp[k*(*nb)+s] = _theta[k*(*nb)+s];
    }
  }
  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*K, 0.25);
    _sigma2[j] = (0.5*(*A))*(0.5*(*A));
    for(s = 0; s < *nb; s++){
      _beta[j*(*nb)+s] = rnorm(0,1);
    }
  }


//  RprintVecAsMat("S", S, *nb, *nb);
                   
  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it
  
  // scratch vectors to update parameters
  double sumdens, cprob, ssq=0.0, xb=0.0, ld;
  //double astar, bstar;
  double max_dens;
  double p0, an, ao, ln, lo, to, tn, os, ns; 
  double log_ao_density, log_an_density, llo,lln,llr,uu;
  double dens[*K], ao_vec[*K], an_vec[*K], scr1[*K], scr2[*nb], scr3[*nb];
  double Mstar[*nb], Sstar[(*nb)*(*nb)], y_tmp[*nobs], sumb[*nb], sumXy[*nb];
  
//  RprintVecAsMat("alpha_grid",alpha_grid, 1, *ndens_alpha);
//  RprintVecAsMat("alpha_density",alpha_density, 1, *ndens_alpha);

  double *tX = R_VectorInit((*nb)*(*nobs),0.0);
  double *XtX = R_Vector(((*nb)*(*nb)));
  double *Xy = R_Vector((*nb));
  double *Xb = R_Vector(*nobs);
  double *outrmvnorm = R_Vector((*nb));
  
  
  mat_transpose(Xmat, tX, *nobs, (*nb));
  matrix_product(tX, Xmat, XtX, (*nb), (*nb), *nobs);
  
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
  double csiga = 0.5, csigl=0.5, csigt=0.5, csigs=0.0001;

  ii = 0;



  for(i=0; i<*niter; i++){
    if((i+1) % 10000 == 0){
      Rprintf("i =========================================== %d\n", i+1);
    }

    // 1) update w
    for(k=0; k<*K; k++){
      if(*alpha_prior_type == 1){//symmetric dirichlet
        alphastar[k] = _alpha;
      }
      if(*alpha_prior_type == 2){//asymmetric dirichlet
        if(k<(*U)){
          alphastar[k] = _alpha1;
        }
        if(k>=(*U)){
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
      REAL(nk_vec)[k] = 0.0;
      nk_vec2[k]=0;
      w_tmp[k] = _w[k];
    }
//  RprintVecAsMat("_w", _w, 1, *K);
  

    // 2) update zi - component labels
	csobs=0;
    for(j=0; j<*n; j++){
      
      for(k=0; k<*K; k++) dens[k] = 0.0;
      for(k=0; k<*K; k++){
//      Rprintf("k = %d\n", k);
	    for(s = 0; s < *nb; s++){
	      dens[k] = dens[k]+dnorm(_beta[j*(*nb)+s], _theta[k*(*nb)+s], sqrt(_lambda2[k]), 1);
	    }
        dens[k] = dens[k] + log(_w[k]);
      }
      
      max_dens = dens[0];
      for(k=1; k<*K; k++){
        if(max_dens < dens[k]) max_dens = dens[k];
      }
      
      sumdens = 0.0;
      for(k=0; k<*K; k++){
        dens[k] = exp(dens[k] - max_dens);
        sumdens = sumdens + dens[k];
      }

      for(k=0; k<*K; k++){
        scr1[k] = dens[k]/sumdens;
      }

//      if(j == 0){
//        Rprintf("j = %d\n", j);
//        RprintVecAsMat("prob",scr1, 1, *K);
//      }
	  uu = runif(0.0,1.0);
	  cprob= 0.0;
      for(k = 0; k < *K; k++){
	    cprob = cprob + scr1[k];
		if (uu < cprob){
		  _z[j] = k+1;
		  break;
		}
	  }
    }
//    RprintIVecAsMat("_z", _z, 1, *n);
   

    if(*alpha_prior_type==2){ // centered prior need to reorder based on cluster size
      R_orderVector(nk_vecOrd, *K, Rf_lang1(nk_vec), TRUE, TRUE);
      for(k=0; k<*K; k++){

        for(j=0; j<*n; j++){
          if(z_tmp[j] == (nk_vecOrd[k]+1)) _z[j] = k+1;
        }
        _lambda2[k] = lambda2_tmp[nk_vecOrd[k]];
        _w[k] = w_tmp[nk_vecOrd[k]];

	    for(s = 0; s < *nb; s++){
	      _theta[k*(*nb)+s] = theta_tmp[nk_vecOrd[k]*(*nb)+s];
	    }
      }
//      rsort_with_index()
    }




    // 3) update alpha
    // This code is for the sparse model for K+
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
      // Use gamma prior rather than pc prior
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


    // This code is for the centered prior on K+
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
              if(k<(*U)){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(k>=(*U)){
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
              if(k<(*U)){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(k>=(*U)){
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
              if(k<(*U)){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(k>=(*U)){
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
              if(k<(*U)){
                ao_vec[k] = ao;
                an_vec[k] = an;
              }
              if(k>=(*U)){
                ao_vec[k] = _alpha2;
                an_vec[k] = _alpha2;
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
    
    
//    Rprintf("alpha = %f\n", _alpha);


    // 4) update beta and sigma2 within same loop
	csobs=0;
    for(j=0; j<*n; j++){

//      Rprintf("j = %d\n", j);
      
      
      // 4a) update beta
	  for(t=0; t<*nobs; t++){
	    y_tmp[t] = y[csobs+t];
	  } 
      csobs = csobs + (*nobs);
      
      matrix_product(tX, y_tmp, Xy, *nb, 1, *nobs);
	  
	  for(s = 0; s < (*nb); s++){

	    sumXy[s] = (1.0/_sigma2[j])*Xy[s] + 
	              (1.0/_lambda2[(_z[j]-1)])*_theta[(_z[j]-1)*(*nb)+s];

	    for(ss = 0; ss < (*nb); ss++){
		  Sstar[s*(*nb)+ss] = (1.0/_sigma2[j])*XtX[s*(*nb)+ss];

		  if(s == ss){
		    Sstar[s*(*nb)+ss] = (1.0/_sigma2[j])*XtX[s*(*nb)+ss]+
		                        (1.0/_lambda2[(_z[j]-1)]);
		  }
		}
	  }


	  cholesky(Sstar, (*nb), &ld);
	  inverse_from_cholesky(Sstar, scr2, scr3, (*nb));
//      RprintVecAsMat("Sstar", Sstar, *nb, *nb);

    
//      RprintVecAsMat("Xy", Xy, 1, *nb);

	  matrix_product(Sstar, sumXy, Mstar, (*nb), 1, (*nb));

//      RprintVecAsMat("Mstar", Mstar, 1, *nb);


	  cholesky(Sstar, (*nb) , &ld);


	  ran_mvnorm(Mstar, Sstar, (*nb), scr2, outrmvnorm);

//      RprintVecAsMat("beta", outrmvnorm, 1, *nb);

	  for(s=0; s < (*nb); s++){
	    _beta[j*(*nb) + s] = outrmvnorm[s];
	  }

      matrix_product(Xmat, outrmvnorm, Xb, *nobs, 1, (*nb));
	  




      // 4b) update sigma2
      os = sqrt(_sigma2[j]);
      ns = rnorm(os, csigs);
      if(ns > 0){
        llo = 0; lln = 0;
        for(jj = 0; jj <*nobs; jj++){
          llo = llo + dnorm(y_tmp[jj], Xb[jj], os, 1);
          lln = lln + dnorm(y_tmp[jj], Xb[jj], ns, 1);
        }
        llo = llo + dunif(os, 0.0, *A, 1);
        lln = lln + dunif(ns, 0.0, *A, 1);
        
        llr = lln - llo;
        uu = runif(0,1);
        if(llr > log(uu)) _sigma2[j] = ns*ns;
      }

/*
      ssq=0.0;
      for(jj = 0; jj < *nobs; jj++){
      	ssq = ssq + (y_tmp[jj] - Xb[jj])*(y_tmp[jj] - Xb[jj]);
      }
            
      astar = 0.5*(*nobs) + 3.0;
      bstar = 0.5*ssq + 2.0;
      
      //bstar is rate and rgamma requires scale hence inverse
      
     _sigma2[j] = 1/rgamma(astar, 1/bstar);
*/      
//      Rprintf("sigma2 = %f\n", _sigma2[j]);
	  
    }



//    // 5) update sigma2    
//    for(k=0; k<*K; k++){
//      ssq=0.0, nk=0, csobs=0;
//      for(j=0; j<*n; j++){
//        if(_z[j] == (k+1)){
//          for(t = 0; t<*nobs; t++){
//      	    xb = 0.0;
//	        for(s = 0; s < *nb; s++){
//	          xb = xb + _beta[k*(*nb)+s]*Xmat[t*(*nb)+s];
//	        }
//            ssq = ssq + (y[csobs+t]-xb)*(y[csobs+t]-xb);
//            nk = nk + 1;
//          }
//        }
//        csobs = csobs + (*nobs);
//      }
//      astar = 0.5*(*nobs)*(nk) + *a;
//      bstar = 0.5*ssq + *b;
//            
//      _sigma2[k] = 1/rgamma(astar, 1/bstar);
////      _sigma2[k] = 0.1*0.1;
//    }

    
    
    // 6) update theta and lambda in the same for-loop
    for(k = 0; k < *K; k++){
//      Rprintf("k = %d\n", k);
      // 6a) update lambda2
//      lo = sqrt(_lambda2[k]);
//      ln = rnorm(lo,csigl);

      lo = _lambda2[k];
      ln = rnorm(lo,csigl);
      
//      Rprintf("lo = %f\n", lo);
//      Rprintf("ln = %f\n", ln);
    
      if(ln > 0){
        llo=0.0, lln=0.0, nk=0;
        for(j=0; j<*n; j++){
          if(_z[j] == k+1){
            for(s=0; s < (*nb); s++){
              llo = llo + dnorm(_beta[j*(*nb)+s], _theta[k*(*nb)+s], sqrt(lo), 1);
              lln = lln + dnorm(_beta[j*(*nb)+s], _theta[k*(*nb)+s], sqrt(ln), 1);
            }
            nk = nk + 1;
          }
        }
    
//    	llo = llo + dunif(lo, 0.0, *A, 1);
//    	lln = lln + dunif(ln, 0.0, *A, 1);

        // Note, this gumbel is parametrized for precision
        llo = llo + log(0.5*(*e_omega)) - 1.5*log(1.0/lo) - (*e_omega)*sqrt(lo);
        lln = lln + log(0.5*(*e_omega)) - 1.5*log(1.0/ln) - (*e_omega)*sqrt(ln);
    
//        Rprintf("llo = %f\n", llo);
//        Rprintf("lln = %f\n", lln);
    	llr = lln - llo;
    	uu = runif(0.0,1.0);
    
    	if(log(uu) < llr) _lambda2[k] = ln;
//    	if(log(uu) < llr) _lambda2[k] = ln;
    	lambda2_tmp[k] = _lambda2[k];
      }
//      Rprintf("lambda2 = %f\n", _lambda2[k]);
      
//      Rprintf("nk = %d\n", nk);

      // 6b) update theta    
      for(s = 0; s < (*nb); s++){
        for(ss = 0; ss < (*nb); ss++){
    	  Sstar[s*(*nb)+ss] = (1/_tau2[k])*S[s*(*nb)+ss];
    	  if(s == ss){ 
    	    Sstar[s*(*nb)+ss] = ((double) nk/(_lambda2[k])) +
    	 							  (1/_tau2[k])*S[s*(*nb)+ss];
    	  }
    	}
    	sumb[s] = 0.0;
      }
    
      cholesky(Sstar, (*nb), &ld);
      inverse_from_cholesky(Sstar, scr2, scr3, (*nb));
    
      for(j=0; j<*n; j++){
        if(_z[j] == k+1){
          for(s=0; s<(*nb); s++){
            sumb[s] = sumb[s] + 
      	               (1.0/_lambda2[k])*_beta[j*(*nb)+s];
          }
        }
      }
      matrix_product(Sstar, sumb, Mstar, (*nb), 1, (*nb));
    
      cholesky(Sstar, (*nb) , &ld);
    
      ran_mvnorm(Mstar, Sstar, (*nb), scr2, outrmvnorm);
      
//      RprintVecAsMat("theta", outrmvnorm, 1, *nb);
    
      for(s=0; s<(*nb); s++){
        _theta[k*(*nb)+s] = outrmvnorm[s];
      }
    }
    
    // 7) update tau2 using a PC prior    
    for(k = 0; k<*K; k++){
      to = _tau2[k];
      tn = rnorm(to, csigt);
      if(tn > 0){

        for(s=0; s<*nb; s++){
          scr2[s] = _theta[k*(*nb) + s];
        }
	    ssq = quform(scr2,S,(*nb));
	    
//        Rprintf("ssq = %f\n", ssq);

        llo =  -0.5*(*nb)*log(to) - (0.5/to)*ssq;
        lln =  -0.5*(*nb)*log(tn) - (0.5/tn)*ssq;

        // Note, this gumbel is parametrized for precision
        llo = llo + log(0.5*(*e_tau)) - 1.5*log(1.0/to) - (*e_tau)*sqrt(to);
        lln = lln + log(0.5*(*e_tau)) - 1.5*log(1.0/tn) - (*e_tau)*sqrt(tn);


    	llr = lln - llo;
    	uu = runif(0.0,1.0);

    	if(log(uu) < llr) _tau2[k] = tn;
      }
    }

//\theta/2 tau^{-3/2}exp(-theta/sqrt(tau))

    
//    // 7) update tau2 (smoothing parameter)    
//    ssq = 0.0;
//    for(k=0; k<*K; k++){
//      for(s=0; s<*nb; s++){
//        scr2[s] = _theta[k*(*nb) + s];
//      }
//	  ssq = quform(scr2,S,(*nb));
//
//      astar = 0.5*(*nb) + *a;
//      bstar = 0.5*ssq + *b;
//    
//      _tau2[k] = 1.0/rgamma(astar, 1.0/bstar);
//    }
//    RprintVecAsMat("tau2", _tau2, 1, *K);
   
 

    // keep iterates
    if((i > (*nburn-1)) & ((i) % *nthin ==0)){
//      Rprintf("ii = %d\n", ii);
      ss = 0;
      for(k=0; k<*K; k++){
        tau2[ii + nout*k] = _tau2[k];
        w[ii + nout*k] = _w[k];
        lambda2[ii + nout*k] = _lambda2[k];
        for(s=0; s<*nb; s++){
          theta[ii + nout*ss] = _theta[k*(*nb)+s];
          ss = ss+1;
        }
      }
      ss=0;
      for(j=0; j<*n; j++){
        sigma2[ii + nout*j] = _sigma2[j];
        for(s=0; s<*nb; s++){
          beta[ii + nout*ss] = _beta[j*(*nb)+s];
          ss = ss+1;
        }
      }
      ss = 0;
      for(j=0; j<*n; j++){
        z[ii + nout*j] = _z[j];

        for(t = 0; t<*nobs; t++){
          xb = 0.0;
	      for(s = 0; s < *nb; s++){
	        xb = xb + _beta[j*(*nb)+s]*Xmat[t*(*nb)+s];	        
	      }
	      line_fit[ii + nout*ss] = xb;
	      ss = ss+1;
        }
      }
      alpha[ii] = _alpha;
      alpha1[ii] = _alpha1;
      alpha2[ii] = _alpha2;
      
//    	matrix_product(Htheta, outrmvnorm, Ht, *nrHt, 1, (*nb));
//    	for(t = 0; t < *nrHt; t++){
//    		fgprime_iter[k*(*nrHt) + t] = Ht[t];
//    	}
    
 
      
      ii = ii + 1;
    }
        
  }
  UNPROTECT(1);
}



                               
SEXP PSPLINE_IFMM_GIBBS(SEXP y, SEXP Xmat, SEXP n, SEXP nobs,            //4
                SEXP S, SEXP nb, SEXP K,                                 //3
                SEXP alpha_prior_type,                                   //1
                SEXP alpha_prior_dist,            						 //1
                SEXP U,               									 //1
                SEXP basemodel, 										 //1		
                SEXP alpha1_val, SEXP alpha2_val,						 //2
                SEXP update_alpha1, SEXP update_alpha2,                  //2
                SEXP alpha_grid, SEXP alpha_density, SEXP ndens_alpha,   //3
                SEXP alpha1_grid, SEXP alpha1_density, SEXP ndens_alpha1,//3
                SEXP alpha2_grid, SEXP alpha2_density, SEXP ndens_alpha2,//3
                SEXP A, SEXP e_tau, SEXP e_omega, 						 //3
                SEXP a_gam, SEXP b_gam,									 //2
                SEXP a1_gam, SEXP b1_gam, 								 //2	
                SEXP a2_gam, SEXP b2_gam, 								 //2
	            SEXP niter, SEXP nburn, SEXP nthin){					 //3	
																		 //36 - total	

  int nprot = 0;	     
         
  int _n = asInteger(n);
  int _nobs = asInteger(nobs);
  int _nb = asInteger(nb);
  int _K = asInteger(K);
  int _alpha_prior_type = asInteger(alpha_prior_type);
  int _alpha_prior_dist = asInteger(alpha_prior_dist);
  int _U = asInteger(U);
  int _basemodel = asInteger(basemodel);
  int _update_alpha1 = asInteger(update_alpha1);
  int _update_alpha2 = asInteger(update_alpha2);
  int _ndens_alpha = asInteger(ndens_alpha);
  int _ndens_alpha1 = asInteger(ndens_alpha1);
  int _ndens_alpha2 = asInteger(ndens_alpha2);
  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);
  double _alpha1_val = asReal(alpha1_val);
  double _alpha2_val = asReal(alpha2_val);
  double _A = asReal(A);
  double _e_tau = asReal(e_tau);
  double _e_omega = asReal(e_omega);
  double _a_gam = asReal(a_gam);
  double _b_gam = asReal(b_gam);
  double _a1_gam = asReal(a1_gam);
  double _b1_gam = asReal(b1_gam);
  double _a2_gam = asReal(a2_gam);
  double _b2_gam = asReal(b2_gam);

  double nout = (_niter-_nburn)/_nthin;

  
  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  Xmat = PROTECT(coerceVector(Xmat, REALSXP)); nprot++;
  S = PROTECT(coerceVector(S, REALSXP)); nprot++;
  alpha_grid = PROTECT(coerceVector(alpha_grid, REALSXP)); nprot++;
  alpha1_grid = PROTECT(coerceVector(alpha1_grid, REALSXP)); nprot++;
  alpha2_grid = PROTECT(coerceVector(alpha2_grid, REALSXP)); nprot++;
  alpha_density = PROTECT(coerceVector(alpha_density, REALSXP)); nprot++;
  alpha1_density = PROTECT(coerceVector(alpha1_density, REALSXP)); nprot++;
  alpha2_density = PROTECT(coerceVector(alpha2_density, REALSXP)); nprot++;
  SEXP BETA = PROTECT(allocMatrix(REALSXP, nout, (_n)*(_nb))); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_n))); nprot++; 
  SEXP THETA = PROTECT(allocMatrix(REALSXP, nout, (_K)*(_nb))); nprot++;
  SEXP LAMBDA2 = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++;
  SEXP TAU2 = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++;
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_K))); nprot++; 
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, (_n))); nprot++; 
  SEXP ALPHA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA1 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP ALPHA2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP LINE_FIT = PROTECT(allocMatrix(REALSXP, nout, (_n)*(_nobs))); nprot++;


  double *BETAout, *SIGMA2out, *THETAout, *LAMBDA2out;
  double *Wout, *TAU2out, *ALPHAout,  *ALPHA1out,  *ALPHA2out, *LINE_FITout;
  int *Zout;
  BETAout = REAL(BETA);
  SIGMA2out = REAL(SIGMA2);
  THETAout = REAL(THETA);
  LAMBDA2out = REAL(LAMBDA2);
  TAU2out = REAL(TAU2);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  ALPHAout = REAL(ALPHA);
  ALPHA1out = REAL(ALPHA1);
  ALPHA2out = REAL(ALPHA2);
  LINE_FITout = REAL(LINE_FIT);

  GetRNGstate();



  pspline_ifmm_gibbs(REAL(y), REAL(Xmat), &_n, &_nobs, 
                     REAL(S), &_nb, &_K, 
                     &_alpha_prior_type, 
                     &_alpha_prior_dist, 
                     &_U,
                     &_basemodel, 
                     &_alpha1_val, &_alpha2_val,
                     &_update_alpha1, &_update_alpha2,
                     REAL(alpha_grid), REAL(alpha_density), &_ndens_alpha,
                     REAL(alpha1_grid), REAL(alpha1_density), &_ndens_alpha1,
                     REAL(alpha2_grid), REAL(alpha2_density), &_ndens_alpha2,
                     &_A, &_e_tau, &_e_omega, 
                     &_a_gam, &_b_gam, 
                     &_a1_gam, &_b1_gam, 
                     &_a2_gam, &_b2_gam,
                     &_niter, &_nburn, &_nthin, 
                     BETAout, SIGMA2out, THETAout, LAMBDA2out, 
                     Wout, Zout, ALPHAout, ALPHA1out, ALPHA2out, TAU2out,
                     LINE_FITout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 11)); nprot++;
  SET_VECTOR_ELT(ans, 0, BETA);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, THETA);
  SET_VECTOR_ELT(ans, 3, LAMBDA2);
  SET_VECTOR_ELT(ans, 4, W);
  SET_VECTOR_ELT(ans, 5, Z);
  SET_VECTOR_ELT(ans, 6, ALPHA);
  SET_VECTOR_ELT(ans, 7, ALPHA1);
  SET_VECTOR_ELT(ans, 8, ALPHA2);
  SET_VECTOR_ELT(ans, 9, TAU2);
  SET_VECTOR_ELT(ans, 10, LINE_FIT);


  SEXP nm = allocVector(STRSXP, 11);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("beta"));
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
  SET_STRING_ELT(nm, 2, mkChar("theta"));
  SET_STRING_ELT(nm, 3, mkChar("lambda2"));
  SET_STRING_ELT(nm, 4, mkChar("w"));
  SET_STRING_ELT(nm, 5, mkChar("z"));
  SET_STRING_ELT(nm, 6, mkChar("alpha"));
  SET_STRING_ELT(nm, 7, mkChar("alpha1"));
  SET_STRING_ELT(nm, 8, mkChar("alpha2"));
  SET_STRING_ELT(nm, 9, mkChar("tau2"));
  SET_STRING_ELT(nm, 10, mkChar("line_fit"));

  UNPROTECT(nprot);
  return(ans);
 
}



