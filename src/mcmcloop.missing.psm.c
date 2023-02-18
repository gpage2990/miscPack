/****************************************************************************************
 * Copyright (c) 2014 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit PPMx models when covariate values are missing
 * The trick is to carry a label vector indicating when a covariate
 * is missing.
 *
 *
 ****************************************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
* nobs = number of response values
* ncon = number of continuous covariates
* ncat = number of categorical covariates
* Mcon = nobs x ncon logical vector indicating missing of continuous covariates.
* Mcat = nobs x ncat logical vector indicating missing of categorical covariates.
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.
* PPM = logical indicating if PPM or PPMx should be used.
* gvec = nobs x 1 integer indicating what beta group individual belongs (which subset of betas).
* nb = nobs x 1 integer that indicates dimension of beta.
* gvecp = npred x 1 integer indicating what beta group predictive individual belongs (which subset of betas).
* nbp = npred x 1 integer that indicates dimension of predictive beta.
* M = double indicating value of M associated with cohesion.
* gcontype = int indicating similarity function for continuous variable
* gcattype = int indicating similarity function for categorical variable
* y = nobs x 1 vector that contains response values
* Xcon = nobs x ncon contiguous double vector that contains continuous covariate values
* Xcat = nobs x ncat contiguous int vector that contains categorical covariate values
* npred = integer indicating number of out of sample predictions
* Xconp = npred x ncon matrix of continuous covariates to make predictions
* Xcatp = npred x ncat matrix of categorical covariates to make predictions
* Mconp = npred x ncon matrix of logicals indicating if continuous covariate is missing
* Mcatp = npred x ncat matrix of logicals indicating if categorical covariate is missing
* simparms = vector containing similarity functions that are ordered in the following way;
* calibrate = integer that indicates to calibration (here it is coarsened i.e. 2)
* consim = 1 if N-N continuous similarity 2 - if it is N-NIG
* alpha = scalar associated with tuning parameter
          if similarity = 3 then g(x_j) =  exp(-alpha var(x_j));

* Output:
* mu = nout x ncov matrix of MCMC iterates
* beta = nout x ncov matrix of MCMC iterates
* sig2 = nout vector of MCMC iterates
* mu0 = nout vector of MCMC iterates
* sig20 = nout vector of MCMC iterates
* Si = nout vector of MCMC iterates
* nclus = nout vector of MCMC iterates
* like = scalar with lpml value
* lpml = scalar with lpml value
* waic = scalar containing waic value
* ispred
* ppred
* ppredclass
*****************************************************************************************/


void mcmc_missing_psm(int *draws, int *burn, int *thin, int *nobs, int *ncon, int *ncat,
      int *Cvec, int *PPM, int *gvec, int *nb, int *gvecp,
      int *nbp, double *M, int *gcontype, int *gcattype,
      double *y, double *Xcon, int *Xcat, int *Mcon, int *Mcat,
      int *npred, double *Xconp, int *Xcatp,
			int *Mconp, int *Mcatp, int *calibrate,
			int*consim, double *alpha, double *simParms, double *modelPriors,
			double *mh,
			double *mu, double *beta, double *sig2, double *mu0, double *sig20,
			int *Si, int *nclus, double *like, double *WAIC, double *lpml,
			double *ispred, double *ppred, int *predclass,
			double *pred2, int *verbose){







	// i - MCMC index
	// ii - MCMC index for saving iterates
	// j - individual index
	// jj - second player index (for double for loops)
	// c - categorical variable index
	// p - number of covariates index
	// pp - second covariate index
	// k - cluster index
	// g - covariate configuration index (groups of covariates depending on missing)
	// t - obs per cluster indx
	// b - covariate index
	// bb - second covariate index

	int i, ii, j, jj, c, p, pp, k, g, t, b, bb;
	int nhtmp, nbtotal;
	int bcon=0, bcat=0;
	int ncov = *ncon + *ncat;

	Rprintf("nobs = %d\n", *nobs);
	Rprintf("ncon = %d\n", *ncon);
	Rprintf("ncat = %d\n", *ncat);
	Rprintf("ncov = %d\n", ncov);
	Rprintf("npred = %d\n", *npred);
	Rprintf("nb = %d\n", *nb);
//	RprintVecAsMat("Xcon", Xcon, *nobs, *ncon);
//	RprintIVecAsMat("Xcat", Xcat, *nobs, *ncat);
//	RprintIVecAsMat("Mcon", Mcon, *nobs, *ncon);
//	RprintIVecAsMat("Mcat", Mcat, *nobs, *ncat);
//	RprintIVecAsMat("nb",nb, 1, *nobs);
	if(*ncat>0) RprintIVecAsMat("Cvec",Cvec, 1, *ncat);

	int cumindx[*nobs];

	double max_C, nout, sumy=0.0,sumx, sumx2;;
	nout = (*draws-*burn)/(*thin);

	max_C = Cvec[0];

	for(p = 0; p < (*ncat); p++){
		if(max_C < Cvec[p]) max_C = Cvec[p];
	}
	if(*ncat == 0) max_C = 1.0;

	int ng;
	ng = 1;
	cumindx[0] = 0;
	nbtotal=0;
	for(j = 0; j < *nobs; j++){
		if(ng < gvec[j]) ng = gvec[j];
		cumindx[j+1] = cumindx[j] + nb[j];
		nbtotal=nbtotal+nb[j];
		sumy += y[j];
	}

	Rprintf("max_C = %f\n", max_C);
	Rprintf("ng = %d\n", ng);
	//Rprintf("sumy = %f\n", sumy);
	//RprintIVecAsMat("cumindx", cumindx, 1, *nobs);

	double *mnmle = R_Vector(*ncon);
	double *s2mle = R_Vector(*ncon);
	for(p = 0; p < *ncon; p++){
	  sumx = 0.0, sumx2=0.0;
	  for(j = 0; j < *nobs; j ++){
	    sumx = sumx + Xcon[j*(*ncon) + p];
	    sumx2 = sumx2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
	  }

	  mnmle[p] = sumx/((double) *nobs);
	  s2mle[p] = sumx2/((double) *nobs) - mnmle[p]*mnmle[p];
	}

	// RprintVecAsMat("mnmle", mnmle, 1, *ncon);
	// RprintVecAsMat("s2mle", s2mle, 1, *ncon);



	double *fullXmat = R_VectorInit((*nobs)*(ncov), 0.0);
	double *fullXmatp = R_VectorInit((*npred)*(ncov), 0.0);
	int fullM[(*nobs)*(ncov)];
	int fullMp[(*npred)*(ncov)];
//	Rprintf("ncov = %d\n", ncov);


	for(j = 0; j < *nobs; j++){
//		Rprintf("j = %d\n", j);
		bcon = 0;
		bcat = 0;
		for(b = 0; b < ncov; b++){
			if(b < *ncon){

				fullXmat[j*ncov+b] = Xcon[j*(*ncon)+bcon];
				fullM[j*ncov+b] = Mcon[j*(*ncon) + bcon];
				bcon = bcon + 1;
			}

			if(b >= *ncon){

				fullXmat[j*ncov+b] = Xcat[j*(*ncat)+bcat];
				fullM[j*ncov+b] = Mcat[j*(*ncat) + bcat];
				bcat = bcat + 1;
			}
		}
	}

//	RprintVecAsMat("fullXmat", fullXmat, *nobs, ncov);
//	RprintIVecAsMat("fullM", fullM, *nobs, ncov);

	for(pp = 0; pp < *npred; pp++){
		bcon = 0;
		bcat = 0;
		for(b = 0; b < ncov; b++){
			if(b < *ncon){

				fullXmatp[pp*ncov+b] = Xconp[pp*(*ncon)+bcon];
				fullMp[pp*ncov+b] = Mconp[pp*(*ncon) + bcon];
				bcon = bcon + 1;
			}

			if(b >= *ncon){

				fullXmatp[pp*ncov+b] = Xcatp[pp*(*ncat)+bcat];
				fullMp[pp*ncov+b] = Mcatp[pp*(*ncat) + bcat];
				bcat = bcat + 1;
			}
		}
	}

//	RprintVecAsMat("fullXmatp", fullXmatp, *nobs, ncov);
//	RprintIVecAsMat("fullMp", fullMp, *npred, ncov);








	// ===================================================================================
	//
	// Memory vectors to hold a single MCMC iterate
	//
	// ===================================================================================

	double mu0_iter = 0.0;
	double sig20_iter =1.0;
	double *beta_iter = R_VectorInit((*nobs)*(ncov),0.0);


	int nclus_iter = 0;

	int kmenos, iaux;

	int Si_iter[*nobs];
	int nh[*nobs];
	int nhc[(*nobs)*(*ncat)];


	double* xcontmp = R_Vector(*nobs);
	double* xcattmp = R_Vector(*nobs);


//	RprintVecAsMat("beta_iter", beta_iter, nb, *nobs);

	for(j = 0; j < *nobs; j++){
//		Si_iter[j] = j+1;
		Si_iter[j] = 1;
		nh[j] = 0;
	}

	for(j = 0; j < *nobs; j++){

		for(k = 0; k < *nobs; k++){

			if(Si_iter[j] == k+1) nh[k] = nh[k] + 1;
		}
	}

	for(j = 0; j < *nobs; j++){

		if(nh[j] > 0) nclus_iter = nclus_iter + 1;

	}


	kmenos = nclus_iter;

//	Rprintf("nclus_iter = %d\n", nclus_iter);
//	RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//	RprintIVecAsMat("nc", nc, 1, *nobs);
//	RprintIVecAsMat("nh", nh, 1, *nobs);



	double *sig2h = R_VectorInit(*nobs, 1.0);
	double *muh = R_VectorInit(*nobs, 0.0);




	// Should I compute the posterior predictive
	double *ppred_iter = R_Vector((*nobs)*(*nobs));
	double *pred2_iter = R_Vector((*nobs)*(*nobs));
	double *ispred_iter = R_VectorInit(*nobs, 0.0);
	int predclass_iter[*nobs];





	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector((*nobs)*(*nobs));
	double *scr2 = R_Vector((*nobs)*(*nobs));



	// stuff that I need to update Si (cluster labels);
	int bcount, btotal;
	double auxm, auxs2,  xb;
	double mudraw, s2draw, sdraw, maxph, denph, cprobh, uu;
	double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw,lgcatt;

	double *ph = R_VectorInit(*nobs, 0.0);
	double *probh = R_VectorInit(*nobs, 0.0);


	//stuff I need for muh
	double sumy_xb,mstar,s2star;

	//suff I need to upate mu0
	double summu;

	//stuff I need to update sig20
	double llo, lln, llr,  os0,ns0;

	//stuff I need to update sig2
	double osig,nsig;

	//stuff I need to update beta
	int bbcount;
	double ld;
	double *sumyxbt = R_VectorInit(*nobs,0.0);
	double *sumXXp = R_VectorInit((*nobs)*(*nobs), 0.0);
	double *Sstar = R_VectorInit((*nobs)*(*nobs), 0.0);
	double *Mstar = R_VectorInit((*nobs)*(*nobs), 0.0);

	// Stuff for out of sample predictions
	double lgcon0, lgcat0, mupred, sig2pred;

	// Stuff to compute lpml, likelihood, and WAIC
	double lpml_iter, elppdWAIC;
	double *CPO = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);



//	Hyper-prior parameters
	// priors for beta
	double mb=modelPriors[0]; double s2b = modelPriors[1];

	// priors for mu0
	double m = modelPriors[2]; double s2 = modelPriors[3];


	// prior values for sig2 and sig02
	double smin=modelPriors[4], smax=modelPriors[5];
	double smin0=modelPriors[6], smax0=modelPriors[7];



	// DP weight parameter
	double Mdp = *M;

	// Similarity function parameters
	// dirichlet denominator parameter
	double *dirweights = R_VectorInit(max_C, simParms[5]);
//	double m0=0.0, s20=0.5, v=1.0, k0=1.0, nu0=1.0;
	double m0=simParms[0];
	double s20=simParms[1];
	double v2=simParms[2];
	double k0=simParms[3];
	double nu0=simParms[4];

	if(*consim==2){
		s20=simParms[2];
		v2=simParms[1];
	}


	Rprintf("M = %f\n", Mdp);
	RprintVecAsMat("dirweights", dirweights, 1,max_C);


	// M-H tuning parameters
	double csigSIG0=mh[0], csigSIG=mh[1];

	ii = 0;

	GetRNGstate();


	// ===================================================================================
	//
	// start of the mcmc algorithm;
	//
	// ===================================================================================

	for(i = 0; i < *draws; i++){

	  if(*verbose){
  		if((i+1) % 1000 == 0){
	  		time_t now;
		  	time(&now);

  			Rprintf("mcmc iter = %d =========================================== \n", i+1);
	  		Rprintf("%s", ctime(&now));
//	  		RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//		  	RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
		  }
	  }

		//////////////////////////////////////////////////////////////////////////////////
		//
		// update the cluster labels using the polya urn scheme of
		// algorithm 8 found in  Radford Neal's
		// "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
		// paper.
		//
		// To accomodate missing here, if a covariate is missing for individual i, then
		// I am simply not evaluating the similarity function for that individual.
		// This is done by dragging around a logical vector that indicates missing or not.
		//
		//////////////////////////////////////////////////////////////////////////////////

//		RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//		Rprintf("nclus_iter = %d\n", nclus_iter);
		for(j = 0; j < *nobs; j++){


//			Rprintf("j = %d =================== \n", j);

//			RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//			Rprintf("nclus_iter = %d\n", nclus_iter);

//			RprintVecAsMat("muh", muh, 1, *nobs);

//			RprintIVecAsMat("nh", nh, 1, *nobs);

//			Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
//			Rprintf("nh[Si_iter[j]-1] = %d\n", nh[Si_iter[j]-1]);

			if(nh[Si_iter[j]-1] > 1){

				// Observation belongs to a non-singleton ...
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] - 1;

			}else{

				// Observation is a member of a singleton cluster ...

				iaux = Si_iter[j];

				if(iaux < nclus_iter){

					// Need to relabel clusters.  I will do this by swapping cluster labels
					// Si_iter[j] and nclus_iter along with cluster specific parameters;


					// All members of last cluster will be assigned subject i's cluster label
					for(jj = 0; jj < *nobs; jj++){

						if(Si_iter[jj] == nclus_iter){

							Si_iter[jj] = iaux;

						}

					}


					Si_iter[j] = nclus_iter;

					// The following steps swaps order of cluster specific parameters
					// so that the newly labeled subjects from previous step retain
					// their correct cluster specific parameters
//					auxt = tau2h[iaux-1];
//					tau2h[iaux-1] = tau2h[nclus_iter-1];
//					tau2h[nclus_iter-1] = auxt;

					auxm = muh[iaux-1];
					muh[iaux-1] = muh[nclus_iter-1];
					muh[nclus_iter-1] = auxm;

					auxs2 = sig2h[iaux-1];
					sig2h[iaux-1] = sig2h[nclus_iter-1];
					sig2h[nclus_iter-1] = auxs2;

//					for(b = 0; b < nb; b++){

//						auxth = thetah[b*(*nobs) + iaux-1];
//						thetah[b*(*nobs) + iaux-1] = thetah[b*(*nobs) + nclus_iter-1];
//						thetah[b*(*nobs) + nclus_iter-1] = auxth;
//					}

					// the number of members in cluster is also swapped with the last
					nh[iaux-1] = nh[nclus_iter-1];
					nh[nclus_iter-1] = 1;

				}


				// Now remove the ith obs and last cluster;
				nh[nclus_iter-1] = nh[nclus_iter-1] - 1;
				nclus_iter = nclus_iter - 1;


			}

//			RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);

//			RprintVecAsMat("muh", muh, 1, *nobs);
//			RprintIVecAsMat("nh", nh, 1, *nobs);

//			Rprintf("nclus_iter = %d\n", nclus_iter);




			// I am summing x beta for different covariate vectors depending on missing
			// Notice that beta_iter contains beta vector for each individual.  This means
			// that there is repetition with-in the vector as individuals with same
			// covariates missing will have the same coefficient vector
//			Rprintf("nb[j] = %d\n", nb[j]);
//			Rprintf("cumindx[j] = %d\n", cumindx[j]);
			xb = 0.0;
			bcount =0;
			for(b = 0; b < ncov; b++){
				if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//					Rprintf("cumindx[j]+bcount = %d\n", cumindx[j]+bcount);
//					xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
					xb = xb + fullXmat[j*ncov+b]*beta_iter[cumindx[j]+bcount];
					bcount = bcount+1;

				}
			}
			btotal = btotal + bcount;
//			Rprintf("xb = %f\n", xb);
//			Rprintf("bcount = %d\n", bcount);

			// The atoms have been relabeled if necessary and now we need to
			// update Si.

//			RprintVecAsMat("btmp", btmp, 1, nb);

			// Begin the cluster probabilities

			for(k=0; k<nclus_iter; k++){

//				Rprintf("k = %d ==================== \n", k);
//				Rprintf("nh[k] = %d\n", nh[k]);


				lgconY = 0.0;
				lgconN = 0.0;
				for(p=0; p<(*ncon); p++){
//					Rprintf("p = %d ====== \n", p) ;
					nhtmp = 0;
					for(jj = 0; jj < *nobs; jj++){
						if(jj != j){
							if(Si_iter[jj] == k+1 & Mcon[jj*(*ncon)+p] == 0){
								xcontmp[nhtmp] = Xcon[jj*(*ncon)+p];
								nhtmp = nhtmp+1;
//								Rprintf("nhtmp = %d\n", nhtmp);
							}
						}
					}

//					Rprintf("nhtmp = %d\n", nhtmp);
//					Rprintf("nh[k] = %d\n", nh[k]);
//					RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp);

					sumx = 0.0;
					sumx2 = 0.0;
					for(t = 0; t < nhtmp; t++){

						sumx = sumx + xcontmp[t];
						sumx2 = sumx2 + xcontmp[t]*xcontmp[t];

					}

//					Rprintf("sumx = %f\n", sumx);
//					Rprintf("sumx2 = %f\n", sumx2);

					if(nhtmp > 0){
						if(*gcontype==1){ // Auxilliary
							if(*consim==1){
								lgconN = lgconN + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							}
							if(*consim==2){
								lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							}

						}
						if(*gcontype==2){ //Double Dipper
							if(*consim==1){
								lgconN = lgconN + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							}
							if(*consim==2){
								lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							}
						}
						if(*gcontype==3){ // Cluster Variance
							lgconN = lgconN + gsimconEV(sumx, sumx2, nhtmp, *alpha,1);
						}

//						Rprintf("lgconN = %f\n", lgconN);
					}

					// now add jth individual back;
//					RprintVecAsMat("xconpred", xconpred, 1, *npred);
//					Rprintf("Mcon[j*(*ncon)+p] = %d\n", Mcon[j*(*ncon)+p]);
					if(Mcon[j*(*ncon)+p] == 0){ // 0 indicates it is not missing
						xcontmp[nhtmp] = Xcon[j*(*ncon)+p];
						sumx = sumx + Xcon[j*(*ncon)+p];
						sumx2 = sumx2 + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
						nhtmp = nhtmp+1;
					}

//					RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp + 1);



//					Rprintf("sumx = %f\n", sumx);
//					Rprintf("sumx2 = %f\n", sumx2);
					if(nhtmp > 0){
						if(*gcontype==1){ // Auxilliary
							if(*consim==1){
								lgconY = lgconY + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							}
							if(*consim==2){
								lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							}
						}
						if(*gcontype==2){ //Double Dipper
							if(*consim==1){
								lgconY = lgconY + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							}
							if(*consim==2){
								lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							}
						}
						if(*gcontype==3){ // variance
							lgconY = lgconY + gsimconEV(sumx, sumx2, nhtmp, *alpha,1);
						}

//						Rprintf("lgcontY = %f\n", lgcont);
//						Rprintf("lgconY = %f\n", lgconY);
//						Rprintf("lgconY = %f\n", lgconY);
					}

				}


//				RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);

				lgcatY=0.0;
				lgcatN=0.0;
				for(p=0; p<(*ncat); p++){
//					Rprintf("p = %d ====== \n", p) ;
					for(c=0;c<Cvec[p];c++){nhc[c]=0;}

					nhtmp = 0;
					for(jj = 0; jj < *nobs; jj++){
//						Rprintf("jj = %d\n", jj);
						if(jj != j){
//							Rprintf("Xcat[jj*(*ncat)+p] = %d\n", Xcat[jj*(*ncat)+p]);
//							Rprintf("Mcat[jj*(*ncat)+p] = %d\n", Mcat[jj*(*ncat)+p]);
//							Rprintf("Si_iter[jj] = %d\n", Si_iter[jj]);

							if(Si_iter[jj]==k+1 & Mcat[jj*(*ncat)+p] == 0){
								nhc[Xcat[jj*(*ncat)+p]] = nhc[Xcat[jj*(*ncat)+p]] + 1; // this needs to be a vector
//                              RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
								xcattmp[nhtmp] = Xcat[jj*(*ncat)+p];
								nhtmp = nhtmp+1;

							}
						}
//						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
					}

//					Rprintf("nhtmp = %d\n", nhtmp);
//					RprintVecAsMat("xcatttmp", xcattmp, 1, nhtmp);
//					Rprintf("sumx = %f\n", sumx);
//					Rprintf("sumx2 = %f\n", sumx2);
//					RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
					if(nhtmp >0){
						if(*gcattype==1){
							lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						}
						if(*gcattype==2){
							lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
						}
						if(*gcattype==3){// Using Entropy instead of variance here
//							lgcatt = gsimconEV(sumx, sumx2, nhtmp, *alpha,1);
							lgcatt = 0.0;
							for(c=0;c<Cvec[p];c++){
                            	if(nhc[c]==0){
                                	lgcatt = lgcatt + 0;
                            	}else{
                                	lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                                                	log((double) nhc[c]/(double) nhtmp)/log(2));
                        		}

							}
						lgcatN = lgcatN + -(*alpha)*lgcatt;
						}
					}

//					Rprintf("Xcat[j*(*ncat)+p] = %f\n", Xcat[j*(*ncat)+p]);
//					Rprintf("Mcat[j*(*ncat)+p] = %d\n", Mcat[j*(*ncat)+p]);

					if(Mcat[j*(*ncat)+p] == 0){
						xcattmp[nhtmp] = Xcat[j*(*ncat)+p];
						nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1;

						nhtmp = nhtmp + 1;
					}

//					RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
//					RprintVecAsMat("xcatttmp", xcattmp, 1, nhtmp+1);
//					Rprintf("nhtmp = %d\n", nhtmp);
//					Rprintf("sumx = %f\n", sumx);
//					Rprintf("sumx2 = %f\n", sumx2);

					if(nhtmp > 0){
						if(*gcattype==1){
							lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						}
						if(*gcattype==2){
							lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
						}
						if(*gcattype==3){// Using Entropy instead of variance here
//							Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, *alpha, 0));
//							lgcatt = gsimconEV(sumx, sumx2, nhtmp, *alpha, 1);
							lgcatt = 0.0;
							for(c=0;c<Cvec[p];c++){
                            	if(nhc[c]==0){
                               	 lgcatt = lgcatt + 0;
                            	}else{
                                	lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                                                	log((double) nhc[c]/(double) nhtmp)/log(2));
                            	}
							}
							lgcatY = lgcatY + -(*alpha)*lgcatt;
						}
					}

				}


//				Rprintf("log(nh[k]) = %f\n", log(nh[k]));
//				RprintVecAsMat("Xcon = ", Xcon, (*ncon), *nobs);
//				Rprintf("Xcon[j] = %f\n", Xcon[j]);



				ph[k] = dnorm(y[j], muh[k] + xb, sqrt(sig2h[k]), 1) +
				        log((double) nh[k]) +
				        lgcatY - lgcatN +
						lgconY - lgconN;

				if(*PPM){
					ph[k] = dnorm(y[j], muh[k] + xb, sqrt(sig2h[k]), 1) +
						    	log((double) nh[k]);  // DP part of cohesion function
				}

				if(*calibrate == 2){

					ph[k] = dnorm(y[j], muh[k] + xb, sqrt(sig2h[k]), 1) +
				            log((double) nh[k]) +
				            (1/((double)*ncon + (double)*ncat))*(lgcatY + lgconY - lgcatN - lgconN);

				}

			}

//			RprintVecAsMat("ph", ph, 1, nclus_iter);

			mudraw = rnorm(mu0_iter, sqrt(sig20_iter));
			sdraw = runif(smin, smax);
//			Rprintf("mudraw = %f\n", mudraw);

			lgcondraw = 0.0;
			for(p=0;p<(*ncon);p++){
				if(Mcon[j*(*ncon)+p] == 0){
					xcontmp[0] = Xcon[j*(*ncon)+p];
//					Rprintf("Xcon[j*(*ncon)+p]=%f\n", Xcon[j*(*ncon)+p]);
					if(*gcontype==1){ // Auxilliary
						if(*consim==1){
							lgcondraw = lgcondraw + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0, 1);
						}
						if(*consim==2){
							lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 0,0, 1);
						}
					}
					if(*gcontype==2){ // Double Dipper
						if(*consim==1){
							lgcondraw = lgcondraw + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p], 1, 1, 0, 1);
						}
						if(*consim==2){
							lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 1, 0, 1);
						}
					}
					if(*gcontype==3){ // Variance
						lgcondraw = lgcondraw + gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,*alpha,1);
					}
//					Rprintf("lgcondraw = %f\n", lgcondraw);
				}
			}

//			Rprintf("lgcondraw = %f\n", lgcondraw);

			lgcatdraw = 0.0;
			for(p=0;p<(*ncat);p++){
				for(c=0;c<Cvec[p];c++){nhc[c] = 0;}
	            if(Mcat[j*(*ncat)+p] == 0){

					nhc[Xcat[j*(*ncat)+p]] = 1;
//					RprintIVecAsMat("nhc =", nhc, 1, Cvec[p]);

//					Rprintf("Xcat[j*(*ncon)+p] = %d\n", Xcat[j*(*ncon)+p]);

					if(*gcattype==1){
						lgcatdraw = lgcatdraw + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
					}
					if(*gcattype==2){
						lgcatdraw = lgcatdraw + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
					}
					if(*gcattype==3){
//						lgcatt = gsimconEV(xcattmp[0], xcattmp[0]*xcattmp[0], 1,*alpha, 1);
//						lgcatdraw = lgcatdraw + lgcatt;
						lgcatdraw = lgcatdraw + -(*alpha)*0;
					}

				}
			}
//			Rprintf("lgcatdraw = %f\n", lgcatdraw);

//			Rprintf("dnorm(y[j],mudraw,sqrt(sig2_iter),1) = %f\n", dnorm(y[j],mudraw,sqrt(sig2_iter),1));
//			Rprintf("log(Mdp) = %f\n", log(Mdp));
//			Rprintf("lgcondraw = %f\n", lgcondraw);
//			Rprintf("lgcatdraw = %f\n", lgcatdraw);


			ph[nclus_iter] = dnorm(y[j],mudraw + xb,sdraw,1) +
			                 log(Mdp) +
			                 lgcondraw +
			                 lgcatdraw;

			if(*PPM){
				ph[nclus_iter] = dnorm(y[j],mudraw + xb,sdraw,1) +
				                 log(Mdp); //DP part
			}

			if(*calibrate==2){
				ph[nclus_iter] = dnorm(y[j],mudraw + xb,sdraw,1) +
				                 log(Mdp) +
				                 (1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);
			}


//			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

			maxph = ph[0];
			for(k = 1; k < nclus_iter+1; k++){
				if(maxph < ph[k]) maxph = ph[k];
			}
//			Rprintf("maxph = %f\n", maxph);


			denph = 0.0;
			for(k = 0; k < nclus_iter+1; k++){

				ph[k] = exp(ph[k] - maxph);
//				ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
				denph = denph + ph[k];

			}

//			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

			for(k = 0; k < nclus_iter+1; k++){

				probh[k] = ph[k]/denph;

			}
//			Rprintf("denph = %f\n", denph);

//			RprintVecAsMat("probh", probh, 1, nclus_iter+1);

			uu = runif(0.0,1.0);
//			Rprintf("uu = %f\n", uu);



			cprobh= 0.0;;
			for(k = 0; k < nclus_iter+1; k++){

				cprobh = cprobh + probh[k];

				if (uu < cprobh){

					iaux = k+1;
					break;
				}
			}


// 			Rprintf("iaux = %d\n \n \n", iaux);

			if(iaux <= nclus_iter){

				Si_iter[j] = iaux;
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] + 1;

			}else{

				nclus_iter = nclus_iter + 1;
				Si_iter[j] = nclus_iter;
				nh[Si_iter[j]-1] = 1;

//				tau2h[Si_iter[j]-1] = tau2draw;
//				lamh[Si_iter[j]-1] = lamdraw;
				muh[Si_iter[j]-1] = mudraw;
				sig2h[Si_iter[j]-1] = sdraw*sdraw;

//				for(b = 0; b < nb; b++){
//					thetah[b*(*nobs) + Si_iter[j]-1] = thetadraw[b];
//				}

//				RprintVecAsMat("thetah after", thetah, nb, *nobs);
			}


// 			Rprintf("Si_iter = %d\n \n \n", Si_iter[j]);

//			RprintVecAsMat("tau2h", tau2h, 1, *nobs);
//			RprintVecAsMat("lamh", lamh, 1, *nobs);
//			RprintVecAsMat("thetah", thetah, nb, *nobs);



		}





//		Rprintf("nclus_iter = %d\n", nclus_iter);
//		RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//		RprintIVecAsMat("nh", nh, 1, nclus_iter);
//		RprintIVecAsMat("nb", nb, 1, *nobs);



		//////////////////////////////////////////////////////////////////////////////////
		//
		// update muh's cluster specific intercepts
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(k = 0; k < nclus_iter; k++){
			sumy_xb = 0.0;
			for(j = 0; j < *nobs; j++){
				if(Si_iter[j] == k+1){
					xb = 0.0;
					bcount =0;
					for(b = 0; b < ncov; b++){
						if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//							Rprintf("cumindx[j]+bcount = %d\n", cumindx[j]+bcount);
//							xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
							xb = xb + fullXmat[j*ncov+b]*beta_iter[cumindx[j]+bcount];
							bcount = bcount+1;
						}
					}
//					Rprintf("xb = %f\n", xb);
					sumy_xb = sumy_xb + (y[j] - xb);
				}

			}
//			Rprintf("sumy_xb = %f\n", sumy_xb);
//			Rprintf("nh[k] = %d\n",  nh[k]);

			s2star = 1/((double) nh[k]/sig2h[k] + 1/sig20_iter);
			mstar = s2star*( (1/sig2h[k])*sumy_xb + (1/sig20_iter)*mu0_iter);

//			Rprintf("sig2_iter = %f\n", sig2_iter);
//			Rprintf("mstar = %f\n", mstar);
//			Rprintf("sqrt(s2star) = %f\n", sqrt(s2star));

			muh[k] = rnorm(mstar, sqrt(s2star));
//			muh[k] = 10.0;
//			Rprintf("muh = %f\n", muh[k]);
		}

//		RprintVecAsMat("muh", muh, 1, nclus_iter);


		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update mu0  prior mean of muh
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		summu = 0.0;
		for(k = 0; k < nclus_iter; k++){

			summu = summu + muh[k];

		}

//		Rprintf("summu = %f\n", summu);

		s2star = 1/(((double) nclus_iter/sig20_iter) + (1/s2));
		mstar = s2star*((1/sig20_iter)*summu + (1/s2)*m);

		mu0_iter = rnorm(mstar, sqrt(s2star));

//		mu0_iter = 0.0;




//		Rprintf("mu0_iter = %f\n", mu0_iter);

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig20  prior variance of muh
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		os0 = sqrt(sig20_iter);
		ns0 = rnorm(os0,csigSIG0);

//		Rprintf("os0 = %f\n", os0);
//		Rprintf("ns0 = %f\n", ns0);

		if(ns0 > 0){

			lln = 0.0;
			llo = 0.0;
			for(k = 0; k < nclus_iter; k++){

				llo = llo + dnorm(muh[k], mu0_iter, os0,1);
				lln = lln + dnorm(muh[k], mu0_iter, ns0,1);
			}

			llo = llo + dunif(os0, smin0, smax0, 1);
			lln = lln + dunif(ns0, smin0, smax0, 1);


//			Rprintf("llo = %f\n", llo);
//			Rprintf("lln = %f\n", lln);

			llr = lln - llo;
			uu = runif(0,1);

//			Rprintf("llr = %f\n", llr);
//			Rprintf("log(uu) = %f\n", log(uu));
			if(log(uu) < llr){
				sig20_iter = ns0*ns0;
			}
		}

//		Rprintf("sig0_iter = %f\n", sqrt(sig20_iter));
//		Rprintf("sig20_iter = %f\n", sig20_iter);

//		s2b0_iter = 10*10;



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig2h  cluster specific variance parameters with metropolis and
		// Uniform prior.
		//
		//////////////////////////////////////////////////////////////////////////////////
//		RprintVecAsMat("y", y, 1, *nobs);
//		RprintVecAsMat("muh", muh, 1, nclus_iter);
		for(k = 0; k < nclus_iter; k++){
			osig = sqrt(sig2h[k]);
			nsig = rnorm(osig,csigSIG);

			if(nsig > 0 & nsig < smax){

				lln = 0.0;
				llo = 0.0;
				for(j = 0; j < *nobs; j++){
					if(Si_iter[j] == k+1){
						xb = 0.0;
						bcount =0;
						for(b = 0; b < ncov; b++){
							if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//								xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
								xb = xb + fullXmat[j*ncov+b]*beta_iter[cumindx[j]+bcount];
								bcount = bcount+1;
							}
						}

						llo = llo + dnorm(y[j], muh[k] + xb, osig,1);
						lln = lln + dnorm(y[j], muh[k] + xb, nsig,1);
					}
				}

//				Rprintf("ms = %f\n", ms);
//				Rprintf("osig = %f\n", osig);
//				Rprintf("nsig = %f\n", nsig);
				llo = llo + dunif(osig, smin, smax, 1);
				lln = lln + dunif(nsig, smin, smax, 1);


//				Rprintf("llo = %f\n", llo);
//				Rprintf("lln = %f\n", lln);

				llr = lln - llo;
				uu = runif(0,1);

//				Rprintf("llr = %f\n", llr);
//				Rprintf("log(uu) = %f\n", log(uu));

				if(log(uu) < llr){
					sig2h[k] = nsig*nsig;
				}
			}

//			sig2h[k] = 1.0;

		}


//		RprintVecAsMat("sig2h", sig2h, 1, nclus_iter);




		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update beta for each cluster configuration.  There are g of these and then
		// they will need to be saved in beta_iter.
		//
		//////////////////////////////////////////////////////////////////////////////////
		btotal = 0;
		for(g = 0; g < ng; g++){

//			Rprintf("g = %d ==== \n", g);

			for(b = 0; b < ncov; b++){
				sumyxbt[b] = 0.0;
				for(bb = 0; bb < ncov; bb++){
					sumXXp[b*ncov+bb] = 0.0;
				}
			}


			for(j = 0; j < *nobs; j++){
//				Rprintf("j = %d === \n", j);
				if(gvec[j] == g+1){
//					Rprintf("j = %d === \n", j);

//					Rprintf("sig2h[Si_iter[j]-1] = %f\n", sig2h[Si_iter[j]-1]);
					bcount =0;
					for(b = 0; b < ncov; b++){

						if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//							Rprintf("nb[j] = %d\n", nb[j]);
							sumyxbt[bcount] = sumyxbt[bcount] +
											(1/sig2h[Si_iter[j]-1])*fullXmat[j*ncov+b]*
								                (y[j] - muh[Si_iter[j]-1]);

						    bbcount = 0;
							for(bb = 0; bb < ncov; bb++){
								if(fullM[j*ncov+bb] == 0){ // fullM == 0 implies NOT missing.
//									Rprintf("bb = %d ==== \n", bb);
//									Rprintf("bbcount = %d\n", bbcount);
									sumXXp[bcount*nb[j]+bbcount] = sumXXp[bcount*nb[j]+bbcount] +
									                    (1/sig2h[Si_iter[j]-1])*
									                    fullXmat[j*ncov+b]*fullXmat[j*ncov+bb];

									bbcount = bbcount+1;
								}
							}

//							Rprintf("bbcount = %d\n", bbcount);
							bcount = bcount + 1;

						}

//						Rprintf("bcount = %d\n", bcount);

					}

				}
//				RprintVecAsMat("sumXXp", sumXXp, bcount, bcount);

			}
//			RprintVecAsMat("sumXXp", sumXXp, bcount, bcount);
//			RprintVecAsMat("sumyxbt", sumyxbt, 1, bcount);



			for(b = 0; b < bcount; b++){
				sumyxbt[b] = sumyxbt[b] + 1/s2b*mb;
				// If this is commented out it implies that I am using a p(beta) propto 1 prior
				// and updating mubeta and s2beta does nothing as they are not included in the
				// prior.
				for(bb = 0; bb < bcount; bb++){
					Sstar[b*bcount + bb] = sumXXp[b*bcount+bb];

				// If this is commented out it implies that I am using a p(beta) propto 1 prior
					if(b==bb) Sstar[b*bcount + bb] = sumXXp[b*bcount+bb] + 1/s2b;
//				Sstar[b*ncov + bb] = sumXXp[b*ncov+bb];

				}
			}


			cholesky(Sstar, bcount, &ld);
			inverse_from_cholesky(Sstar, scr1, scr2, bcount); //Sstar is now an inverse;
//			RprintVecAsMat("Sstar", Sstar, bcount, bcount);


			matrix_product(Sstar, sumyxbt, Mstar, bcount, 1, bcount);

//			RprintVecAsMat("Mstar", Mstar, 1, bcount);


			cholesky(Sstar, bcount , &ld);

			ran_mvnorm(Mstar, Sstar, bcount, scr1, scr2);

//			RprintVecAsMat("beta_iter", scr2, 1, bcount);


			for(j = 0; j < *nobs; j++){
				if(gvec[j] == g+1){
					for(b = 0; b < bcount; b++){
						beta_iter[cumindx[j]+b]=scr2[b];
					}
				}
			}

		}

//		RprintIVecAsMat("gvec",gvec,1,*nobs);
//		RprintVecAsMat("beta_iter", beta_iter, 1, (*nobs)*(ncov));


		//////////////////////////////////////////////////////////////////////////////////
		//
		// evaluating likelihood that will be used to calculate LPML?
		// (see page 81 Christensen Hansen and Johnson)
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(i > (*burn-1) & i % (*thin) == 0){
			for(j = 0; j < *nobs; j++){
//				Rprintf("j = %d\n", j);

				mudraw = muh[Si_iter[j]-1];
				s2draw = sig2h[Si_iter[j]-1];
//				Rprintf("mudraw = %f\n", mudraw);
				xb = 0.0;
				bcount =0;
				for(b = 0; b < ncov; b++){
					if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//						xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
						xb = xb + fullXmat[j*ncov+b]*beta_iter[cumindx[j]+bcount];
						bcount = bcount+1;
					}
				}

				like_iter[j] = dnorm(y[j], mudraw + xb, sqrt(s2draw), 0);
//				Rprintf("like_iter = %f\n", like_iter[j]);
//				Rprintf("like_iter = %40.9f\n", like_iter[j]);

				// These are needed for WAIC
				mnlike[j] = mnlike[j] + (like_iter[j])/(double) nout;
				mnllike[j] = mnllike[j] + log(like_iter[j])/(double) nout;

				CPO[j] = CPO[j] + (1/(double) nout)*
				                  (1/like_iter[j]);

//				Rprintf("CPO = %f\n", CPO[j]);



			}

//			if(i == (*draws-1)) Rprintf("xb = %f\n", xb);

		}


		lpml_iter=0.0;
		if(i == (*draws-1)){

			for(jj = 0; jj < *nobs; jj++){

				lpml_iter = lpml_iter + log(1/CPO[jj]);

			}

			lpml[0] = lpml_iter;

		}

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Computing WAIC  (see Gelman article in lit review folder)
		// An issue with WAIC is that considering spatially structure data
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		elppdWAIC = 0.0;
		if(i == (*draws - 1)){
//			RprintVecAsMat("mnlike",mnlike, 1,*nobs);
//			RprintVecAsMat("mnllike", mnllike, 1, *nobs);

			for(j = 0; j < *nobs; j++){

				elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
			}

			WAIC[0] = -2*elppdWAIC;

//			Rprintf("WAIC = %f\n", -2*elppdWAIC);
		}



		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// out of sample prediction using posterior predictive?
		//
		////////////////////////////////////////////////////////////////////////////////////////////
//		Rprintf("npred = %d\n", *npred);
//		Rprintf("nclus_iter = %d\n", nclus_iter);
		if(i > (*burn-1) & i % (*thin) == 0){
			for(pp = 0; pp < *npred; pp++){

//				Rprintf("pp = %d ===================================================== \n", pp+1);
//				Rprintf("nclus_iter = %d\n", nclus_iter);

				for(k = 0; k < nclus_iter; k++){
//					Rprintf("k = %d  ========== \n", k);


					lgconN=0.0, lgconY=0.0;
					for(p=0; p<(*ncon); p++){
//						Rprintf("p = %d ====== \n", p) ;
						nhtmp = 0;
						sumx = 0.0;
						sumx2 = 0.0;
						for(j = 0; j < *nobs; j++){
							if(Si_iter[j] == k+1 & Mcon[j*(*ncon)+p] == 0){ // Mcon == 0 implies not missing
								xcontmp[nhtmp] = Xcon[j*(*ncon)+p]; //create cluster specific x-vector
								nhtmp = nhtmp+1;
								sumx = sumx + Xcon[j*(*ncon)+p];
								sumx2 = sumx2 + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
//								Rprintf("nhtmp = %d\n", nhtmp);
							}
						}

//						Rprintf("nhtmp = %d\n", nhtmp);
//						Rprintf("nh[k] = %d\n", nh[k]);
//						RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp);

//						Rprintf("sumx = %f\n", sumx);
//						Rprintf("sumx2 = %f\n", sumx2);




						if(nhtmp > 0){
							if(*gcontype==1){// Auxilliary
								if(*consim==1){// NN
									lgconN = lgconN + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								}
							}
							if(*gcontype==2){// Double Dipper
								if(*consim==1){// NN
									lgconN = lgconN + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								}
							}
							if(*gcontype==3){
								lgconN = lgconN + gsimconEV(sumx, sumx2, nhtmp,*alpha, 1);
							}
						}
//						Rprintf("lgconN = %f\n", lgconN);

					// now add ppth prediction to cluster;
//						Rprintf("xconpred[pp] = %f\n", xconpred[pp]);
						if(Mconp[pp*(*ncon)+p] == 0){
							xcontmp[nhtmp] = Xconp[pp*(*ncon)+p];
							sumx = sumx + Xconp[pp*(*ncon)+p];
							sumx2 = sumx2 + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
							nhtmp = nhtmp + 1;
						}

//						RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp + 1);

//						Rprintf("sumx = %f\n", sumx);
//						Rprintf("sumx2 = %f\n", sumx2);
						if(nhtmp > 0){
							if(*gcontype==1){ // Auxilliary
								if(*consim==1){// NN
									lgconY = lgconY + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								}
							}
							if(*gcontype==2){ // Double Dipper
								if(*consim==1){// NN
									lgconY = lgconY + gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								}
							}
							if(*gcontype==3){ // Variance
								lgconY = lgconY + gsimconEV(sumx, sumx2, nhtmp,*alpha, 1);
							}

//							Rprintf("lgconY = %f\n", lgconY);
						}

					}

//					Rprintf("lgconN = %f\n", lgconN);
//					Rprintf("lgconY = %f\n", lgconY);

					lgcatY=0.0, lgcatN=0.0;
					for(p=0; p<(*ncat); p++){
//						Rprintf("p = %d ====== \n", p) ;
						for(c=0;c<Cvec[p];c++){nhc[c]=0;}
						nhtmp=0;
						for(j = 0; j < *nobs; j++){
//							Rprintf("j = %d\n", j);
//							Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
//							Rprintf("Xcat[j*(*ncat)+p] = %d\n", Xcat[j*(*ncat)+p]);

							if(Si_iter[j]==k+1 & Mcat[j*(*ncat)+p] == 0){
								nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1; // this needs to be a vector
								xcattmp[nhtmp] = Xcat[j*(*ncat)+p];
								nhtmp = nhtmp+1;

//                              	RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
							}
						}


//						Rprintf("nhtmp =%d\n", nhtmp);
//						RprintIVecAsMat("xcatttmp", xcattmp, 1, nhtmp);

//						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
						if(nhtmp > 0){
							if(*gcattype==1){// Auxilliary
								lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							}
							if(*gcattype==2){// Double Dipper
								lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
							}
							if(*gcattype==3){
//								lgcatt = gsimconEV(sumx, sumx2, nhtmp,*alpha, 1);
								lgcatt = 0.0;
								for(c=0;c<Cvec[p];c++){
                            		if(nhc[c]==0){
                                		lgcatt = lgcatt + 0;
                            		}else{
                                		lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                                    	            	log((double) nhc[c]/(double) nhtmp)/log(2));
                        			}
								}
								lgcatN = lgcatN + -(*alpha)*lgcatt;
							}
//							Rprintf("lgcatN = %f\n", lgcatN);
						}



//						Rprintf("Xcat[j*(*ncat)+p] = %d\n", Xcat[j*(*ncat)+p]);
						if(Mcatp[pp*(*ncat)+p] == 0){
								nhc[Xcatp[pp*(*ncat)+p]] = nhc[Xcatp[pp*(*ncat)+p]] + 1;
								nhtmp=nhtmp + 1;
						}
//						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);

						if(nhtmp > 0){
							if(*gcattype==1){// Auxilliary
								lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							}
							if(*gcattype==2){// Double Dipper
								lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
							}
							if(*gcattype==3){// Use entropy
//								lgcatt = gsimconEV(sumx, sumx2, nhtmp, *alpha,  1);
								lgcatt = 0.0;
								for(c=0;c<Cvec[p];c++){
                            		if(nhc[c]==0){
                                		lgcatt = lgcatt + 0;
                            		}else{
                                		lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                                    	            	log((double) nhc[c]/(double) nhtmp)/log(2));
                        			}
								}
								lgcatY = lgcatY + -(*alpha)*lgcatt;
							}
						}
					}

//					Rprintf("lgcatN = %f\n", lgcatN);
//					Rprintf("lgcatY = %f\n", lgcatY);

					ph[k] = log((double) nh[k]) +
				         	lgcatY - lgcatN +
							lgconY - lgconN;


					if(*PPM) ph[k] = log((double) nh[k]);


					if(*calibrate == 2){

						ph[k] = log((double) nh[k]) +
						           (1/((double)*ncon + (double)*ncat))*
						            (lgcatY + lgconY - lgcatN - lgconN);

//						Rprintf("ph[k] = %f\n", ph[k]);

					}

				}

//				Rprintf("nclus_iter = %d\n", nclus_iter);

//				RprintVecAsMat("ph", ph, 1, nclus_iter);

				lgcon0=0.0;
				for(p=0;p<*ncon;p++){ // 0 means that data are not missing

//					Rprintf("p = %d\n", p);
//					Rprintf("Mconp[pp*(*ncon)+p] = %d\n", Mconp[pp*(*ncon)+p]);

					if(Mconp[pp*(*ncon)+p] == 0){
						xcontmp[0] = Xconp[pp*(*ncon)+p];
//						Rprintf("xcontmp = %f\n", xcontmp[0]);
//						Rprintf("*gcontype = %d\n", *gcontype);
//						Rprintf("*consim = %d\n", *consim);

						if(*gcontype==1){ // Auxilliary
							if(*consim==1){
								lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0,1);
							}
							if(*consim==2){
								lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,0,1);
							}
						}
						if(*gcontype==2){ // Double Dipper
							if(*consim==1){
								lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,0,1);
							}
							if(*consim==2){
								lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 0,1);
							}
						}
						if(*gcontype==3){ // Variance
							lgcon0 = lgcon0 + gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,*alpha,1);
						}

					}
				}





				lgcat0 = 0.0;
				for(p=0;p<(*ncat);p++){
//					Rprintf("p = %d\n", p);
					for(c=0;c<Cvec[p];c++) nhc[c] = 0;

//					Rprintf("Mcatp[pp*(*ncat)+p] = %d\n", Mcatp[pp*(*ncat)+p]);
//					Rprintf("Xcatp[pp*(*ncat)+p] = %d\n", Xcatp[pp*(*ncat)+p]);
	        	    if(Mcatp[pp*(*ncat)+p] == 0){

						nhc[Xcatp[pp*(*ncat)+p]] = 1;
//						RprintIVecAsMat("nhc =", nhc, 1, Cvec[p]);


						if(*gcattype==1){
							lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						}
						if(*gcattype==2){
							lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
						}
						if(*gcattype==3){
//							lgcatt = gsimconEV(xcattmp[0], xcattmp[0]*xcattmp[0], 1,*alpha, 1);
//							lgcatdraw = lgcatdraw + lgcatt;
							lgcat0 = lgcat0 +  -(*alpha)*0;
						}

					}

				}


				ph[nclus_iter] = log((double) Mdp) + lgcon0 + lgcat0;

				if(*calibrate==2){
					ph[nclus_iter] = log(Mdp) +
									(1/((double)*ncon + (double)*ncat))*(lgcon0 + lgcat0);
				}

				if(*PPM) ph[nclus_iter] = log(Mdp);

//				RprintVecAsMat("ph", ph, 1, nclus_iter+1);

				maxph = ph[0];
				for(k = 1; k < nclus_iter+1; k++){
					if(ph[k] > maxph) maxph=ph[k];
				}
//				Rprintf("maxph = %f\n", maxph);

				denph = 0.0;
				for(k = 0; k < nclus_iter+1; k++){

					ph[k] = exp(ph[k] - maxph);
//					ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
					denph = denph + ph[k];

				}

//				RprintVecAsMat("ph", ph, 1, nclus_iter+1);


				for(k = 0; k < nclus_iter+1; k++){

					probh[k] = ph[k]/denph;

				}
//				Rprintf("denph = %f\n", denph);

//				RprintVecAsMat("probh", probh, 1, nclus_iter+1);

				uu = runif(0.0,1.0);
//				Rprintf("uu = %f\n", uu);

				cprobh= 0.0;

				for(k = 0; k < nclus_iter+1; k++){

					cprobh = cprobh + probh[k];

					if (uu < cprobh){

						iaux = k+1;
						break;
					}

				}

//				Rprintf("iaux = %d\n", iaux);
//				RprintVecAsMat("muh", muh, 1, nclus_iter);

				if(iaux <= nclus_iter){
					mupred = muh[(iaux-1)];
					sig2pred = sig2h[(iaux-1)];
				}else{
					mupred = rnorm(mu0_iter,sqrt(sig20_iter));
					sig2pred = runif(smin, smax);
					sig2pred = sig2pred*sig2pred;
				}

				xb = 0.0;
				for(j = 0; j < *nobs; j++){
					if(gvec[j] == gvecp[pp]){
						bcount =0;
						xb = 0.0;
						for(b = 0; b < ncov; b++){
							if(fullMp[pp*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//								Rprintf("cumindx[j]+bcount = %d\n", cumindx[j]+bcount);
//								xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
//								Rprintf("fullXmatp[pp*ncov+b] = %f\n", fullXmatp[pp*ncov+b]);
//								Rprintf("beta_iter[cumindx[j]+bcount] = %f\n", beta_iter[cumindx[j]+bcount]);
								xb = xb + fullXmatp[pp*ncov+b]*beta_iter[cumindx[j]+bcount];
								bcount = bcount+1;

							}
						}
						break;
					}

				}
//				xb = 0;
//				Rprintf("j = %d\n", j);
//				Rprintf("mupred = %f\n", mupred);
//				Rprintf("xb = %f\n", xb);
//				Rprintf("sqrt(sig2pred) = %f\n", sqrt(sig2pred));

				ppred_iter[pp] = rnorm(mupred + xb, sqrt(sig2pred));
				predclass_iter[pp] = iaux;

//				Rprintf("ppred = %f\n", ppred_iter[pp]);
				mupred = 0.0;
				for(k = 0; k < nclus_iter; k++){
					mupred = mupred + probh[k]*muh[k];
				}
				mupred = mupred + probh[nclus_iter]*rnorm(mu0_iter,sqrt(sig20_iter));

				pred2_iter[pp] = mupred + xb;

			}

		}


		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// in sample prediction to assess model fit
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		if(i > (*burn-1) & i % (*thin) == 0){

			for(j = 0; j < *nobs; j++){

				xb = 0.0;
				bcount =0;
				for(b = 0; b < ncov; b++){
					if(fullM[j*ncov+b] == 0){ // fullM == 0 implies NOT missing.
//						Rprintf("cumindx[j]+bcount = %d\n", cumindx[j]+bcount);
//						xb = xb + fullXmat[j*ncov+b]*beta_iter[j*nb[j]+bcount];
						xb = xb + fullXmat[j*ncov+b]*beta_iter[cumindx[j]+bcount];
						bcount = bcount+1;

					}
				}

				ispred_iter[j] = rnorm(muh[Si_iter[j]-1] + xb, sqrt(sig2h[Si_iter[j]-1]));


			}


//			RprintVecAsMat("ispred_iter", ispred_iter, 1, *nobs);

		}



		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// Save MCMC iterates															//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){

			mu0[ii] = mu0_iter;
			sig20[ii] = sig20_iter;

			nclus[ii] = nclus_iter;

			for(j = 0; j < *nobs; j ++){

				mu[ii*(*nobs) + j] = muh[Si_iter[j]-1];
				sig2[ii*(*nobs) + j] = sig2h[Si_iter[j]-1];
				Si[ii*(*nobs) + j] = Si_iter[j];

				like[ii*(*nobs) + j] = like_iter[j];
				ispred[ii*(*nobs) + j] = ispred_iter[j];

			}

			for(b = 0; b < nbtotal; b++){
				beta[ii*nbtotal + b] = beta_iter[b];
			}
			for(pp = 0; pp < *npred; pp++){

				ppred[ii*(*npred) + pp] = ppred_iter[pp];
				predclass[ii*(*npred) + pp] = predclass_iter[pp];
				pred2[ii*(*npred) + pp] = pred2_iter[pp];

			}


			ii = ii+1;

		}



	}




	PutRNGstate();




}
