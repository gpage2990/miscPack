# This is an MCMC wrapper that fits missing PPMx model with the possibility of
# employing the pattern submodel technique in the likelihood if so desired.
# y<-dat[,2];Xmat <- dat[,c(3:11,13)]; Xpred=NULL;psm=FALSE;M=1;gcontype=1;gcattype=1;
# consim=1;calibrate=0;alpha=1;PPM=FALSE;draws=11000;burn=1000;thin=10
# simparms=c(0.0, 1.0, 1.0, 1.0, 1.0, 0.1);
# modelPriors=c(0,10*10,0,10*10, 0,10.0,0,20);
ppmx.missing  <- function(y, X, Xpred=NULL, psm=FALSE,PPM=FALSE,
                    M=1,
                    gcontype=1,gcattype=1,
                    consim=1, calibrate=0,alpha=1,
                    simParms=c(0.0, 1.0, 1.0, 1.0, 1.0, 0.1),
                    modelPriors=c(0.0,10.0^2, 0.0, 10.0^2, 0.0, 10.0, 0.0, 10.0),
                    mh=c(0.5, 0.5),
                    draws=11000,burn=1000,thin=10, verbose=FALSE){

# consim: argument indicating which similarity function to use for continuous covariate
#         1 - NN
#         2 - NNIG
# calibrate: argument indicatiing if similiarity should be calibrated
#         1 - No included
#         2 - coarsening (raising similarity value to 1/p power)
# alpha: parameter associated with variance similarity g(x_j) = alpha exp(-var(x_j));
#
# simparms: the order here is c(mu0, s20,   v,   k0, nu0, a0)
#
# modelPriors: vector containing (mub, s2b, m0,  s20	  sig2min, sig2max, sig20min, sig20max)
#

  out <- NULL

  Xmat <- X

  if(!is.data.frame(Xmat)) Xmat <- data.frame(Xmat)
  if(!is.data.frame(Xpred) & !is.null(Xpred)) {Xpred <- data.frame(Xpred); colnames(Xpred) <- colnames(Xmat)}

  ## Function that allows me to standardize each column of data matrix ignoring missing values that are set
  ## to 999
  stdM <- function(x){
	  xmi <- x[!is.na(x)]
	  mnx <- mean(xmi)
	  sdx <- sd(xmi)
	  stdx <- (xmi - mnx)/sdx
	  x[!is.na(x)] <- stdx
	  x
  }
  # Function that relabels categorical variables to begin with 0
  relab <- function(x) as.numeric(as.factor(as.character(x))) - 1

  nout <- (draws-burn)/thin

  nobs <- length(y)

  npred <- nrow(Xpred)

  Xall <- rbind(Xmat, Xpred)

  Mmat <- matrix(0,nrow=nrow(Xall), ncol=ncol(Xall))
  Mmat[is.na(Xall)] <- 1
  Xall[is.na(Xall)] <- 999

  classes <- sapply(Xall, class)
  catvars <- classes %in% c("factor","character")
  if(sum(!catvars) > 0){
  	Xconstd <- apply(Xall[,!catvars, drop=FALSE], 2, stdM)
	  Xcon <- Xconstd[1:nobs,,drop=FALSE];
	  Mcon <- Mmat[1:nobs,!catvars,drop=FALSE]
    ncon <- ncol(Xcon)

  }else{
    Xcon <- cbind(rep(0,nobs));
    Mcon <- cbind(rep(0,nobs))
    ncon <- 0
  }

  if(sum(catvars) > 0){
    # Change the factors or characters into integers with category starting at 0.
    Xcatall <- sapply(Xall[, catvars,drop=FALSE], relab)
    Cvec <- apply(Xcatall,2,function(x)length(unique(x)))

	Xcat <- Xcatall[1:nobs,,drop=FALSE];
	Mcat <- Mmat[1:nobs,catvars,drop=FALSE]
	ncat <- ncol(Xcat)

	}else{
	  Xcat <- cbind(rep(0,nobs));
	  Mcat <- cbind(rep(0,nobs))
	  Cvec <- 0
	  ncat <- 0
	}


	cat("ncon = ", ncon, "\n")
	cat("ncat = ", ncat, "\n")



	if(!is.null(npred)){
	  if(sum(!catvars) > 0){
      Xconp <- Xconstd[-c(1:nobs),,drop=FALSE];
      Mconp <- Mmat[-c(1:nobs),!catvars,drop=FALSE]
	  }else{
		  Xconp <- cbind(rep(0,npred)); Mconp <- cbind(rep(0,npred))
	  }

	  if(sum(catvars) > 0){
	    Xcatp <- Xcatall[-c(1:nobs),,drop=FALSE];
	    Mcatp <- Mmat[-c(1:nobs),catvars,drop=FALSE]
	  }else{
		  Xcatp <- cbind(rep(0,npred)); Mcatp <- cbind(rep(0,npred))
	  }

 	} else {
 		npred <- 0
		Xconp <- cbind(rep(0,1)); Mconp <- cbind(rep(0,1))
		Xcatp <- cbind(rep(0,1)); Mcatp <- cbind(rep(0,1))
	}




	if(psm){
		mcmc.out <- NULL
		## Function that identifies individuals covariate configuration group
		## Missing values are identified with 999.
		covgrp <- function(Xcon=NULL, Xcat=NULL){

			out <- NULL
			X <- cbind(Xcon, Xcat)
			ncov <- dim(X)[2]
			Xlg <- X == 999
			unXlg <- unique(Xlg)
			g <- numeric()
			for(j in 1:dim(unXlg)[1]){
				g[apply(unXlg[j,] == t(Xlg),2,sum) == ncov] <- j
			}

			out$g <- g
			out$beta_dim <- apply(1-unXlg,1,sum)
			out
		}


		gvec_all <- covgrp(Xall[!catvars], Xall[catvars])$g
		gvec <- gvec_all[1:nobs]
		gvecp <- gvec_all[-(1:nobs)]

		nb_all <- apply(cbind(Xall) != 999, 1, sum)
		nb <- nb_all[1:nobs]
		nbp <- nb_all[-(1:nobs)]

		beta <- matrix(0, nrow=nout, ncol=sum(nb))
		mu <- sig2 <- Si <- like <- ispred <- matrix(1,nrow=nout,ncol=nobs)
		mu0 <- sig20 <- nclus <- rep(1,nout)
		ppred <- predclass <- pred2 <- matrix(1, nrow=nout, ncol=npred)
		WAIC <- lpml <- rep(1,1)


		C.out <- .C("mcmc_missing_psm",
					as.integer(draws), as.integer(burn), as.integer(thin),
					as.integer(nobs),as.integer(ncon), as.integer(ncat),
					as.integer(Cvec), as.integer(PPM),
					as.integer(gvec), as.integer(nb),
					as.integer(gvecp), as.integer(nbp),
					as.double(M),as.integer(gcontype), as.integer(gcattype),
					as.double(y),
					as.double(t(Xcon)),as.integer(t(Xcat)),
					as.integer(t(Mcon)), as.integer(t(Mcat)),
					as.integer(npred),
					as.double(t(Xconp)),as.integer(t(Xcatp)),
					as.integer(t(Mconp)), as.integer(t(Mcatp)),
					as.integer(calibrate),as.integer(consim),as.double(alpha),
					as.double(simParms), as.double(modelPriors),
					as.double(mh),
					mu.out=as.double(mu), beta.out=as.double(beta), sig2.out=as.double(sig2),
					mu0.out=as.double(mu0), sig20.out=as.double(sig20),
					Si.out=as.integer(Si), nclus.out=as.integer(nclus),
					like.out=as.double(like), WAIC.out=as.double(WAIC),
					lpml.out=as.double(lpml),ispred.out=as.double(ispred),
					ppred.out=as.double(ppred),predclass.out=as.integer(predclass),
					pred2.out=as.double(pred2), as.integer(verbose))


		mcmc.out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
		mcmc.out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
		mcmc.out$beta <- matrix(C.out$beta.out, nrow=nout, byrow=TRUE)
		mcmc.out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
		mcmc.out$like <- matrix(C.out$like.out, nrow=nout, byrow=TRUE)
		mcmc.out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
		mcmc.out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
		mcmc.out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
		mcmc.out$pred2 <- matrix(C.out$pred2.out, nrow=nout, byrow=TRUE)
		mcmc.out$mu0 <- C.out$mu0.out
		mcmc.out$sig20 <- C.out$sig20.out
		mcmc.out$nclus <- C.out$nclus.out
		mcmc.out$WAIC <- C.out$WAIC.out
		mcmc.out$lpml <- C.out$lpml.out
	}

	if(!psm){

		mcmc.out <- NULL

		mu <- sig2 <- Si <- like <- ispred <- matrix(1,nrow=nout,ncol=nobs)
		mu0 <- sig20 <- nclus <- rep(1,nout)
		ppred <- predclass <- pred2 <- matrix(1, nrow=nout, ncol=npred)
		WAIC <- lpml <- rep(1,1)

		C.out <- .C("mcmc_missing",
              	as.integer(draws), as.integer(burn), as.integer(thin),
              	as.integer(nobs),as.integer(ncon), as.integer(ncat),
              	as.integer(Cvec), as.integer(PPM), as.double(M),
              	as.integer(gcontype), as.integer(gcattype),
              	as.double(y),
 	            as.double(t(Xcon)),as.integer(t(Xcat)),
               	as.integer(t(Mcon)), as.integer(t(Mcat)),
             	as.integer(npred),
              	as.double(t(Xconp)),as.integer(t(Xcatp)),
              	as.integer(t(Mconp)), as.integer(t(Mcatp)),
              	as.integer(calibrate), as.integer(consim),as.double(alpha),
              	as.double(simParms),  as.double(modelPriors),
              	as.double(mh),
              	mu.out=as.double(mu), sig2.out=as.double(sig2),
              	mu0.out=as.double(mu0), sig20.out=as.double(sig20),
              	Si.out=as.integer(Si), nclus.out=as.integer(nclus),
              	like.out=as.double(like), WAIC.out=as.double(WAIC),
              	lpml.out=as.double(lpml),ispred.out=as.double(ispred),
              	ppred.out=as.double(ppred),predclass.out=as.integer(predclass),
		            pred2.out=as.double(pred2), as.integer(verbose))


 		mcmc.out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
  		mcmc.out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  		mcmc.out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
  		mcmc.out$like <- matrix(C.out$like.out, nrow=nout, byrow=TRUE)
  		mcmc.out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
  		mcmc.out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
  		mcmc.out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
  		mcmc.out$pred2 <- matrix(C.out$pred2.out, nrow=nout, byrow=TRUE)
  		mcmc.out$mu0 <- C.out$mu0.out
  		mcmc.out$sig20 <- C.out$sig20.out
  		mcmc.out$nclus <- C.out$nclus.out
  		mcmc.out$WAIC <- C.out$WAIC.out
  		mcmc.out$lpml <- C.out$lpml.out

	}

	out$Mmat <- Mmat
	if(psm){
	  out$gvec <- gvec_all
	  out$nb <- nb_all
	}
	out$mcmc <- mcmc.out
	out


}
