
#-----------------------------------------------------------------------------------------------------


## R wrapper that accesses C code to fit curve cluster model

varpar_curve_ppmx <- function(y, z, subject,
                       Xcon=NULL,Xcat=NULL,
                       Xconp=NULL,Xcatp=NULL,
                       PPM, M,
                       q=3, rw_order=2, balanced=1,
                       nknots,npredobs,
                       modelPriors,
                       similarity_function=1,
                       consim,
                       simParms,
                       mh=c(1,1,1,1),
                       draws=1100,burn=100,thin=1){



  ## consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  ## consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)

  out <- NULL
  nobs <- tapply(subject, subject, length)
  nsubject <- length(nobs)

  if(is.null(Xcon)){
    ncon <- 0; Xcon <- cbind(rep(0,nsubject))
  }else{
    ncon <- ncol(Xcon)
  }
  #  cat("ncon = ", ncon, "\n")

  if(is.null(Xcat)){
    ncat <- 0; Xcat <- cbind(rep(0,nsubject))
    Cvec <- 0
  }else{
    ncat <- ncol(Xcat)
    Cvec <- apply(Xcat,2,function(x)length(unique(x)))
  }
  # cat("ncat = ", ncat, "\n")


  if(is.null(Xconp) & is.null(Xcatp)){
    npred <- 0
  } else {
    npred <- ifelse(is.null(Xconp),nrow(Xcatp), nrow(Xconp))
  }

  #  cat("npred = ", npred, "\n")

  if(is.null(Xconp)){
    Xconp <- cbind(rep(0,npred))
  }
  if(is.null(Xcatp)){
    Xcatp <- cbind(rep(0,npred))
  }


  ndx <- nknots  # number of segments (inner knots)

  # Best to feed in each subjects basis created from Cote's code rather than
  # creating code in C. Not the Hmat matrix depends on the balanced argument
  #
  # balanced 1 - All subjects have same number of measurements and time at which
  #              they were measure and so have same design matrix
  #		       0 - Subjects do not have same number of measurements and they are
  #		           measured at different time points so each subject needs their
  #              own design matrix

  # Unique time points
  tt <- sort(unique(z))

  # Create the basis
#  G <- bbase(tt, ndx = ndx, bdeg = q)

  B <- scaleGMRF::bspline(x = tt, K = ndx, m = 0, M = 1)

  # Linear decomposition
  Z <- cbind(
    rep(1, ndx),
    seq_len(ndx))
  A <- qr.Q(qr(Z), complete = TRUE)

  Hmat <- B%*%A

  # For prediction
  # Don't worry about prediction at this point
  #trate <- mean(diff(z))
  #tpred <- seq(min(z), max(z)+npredobs*trate, length = length(tt)+npredobs)
  #Gp <- predict.bbase(G, tpred)
  #Hmat_pred <- cbind(1,tpred,Gp)
  Hmat_pred <- Hmat


  Q_rw <- as.matrix(precmat.RWn(ndx, order=rw_order))
  Q_rw <- Q_rw * mean(diag(gen_inv(Q_rw, rank_def = rw_order))) #scaling

  #Since theta=A%*%beta, and beta=t(A)%*%theta (see line 23),
  # then the T_matrix (without the variance parameters) is:
  T_matrix <- t(A)%*%gen_inv(Q_rw,rank_def=rw_order)%*%A

  K <- T_matrix[3:ndx,3:ndx] # not sure about the scaling of this matrix
  Kinv <- gen_inv(K) # proper precision matrix of RW2 in the new space


#  D <- diff(diag(ncol(G)), differences = rw_order)

  nb <- ncol(Hmat)

  # Will have to deal with unbalanced designs later
#  if(!balanced) {
#    Hmat <- NULL
#    Hmat_pred <- NULL
#    cnobs <- 0
#    for(j in 1:nsubject){
#      ttmp <- z[(cnobs+1):(cnobs + nobs[j])]
#      bsb <- predict.bbase(G, ttmp)
#      Hmat <- c(Hmat, c(t(bsb)))
#
#      ttmp <- c(z[(cnobs+1):(cnobs + nobs[j])],
#                (z[(cnobs + nobs[j])]+1):(z[(cnobs + nobs[j])]+npredobs))
#      bsb <- predict.bbase(G, ttmp)
#      Hmat_pred <- c(Hmat_pred, c(t(bsb)))
#
#      cnobs <- cnobs + nobs[j]
#
#    }
#  }



  nout <- (draws-burn)/thin;

  Si <- llike  <- lamint <- lamslope <- lamcurve <- matrix(1,nrow=nout,ncol=nsubject)
  fitlines <- matrix(0, nrow=nout, ncol=sum(nobs))
  sig2 <- matrix(1,nrow=nout,ncol=nsubject)
  iCPO <- tau2int <- tau2slope <- tau2curve <- matrix(1,nrow=nout,ncol=nsubject)

  beta <- theta <-  matrix(0, nrow = nout, ncol = nb*nsubject)
  mu <- matrix(0, nrow=nout, ncol=nb)
  nclus <- rep(1, nout)
  predclass <- matrix(1, nrow=nout, ncol=npred)
  ppred <- matrix(0, nrow=nout, ncol=npred*npredobs)
  predlines <- matrix(0, nrow=nout, ncol=sum(nobs+npredobs))
  lpml <- WAIC <- rep(0,1)


  C.out <- .C("varpar_curvecluster",
              as.integer(draws), as.integer(burn),as.integer(thin),
              as.integer(nsubject), as.integer(nobs),
              as.double(t(y)), as.double(z),
              as.double(t(Kinv)), as.integer(nb), as.double(t(K)),
              as.integer(ncon), as.integer(ncat), as.integer(Cvec),
              as.double(t(Xcon)), as.integer(t(Xcat)),
              as.integer(PPM), as.double(M),
              as.integer(similarity_function), as.integer(consim),
              as.integer(npred),as.integer(npredobs),
              as.double(t(Xconp)),as.integer(t(Xcatp)),
              as.double(simParms),
              as.double(modelPriors),
              as.double(t(Hmat)), as.integer(balanced), as.double(mh),
              as.double(t(Hmat_pred)),
              beta.out= as.double(beta),theta.out= as.double(theta),
              mu.out = as.double(mu), sig2.out= as.double(sig2),
              lam_int.out= as.double(lamint), lam_slope.out= as.double(lamslope),
              lam_curve.out= as.double(lamcurve),
              tau2_int.out= as.double(tau2int), tau2_slope.out= as.double(tau2slope),
              tau2_curve.out= as.double(tau2curve),
              Si.out= as.integer(Si),nclus.out=as.integer(nclus),
              ppred.out=as.double(ppred),
              predclass.out=as.integer(predclass),
              fitlines.out=as.double(fitlines),
              predlines.out=as.double(predlines),
              llike.out=as.double(llike),
              lpml.out=as.double(lpml),WAIC.out=as.double(WAIC))

  out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
  out$nclus <- matrix(C.out$nclus.out, nrow=nout, byrow=TRUE)
  out$beta <- array(C.out$beta.out, c(nsubject, nb, nout))
  out$theta <- array(C.out$theta.out, c(nsubject, nb, nout))
  out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
  out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  out$tau2_int <- matrix(C.out$tau2_int.out, nrow=nout, byrow=TRUE)
  out$tau2_slope <- matrix(C.out$tau2_slope.out, nrow=nout, byrow=TRUE)
  out$tau2_curve <- matrix(C.out$tau2_curve.out, nrow=nout, byrow=TRUE)
  out$lam_int <- matrix(C.out$lam_int.out, nrow=nout, byrow=TRUE)
  out$lam_slope <- matrix(C.out$lam_slope.out, nrow=nout, byrow=TRUE)
  out$lam_curve <- matrix(C.out$lam_curve.out, nrow=nout, byrow=TRUE)
  out$llike <- matrix(C.out$llike.out, nrow=nout, byrow=TRUE)
  out$fitted <- matrix(C.out$fitlines.out, nrow=nout, byrow=TRUE)
  #	out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
  #	out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
  out$pred_curve <- matrix(C.out$predlines.out, nrow=nout, byrow=TRUE)
  out$WAIC <- C.out$WAIC.out
  out$lpml <- C.out$lpml.out


  out$Hmat <- Hmat
  out$Hmat_pred <- Hmat_pred
  out$number.of.basis <- nb
  out$K <- K
  out
}



