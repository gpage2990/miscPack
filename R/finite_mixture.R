gaussian_mixture <- function(y, Xmat=NULL, N, m, v, a, b, alpha, niter, nburn, nthin){

  n <- length(y)
  p <- ncol(Xmat)
  if(is.null(Xmat)){
    run <- .Call("GIBBS",
                    as.double(y), as.integer(n), as.integer(N),
                    as.double(m), as.double(v), as.double(a), as.double(b),
			  	          as.double(alpha),
                    as.integer(niter), as.integer(nburn), as.integer(nthin))
  } else {
    m <- 0;
    run <- .Call("GIBBS_REG",
                 as.double(y), as.double(t(Xmat)),
                 as.integer(n), as.integer(N), as.integer(p),
                 as.double(m), as.double(v), as.double(a), as.double(b),
                 as.double(alpha),
                 as.integer(niter), as.integer(nburn), as.integer(nthin))
  }

  run
}



dp_mixture <- function(y, N, m, v, a, b, alpha,
                       ndens_y=0,
                       subsample=FALSE, ss_size=100,
                       niter, nburn, nthin){

  n <- length(y)
  ygrid <- seq(min(y)-1, max(y)+1, length=ndens_y)

  run <- .Call("DP_GIBBS",
               as.double(y), as.integer(n), as.integer(N),
               as.double(m), as.double(v), as.double(a), as.double(b),
               as.double(alpha), as.integer(ndens_y), as.double(ygrid),
               as.integer(subsample), as.integer(ss_size),
               as.integer(niter), as.integer(nburn), as.integer(nthin))


  run
}




informed_mixture <- function(y, K=25,
                             alpha_prior_type="centered",
                             alpha_prior_dist="pc",
                             mu_sigma_prior = 1,
                             mylambda = NULL,
                             basemodel=0,
                             U=4,
                             alpha1_val = 1.0,
                             alpha2_val = 1e-5,
                             update_alpha1=TRUE,
                             update_alpha2=FALSE,
                             m0=0, s20=1, A=10, A0=10,
                             a0=1, b0=1,
                             k0=0.01, nu0 = 4,
                             a_gam=1,  b_gam=1,
                             a1_gam=1, b1_gam=1,
                             a2_gam=1, b2_gam=1,
                             tail.prob=0.01,
                             alpha_grid=NULL, alpha_density=NULL,
                             alpha0=1e-5, alphaMax=1,
                             ndens_y=0, hierarchy="NO",
                             niter, nburn, nthin){
  # note that alpha0 and alphaMax are needed only if alpha_grid and alpha_density
  # are not provided.

  # There are three prior specifications I am considering for (muk, sigma2k)
  # 1 - (muk, sigmak) ~ N(mu0, sigma20) x UN(0, A)
  # 2 - (muk, sigma2k) ~ N(mu0, sigma20) x IG(a0, b0)
  # 3 - (muk, sigma2k) ~ N(mu0, 1/k0sigma2k) x IG(1/2nu0, 1/2(nu0/s20)) This one matches Rafaelle

  n <- nrow(cbind(y))
  ygrid <- seq(min(y)-1, max(y)+1, length=ndens_y)

  if(!is.null(mylambda)) cat("employing user supplied lambda", mylambda, "\n")

  if(alpha_prior_type=="sparse" & alpha_prior_dist=="pc"){
    if(is.null(alpha_grid)){
      if(basemodel==0){
        alpha_grid <- exp(seq(log(alpha0), log(alphaMax), length.out = 1e4))

        if(is.null(mylambda)){
          cat("working on finding lambda", "\n")
          mylambda <- lambda_finder(U=U, tail.prob=tail.prob, alpha0=alpha0, agrid=alpha_grid,
                                     K=K, basemodel=0, n.obs=n, n.samples=1e4,
                                     alphamax=alphaMax,truncation=TRUE)

          cat("found lambda = ", mylambda, "\n")
        }
        alpha_density <- pcprior.a2(agrid, lam=mylambda, K=K,
                                              alpha0=alpha0, agrid=agrid,
                                              alphamax=alphaMax, TR=TRUE)
      }
      # This is obsolete as we want to the centered prior specification
      if(basemodel==1){
        alpha_grid <- exp(seq(log(0+1e-15), log(1), length.out = 10000))
#       lambda <- lambda_givenU(U=U,tail.prob=tail.prob,K=K,alpha0=1,
#                                n.obs=n,n.samples=10000,agrid=alpha_grid,
#                                interval=range(alpha_grid))
        cat("working on finding lambda", "\n")
        if(is.null(mylambda)){
          mylambda <- lambda_finder(U=U, tail.prob=tail.prob, alpha0=1, agrid=alpha_grid,
                                    K=K, basemodel=1, n.obs=n, n.samples=10000,
                                    alphamax=1e-15,truncation=TRUE)
          cat("found lambda", "\n")
          cat("lambda = ", mylambda, "\n")
        }
        alpha_density <- pcprior.a2(a.eval=agrid,
                                        lam=mylambda,
                                        K=K,
                                        alpha0=1,
                                        agrid=alpha_grid,
                                        alphamax=1e-15,
                                        TR=TRUE)
      }
    }
  }

  if(alpha_prior_type == "centered" & alpha_prior_dist=="pc"){
    alpha_grid <- alpha_density <- 1
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1

    if(update_alpha1==TRUE){

      # set initial search interval for lambda
      lam.ini <- c(1e-3, 3)
      ndens_alpha1 <- 1e3
      a1.max <- U-1e-3  # the max value here is U.  We just move below for numerical
      alpha1_grid <- exp(seq(log(1e-5), log(a1.max), length.out = 1e4))
      alpha_grid <- alpha1_grid
      U.lwr <- U - 1
      if(is.null(mylambda)){
        cat("working on finding lambda", "\n")
        mylambda <- miscPack:::lambda_finder_pcprior.asym.a1(U = U,
                                                  U.lwr = U.lwr,
                                                  tail.prob = tail.prob,
                                                  K = K,
                                                  alpha1.base=U,
                                                  alpha2.base=1e-5,
                                                  alpha2.fixed=alpha2_val,
                                                  agrid=alpha1_grid,
                                                  n.obs=n,
                                                  n.samples=1e4,
                                                  TR=TRUE,
                                                  lambda.interval=lam.ini,
                                                  tol=1e-5)

        cat("found lambda = ", mylambda, "\n")
      }
      alpha1_density <- miscPack:::pcprior.asym.a1(a.eval=alpha1_grid, lam=mylambda,
                                                   U=U, K=K, alpha1.base = U,
                                                   alpha2.base = 1e-5,
                                                   alpha2.fixed = alpha2_val,
                                                   agrid=alpha1_grid,
                                                   TR=TRUE)
      alpha_density <- alpha1_density

    }
    if(update_alpha2==TRUE){

    }
  }
  if(alpha_prior_dist=="gamma"){
    alpha_grid <- alpha_density <- 1
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1
  }
  if(alpha_prior_type=="sparse"){
    update_alpha1 <- update_alpha2 <- 0
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1
  }
  ndens_alpha <- length(alpha_grid)
  ndens_alpha1 <- length(alpha1_density)
  ndens_alpha2 <- length(alpha2_density)

#  print(summary(alpha_grid))
#  print(summary(alpha1_grid))
#  print(summary(alpha2_grid))

  if(ncol(cbind(y)) > 1){
    p <- ncol(cbind(y))
    #  plot(alpha_grid, alpha_density, type='l')
    M0 <- m0
    Sigma0 <- s20*diag(p)
    K0 <- k0*diag(p)
    run <- .Call("MVT_IFMM_GIBBS",
                 as.double(t(y)), as.integer(n), as.integer(p), as.integer(K),
                 as.integer(ifelse(alpha_prior_type=="sparse",1 ,2)),
                 as.integer(ifelse(alpha_prior_dist=="pc",1 ,2)),
                 as.integer(mu_sigma_prior),
                 as.integer(U),
                 as.integer(basemodel),
                 as.double(alpha1_val), as.double(alpha2_val),
                 as.integer(update_alpha1), as.integer(update_alpha2),
                 as.double(ygrid), as.integer(ndens_y),
                 as.integer(ifelse(hierarchy=="NO", 0, 1)),
                 as.double(alpha_grid), as.double(alpha_density), as.integer(ndens_alpha),
                 as.double(alpha1_grid), as.double(alpha1_density), as.integer(ndens_alpha1),
                 as.double(alpha2_grid), as.double(alpha2_density), as.integer(ndens_alpha2),
                 as.double(M0), as.double(Sigma0),  as.integer(nu0), as.double(K0),
                 as.double(A), as.double(a0), as.double(b0),
                 as.double(a_gam), as.double(b_gam),
                 as.double(a1_gam), as.double(b1_gam),
                 as.double(a2_gam), as.double(b2_gam),
                 as.integer(niter), as.integer(nburn), as.integer(nthin))

  } else {

    #  plot(alpha_grid, alpha_density, type='l')
    run <- .Call("IFMM_GIBBS",
                 as.double(y), as.integer(n), as.integer(K),
                 as.integer(ifelse(alpha_prior_type=="sparse",1 ,2)),
                 as.integer(ifelse(alpha_prior_dist=="pc",1 ,2)),
                 as.integer(mu_sigma_prior),
                 as.integer(U),
                 as.integer(basemodel),
                 as.double(alpha1_val), as.double(alpha2_val),
                 as.integer(update_alpha1), as.integer(update_alpha2),
                 as.double(ygrid), as.integer(ndens_y),
                 as.integer(ifelse(hierarchy=="NO", 0, 1)),
                 as.double(alpha_grid), as.double(alpha_density), as.integer(ndens_alpha),
                 as.double(alpha1_grid), as.double(alpha1_density), as.integer(ndens_alpha1),
                 as.double(alpha2_grid), as.double(alpha2_density), as.integer(ndens_alpha2),
                 as.double(m0), as.double(s20), as.double(A), as.double(A0),
                 as.double(a0), as.double(b0), as.double(nu0), as.double(k0),
                 as.double(a_gam), as.double(b_gam),
                 as.double(a1_gam), as.double(b1_gam),
                 as.double(a2_gam), as.double(b2_gam),
                 as.integer(niter), as.integer(nburn), as.integer(nthin))


  }

  if(alpha_prior_type=="centered") run$alpha <- NULL
  if(alpha_prior_type!="centered"){run$alpha1 <- NULL; run$alpha2 <- NULL}
  if(hierarchy=="NO"){run$mu0 <- NULL; run$sigma20 <- NULL}
  run$kp <- apply(run$z,1,function(x)length(unique(x)))
  run$ygrid <- ygrid
  if(alpha_prior_dist== "pc"){
    run$lambda_val <- mylambda
    run$alpha_grid <- alpha_grid
    run$alpha_density <- alpha_density
  }
  run
}




pspline_mixture <- function(y, t, ids, K,
                            alpha_prior_type="centered",
                            alpha_prior_dist="pc",
                            mylambda = NULL,
                            basemodel=0,
                            U=4,
                            alpha1_val = 1.0,
                            alpha2_val = 1e-5,
                            update_alpha1=TRUE,
                            update_alpha2=FALSE,
                            tail.prob=0.9,
                            ndx = 10, q = 3, ndens_y = 1,
                            A=1, Akap=1,
                            U_tau=1, a_tau=0.1,
                            U_omega=1, a_omega=0.1,
                            a_gam=1, b_gam = 1,
                            a1_gam =1, b1_gam=1,
                            a2_gam=1, b2_gam=2,
                            alpha_grid, alpha_density,
                            niter, nburn, nthin){
  # A sampling of the inputs are the following
  # U_tau and a_tau help ilicit prior for tau2k.  U_tau is an upperbound
  # for tauk and a_tau is the probability of exceeding it.  i.e., we have
  # Pr(tauk > U_tau) = a_tau
  # Then from U_tau and a_tau, e_tau (parameter of the Gumbel distributio)
  # is derived using e_tau = -log(a_tau)/U_tau

  e_tau = -log(a_tau)/U_tau
  e_omega = -log(a_omega)/U_omega
  cat("e_tau = ", e_tau, "\n");
  cat("e_omega = ", e_omega, "\n");
  n <- length(unique(ids))
  nobs <- tapply(ids, ids, length)


  # Best to feed in each subjects basis created from Cote's code rather than
  # creating code in C. Note that the G matrix depends on the balanced argument

  balanced <- ifelse(length(unique(nobs))==1, 1, 0)

  # balanced 1 - All subjects have same number of measurements and time at which
  #              they were measure and so have same design matrix
  #		       0 - Subjects do not have same number of measurements and they are
  #		           measured at different time points so each subject needs their
  #              own design matrix

  # Unique time points
  tt <- sort(unique(t))

  # Create the basis
  G <- bbase(tt, ndx = ndx, bdeg = q)

  # Smoothing matrix
  D <- diff(diag(ncol(G)), differences = 1)
  S <- crossprod(D)
#  S[1,1] <- 2 # to make invertible
  S <- S + 0.01*diag(ncol(S))
#  print(solve(S))
  nb <- ncol(S)

  Xmat <- G
  weights <- apply(Xmat,2,sum)/sum(apply(Xmat,2,sum))

  # As of now, I am not permitting unbalanced designs.
  if(!balanced) {
    Xmat <- NULL
    cnobs <- 0
    for(j in 1:n){
      ttmp <- df$t[(cnobs+1):(cnobs + nobs[j])]
      bsb <- miscPack::predict.bbase(G, ttmp)
      Xmat <- c(Xmat, c(t(bsb)))
      cnobs <- cnobs + nobs[j]

    }
  }

  if(alpha_prior_type=="sparse" & alpha_prior_dist=="pc"){
    if(is.null(alpha_grid)){
      if(basemodel==0){
        alpha_grid <- exp(seq(log(alpha0), log(alphaMax), length.out = 1e4))

        if(is.null(mylambda)){
          cat("working on finding lambda", "\n")
          mylambda <- lambda_finder(U=U, tail.prob=tail.prob, alpha0=alpha0, agrid=alpha_grid,
                                  K=K, basemodel=0, n.obs=n, n.samples=1e4,
                                  alphamax=alphaMax,truncation=TRUE)

          cat("found lambda = ", mylambda, "\n")
        }
        alpha_density <- pcprior.a2(agrid, lam=mylambda, K=K,
                                    alpha0=alpha0, agrid=agrid,
                                    alphamax=alphaMax, TR=TRUE)
      }

      # This is obsolete as we went to the centered prior specification
      if(basemodel==1){
        alpha_grid <- exp(seq(log(0+1e-15), log(1), length.out = 10000))
        #       lambda <- lambda_givenU(U=U,tail.prob=tail.prob,K=K,alpha0=1,
        #                                n.obs=n,n.samples=10000,agrid=alpha_grid,
        #                                interval=range(alpha_grid))
        if(is.null(mylambda)){
          cat("working on finding lambda", "\n")
          mylambda <- lambda_finder(U=U, tail.prob=tail.prob, alpha0=1, agrid=alpha_grid,
                                  K=K, basemodel=1, n.obs=n, n.samples=10000,
                                  alphamax=1e-15,truncation=TRUE)
          cat("found lambda", "\n")
          cat("lambda = ", mylambda, "\n")
        }
        alpha_density <- pcprior.a2(a.eval=agrid,
                                    lam=mylambda,
                                    K=K,
                                    alpha0=1,
                                    agrid=alpha_grid,
                                    alphamax=1e-15,
                                    TR=TRUE)
      }
    }
  }

  if(alpha_prior_type == "centered" & alpha_prior_dist=="pc"){
    alpha_grid <- alpha_density <- 1
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1

    if(update_alpha1==TRUE){

      # set initial search interval for lambda
      lam.ini <- c(1e-3, 3)
      ndens_alpha1 <- 1e3
      a1.max <- U-1e-3  # the max value here is U.  We just move below for numerical
      alpha1_grid <- exp(seq(log(1e-5), log(a1.max), length.out = 1e4))
      U.lwr <- U - 1
      if(!is.null(mylambda)) cat("employing user supplied lambda", mylambda, "\n")
      if(is.null(mylambda)){
        cat("working on finding lambda", "\n")
        mylambda <- miscPack:::lambda_finder_pcprior.asym.a1(U = U,
                                                             U.lwr = U.lwr,
                                                             tail.prob = tail.prob,
                                                             K = K,
                                                             alpha1.base=U,
                                                             alpha2.base=1e-5,
                                                             alpha2.fixed=alpha2_val,
                                                             agrid=alpha1_grid,
                                                             n.obs=n,
                                                             n.samples=1e4,
                                                             TR=TRUE,
                                                             lambda.interval=lam.ini,
                                                             tol=1e-5)

        cat("found lambda = ", mylambda, "\n")
      }
      alpha1_density <- miscPack:::pcprior.asym.a1(a.eval=alpha1_grid, lam=mylambda,
                                                   U=U, K=K, alpha1.base = U,
                                                   alpha2.base = 1e-5,
                                                   alpha2.fixed = alpha2_val,
                                                   agrid=alpha1_grid,
                                                   TR=TRUE)
      # These two are necessary to find the initial
      # log_ao_density
      alpha_density <- alpha1_density
      alpha_grid <- alpha1_grid
      Kplus.draw <- impliedpcprior.asym.a1.Kplus(lam=mylambda,
                                                 U = U,
                                                 K = K,
                                                 alpha1.base=U,
                                                 alpha2.base=1e-5,
                                                 alpha2.fixed=1e-5,
                                                 agrid=alpha1_grid,
                                                 n.obs = n,
                                                 n.samples = 1e4,
                                                 TR=TRUE)$Kplus
      cat("Pr(K+ < U)", mean(Kplus.draw<U), "\n")
      cat("5 number summary of implied prior of K+", "\n", summary(Kplus.draw), "\n")

    }
    if(update_alpha2==TRUE){

    }
  }
  if(alpha_prior_dist=="gamma"){
    alpha_grid <- alpha_density <- 1
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1
  }
  if(alpha_prior_type=="sparse"){
    update_alpha1 <- update_alpha2 <- 0
    alpha1_grid <- alpha1_density <- 1
    alpha2_grid <- alpha2_density <- 1
  }
  ndens_alpha <- length(alpha_grid)
  ndens_alpha1 <- length(alpha1_density)
  ndens_alpha2 <- length(alpha2_density)


  # print(summary(alpha1_grid))
  # print(summary(alpha1_density))

  geom.mean <- function(x) exp(mean(log(diag(ginv(x)))))
  gm <- geom.mean(S)
  S.scaled <- gm*S

  nobs <- nobs[1]
  run <- .Call("PSPLINE_IFMM_GIBBS",
          as.double(y), as.double(t(Xmat)), as.integer(n), as.integer(nobs),           #4
          as.double(t(S)), as.integer(nb), as.integer(K),                              #3
          as.integer(ifelse(alpha_prior_type=="sparse",1 ,2)),                         #1
          as.integer(ifelse(alpha_prior_dist=="pc",1 ,2)),                             #1
          as.integer(U),                                                               #1
          as.integer(basemodel), as.double(weights),                                   #2
          as.double(alpha1_val), as.double(alpha2_val),                                #2
          as.integer(update_alpha1), as.integer(update_alpha2),                        #2
          as.double(alpha_grid), as.double(alpha_density), as.integer(ndens_alpha),    #3
          as.double(alpha1_grid), as.double(alpha1_density), as.integer(ndens_alpha1), #3
          as.double(alpha2_grid), as.double(alpha2_density), as.integer(ndens_alpha2), #3
          as.double(A), as.double(Akap),                                               #2
          as.double(e_tau), as.double(e_omega),                                        #2
          as.double(a_gam), as.double(b_gam),                                          #2
          as.double(a1_gam), as.double(b1_gam),                                        #2
          as.double(a2_gam), as.double(b2_gam),                                        #2
          as.integer(niter), as.integer(nburn), as.integer(nthin))                     #3
                                                                                       #38 - total

  run$kp <- apply(run$z, 1, function(x)length(unique(x)))
  run$Xmat <- Xmat
  run$U <- U
  if(alpha_prior_dist== "pc") run$lambda_pc_val <- mylambda
  run$alpha_grid <- alpha_grid

  run
}


pspline_mixture_old <- function(y, t, ids, K,
                            alpha_prior, basemodel, alpha0,
                            ndx = 10, q = 3,
                            A=1, U_tau=1, a_tau=0.1, U_omega=1, a_omega=0.1,
                            a_gam=1,
                            alpha_grid, alpha_density,
                            niter, nburn, nthin){
  # A sampling of the inputs are the following
  # U_tau and a_tau help ilicit prior for tau2k.  U_tau is an upperbound
  # for tauk and a_tau is the probability of exceeding it.  i.e., we have
  # Pr(tauk > U_tau) = a_tau
  # Then from U_tau and a_tau, e_tau (parameter of the Gumbel distribution)
  # is derived using e_tau = -log(a_tau)/U_tau

  e_tau = -log(a_tau)/U_tau
  e_omega = -log(a_omega)/U_omega
  cat("e_tau = ", e_tau, "\n");
  cat("e_omega = ", e_omega, "\n");
  n <- length(unique(ids))
  nobs <- tapply(ids, ids, length)
  ndens_alpha <- length(alpha_grid)

  # Best to feed in each subjects basis created from Cote's code rather than
  # creating code in C. Note that the Hmat matrix depends on the balanced argument

  balanced <- ifelse(length(unique(nobs))==1, 1, 0)

  # balanced 1 - All subjects have same number of measurements and time at which
  #              they were measure and so have same design matrix
  #		       0 - Subjects do not have same number of measurements and they are
  #		           measured at different time points so each subject needs their
  #              own design matrix

  # Unique time points
  tt <- sort(unique(t))

  # Create the basis
  G <- bbase(tt, ndx = ndx, bdeg = q)

  # Smoothing matrix
  D <- diff(diag(ncol(G)), differences = 1)
  S <- crossprod(D)
  S[1,1] <- 2 # to make invertible essentially assumes beta1 ~ N(0,tau2)

  nb <- ncol(S)
  Xmat <- G

  # As of now, I am not permitting unbalanced designs.
  if(!balanced) {
    Xmat <- NULL
    cnobs <- 0
    for(j in 1:n){
      ttmp <- df$t[(cnobs+1):(cnobs + nobs[j])]
      bsb <- predict.bbase(G, ttmp)
      Xmat <- c(Xmat, c(t(bsb)))
      cnobs <- cnobs + nobs[j]

    }
  }

  nobs <- nobs[1]
  run <- .Call("PSPLINE_IFMM_GIBBS",
               as.double(y), as.double(t(Xmat)), as.integer(n), as.integer(nobs),
               as.double(t(S)), as.integer(nb), as.integer(K),
               as.integer(alpha_prior), as.integer(basemodel), as.double(alpha0),
               as.double(alpha_grid), as.double(alpha_density), as.integer(ndens_alpha),
               as.double(A), as.double(e_tau), as.double(e_omega), as.double(a_gam),
               as.integer(niter), as.integer(nburn), as.integer(nthin))


  run
}



