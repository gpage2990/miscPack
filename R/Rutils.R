# Functions needed to fit Derivative Curve Clustering Methodology

tpower <- function(x, t, p) {
  # Function for truncated p-th power function
  return((x - t) ^ p * (x > t))
}
# Function for constructing a B-spline basis
# x: covariate
# ndx: number of segments (related to the number of knots)
# bdeg: degree of the polynomial (e.g. 3: Cubic B-splines)
bbase <- function(x, ndx, bdeg = 3, eps = 1e-5) {
  xl = min(x)
  xr = max(x)
  dx <- (xr - xl)/ndx
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B[B < eps] <- 0
  attr(B,"knots") <- knots
  attr(B,"bdeg") <- bdeg
  attr(B,"eps") <- eps
  class(B) <- c("bbase")
  B
}
# Prediction function
predict.bbase <- function(object, newx) {
  knots <- attr(object,"knots")
  bdeg <- attr(object,"bdeg")
  eps <- attr(object,"eps")

  dx <- diff(knots)[1]
  P <- outer(newx, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1)*dx^bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B[B < eps] = 0
  B
}


# Functions that are useful
relabel <- function(z) {
  # Get the unique labels in the input vector
  unique_labels <- unique(z)

  # Create a mapping from old labels to new labels
  old_to_new_labels <- 1:length(unique_labels)
  names(old_to_new_labels) <- unique_labels

  # Use the match() function to update labels
  new_labels <- old_to_new_labels[match(z, unique_labels)]

  match(z, unique_labels)
}


# in Eilers 1996 viene useto il package splines e la seguente funzione
bspline = function(x, xl, xr, ndx, bdeg) {
  dx = (xr - xl) / ndx
  knots = seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B = spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
  B
}


# Kullback-Leibler distance between two Symmetric Dirichlet:
# KLD[Dir(alpha)|Dir(alpha0)]
kld <- function(alpha, alpha0){

  c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
    lgamma(sum(alpha0)) - sum(lgamma(alpha0)))

  c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))

  c1 + c2
}


# draw a sample of a using MC and a approximated
# KLD inverse function

# alpha0 is the value of alpha at the base model
# alphamax is the value of alpha that is most far apart from alpha0 (most flex model)
draw.a <- function(lam,
                   K,
                   basemodel=0,
                   alpha0,
                   n.samples,
                   agrid,
                   alphamax,
                   TR=T){
  a <- agrid
  if (basemodel==0) b.texp <- sqrt(2*kld(rep(alphamax,K),rep(alpha0,K)))
  if (basemodel==1) b.texp <- Inf # sqrt(2*kld(rep(alphamax,K),rep(alpha0,K)))

  # evaluate the kullback-Leibler distance
  kld_val <- sapply(a, function(x) sqrt(2*kld(rep(x,K),rep(alpha0,K))))
  #a.to.kld <- splinefun(a, kld_val)
  kld.to.a <- splinefun(kld_val, a)

  if(TR) {
    sample.kld <- rtexp(n.samples, rate = lam, a=0, b=b.texp)
  } else {
    sample.kld <- rexp(n.samples, rate = lam)
  }

  sample.a <- kld.to.a(sample.kld)
  sample.a[sample.a<0] <- 0
  sample.a
}

# sample from the implied distribution of K+
draw.Kplus <- function(lam, K,
                       basemodel=0,
                       alpha0,
                       n.obs,
                       n.samples,
                       agrid,
                       alphamax,
                       TR=T){
  sample.a <- draw.a(lam=lam, K=K,
                     basemodel = basemodel,
                     alpha0=alpha0,
                     n.samples=n.samples,
                     agrid=agrid,
                     alphamax = alphamax,
                     TR=TR)
  sample.w <- t(sapply(sample.a, function(x) rdirichlet(1, rep(x, K))))
  sample.z <- matrix(NA, nrow=n.samples, ncol=K)
  for(i in 1:n.samples) {
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus, a=sample.a)
}


# sample from the implied distribution of K+
draw.Kplus.given.a <- function(a.eval, K,
                               n.obs, n.samples,
                               agrid){
  sample.a <- a.eval
  sample.w <- t(sapply(sample.a, function(x) rdirichlet(1, rep(x, K))))
  sample.z <- matrix(NA, nrow=n.samples, ncol=K)
  for(i in 1:n.samples) {
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus, a=sample.a)
}

# BUILD splinefun pc prior density for a \in (0,1]
pcprior.a <- function(lam=1,
                      K=5,
                      alpha0,
                      agrid,
                      alphamax,
                      TR=T){

  # agrid is a sequence of values for which kullback-Leibler distance will be evaluated
  # evaluate the kullback-Leibler distance
  kld_val <- sapply(agrid, function(x) kld(rep(x,K),rep(alpha0,K)))


  # expontential density (PC prior is based on exponential of KL dist(ance)
  if(TR){
    dmax <- sqrt(2*kld(rep(alphamax,K),rep(alpha0,K)))
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of function KLD(a) a < 1
  a.to.kld <- splinefun(agrid, sqrt(2*kld_val))

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(agrid, fexp(lambda = lam, d=a.to.kld(agrid)) * abs(a.to.kld(agrid, deriv=1)))

  return(pcprior)

}

# EVALUATE pc prior density for a \in (0,1]
pcprior.a2 <- function(a.eval, lam=1, K=5,
                       alpha0,
                       agrid,
                       alphamax,
                       TR=T){

  # agrid is a sequence of values for which kullback-Leibler distance will be evaluated
  # evaluate the kullback-Leibler distance
  kld_val <- sapply(agrid, function(x) kld(rep(x,K),rep(alpha0,K)))


  # expontential density (PC prior is based on exponential of KL dist(ance)
  if(TR){
    dmax <- sqrt(2*kld(rep(alphamax,K),rep(alpha0,K)))
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of function KLD(a) a < 1
  a.to.kld <- splinefun(agrid, sqrt(2*kld_val))

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(agrid, fexp(lambda = lam, d=a.to.kld(agrid)) * abs(a.to.kld(agrid, deriv=1)))

  tmp <- pcprior(a.eval)
  return(ifelse(tmp < 0, 0, tmp))

}

# This produces the density of the induced prior on "log(alpha)"
pcprior.log_a <- function(log_a.eval, lam, K=5,
                          alpha0,
                          agrid,
                          alphamax,
                          TR=T){

  # evaluate the kullback-Leibler distance
  kld_val <- sapply(agrid, function(x) kld(rep(x,K),rep(alpha0,K)))

  if(TR){
    dmax <- sqrt(2*kld(rep(alphamax,K),rep(alpha0,K)))
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  a.to.kld <- splinefun(agrid, sqrt(2*kld_val)) # distance function, kld

  pcprior <- splinefun(agrid, fexp(lambda = lam, d=a.to.kld(agrid)) * abs(a.to.kld(agrid, deriv=1)))

  pcprior.spline <- splinefun(agrid,pcprior(agrid))

  pcprior.spline.logscale <- splinefun(log(agrid), pcprior.spline(agrid)*agrid)
  tmp <- pcprior.spline.logscale(log_a.eval)
  return(ifelse(tmp < 0, 0, tmp))
}


######### lambda finder #########

### old function (manual search)
# lambda_finder <- function(U, tail.prob, alpha0, agrid, K,
#                           basemodel=0,
#                           n.obs,
#                           n.samples,
#                           alphamax,
#                           truncation,
#                           lam_ini=2,
#                           lam_tun=0.01){
#   #  lam_vec <- seq(0.01, 2, by=0.01)
#   #  tmp <- numeric()
#   #  for(kk in 1:length(lam_vec)){
#   #    samp <- draw.Kplus(lam=lam_vec[kk], alpha0=alpha0, K=K, n.obs=n.obs, n.samples=n.samples,
#   #                       agrid=agrid, alphamax=alphamax, TR=truncation)$Kplus
#   #    tmp[kk] <- quantile(samp, 1-tail.prob)
#   #    if(kk > 1){
#   #      if(( (tmp[kk] - U < 1) & (tmp[kk-1]-tmp[kk]>0) ) |
#   #           (tmp[kk] - U > 1)  & (tmp[kk-1]-tmp[kk]>0)) break
#   #    }
#   #  }
#   #  mean(lam_vec[which(abs(tmp-U) == min(abs(tmp - U)))])
#
#   #  bmp <- rdirichlet(1000, rep(alpha0, K))
#   #  bmKp <- apply(sapply(1:1000, function(x) sample(1:K, n.obs, replace=TRUE, prob=bmp[x,])),
#   #                2, function(x) length(unique(x)))
#   #  bmq <- quantile(bmKp, 1-tail.prob)
#   #  if(bmq >= U & basemodel==0){
#   #    stop(paste0("The 1-",tail.prob," quantile of K+ for alpha0 = ", alpha0,
#   #                " and sample size of ", n.obs,
#   #                " is lower bounded at ", bmq, ".", "\n",
#   #               " consider using a smaller alpha0 value"))
#   #  }
#   if(basemodel==1){
#     lam_tmp <- lam_ini
#     qq <- quantile(draw.Kplus(lam=lam_tmp,
#                               alpha0=alpha0,
#                               K=K,
#                               basemodel=basemodel,
#                               n.obs=n.obs,
#                               n.samples=n.samples,
#                               agrid=agrid,
#                               alphamax=alphamax,
#                               TR=truncation)$Kplus, tail.prob)
#     while(abs(qq - U) != 0){
#       samp <- draw.Kplus(lam=lam_tmp,
#                          alpha0=alpha0,
#                          K=K,
#                          basemodel=basemodel,
#                          n.obs=n.obs,
#                          n.samples=n.samples,
#                          agrid=agrid,
#                          alphamax=alphamax,
#                          TR=truncation)$Kplus
#       qq <- quantile(samp, tail.prob)
#       if(qq > U){
#         lam_tmp <- lam_tmp - lam_tun
#       } else {
#         lam_tmp <- lam_tmp + lam_tun
#       }
#       if(lam_tmp == lam_tun) break
#     }
#   }
#   if(basemodel==0){
#     lam_tmp <- lam_ini
#     qq <- quantile(draw.Kplus(lam=lam_tmp, alpha0=alpha0, K=K,
#                               basemodel=basemodel,
#                               n.obs=n.obs,
#                               n.samples=n.samples,
#                               agrid=agrid,
#                               alphamax=alphamax,
#                               TR=truncation)$Kplus, 1-tail.prob)
#     while(abs(qq - U) != 0){
#       samp <- draw.Kplus(lam=lam_tmp, alpha0=alpha0, K=K,
#                          basemodel=basemodel,
#                          n.obs=n.obs,
#                          n.samples=n.samples,
#                          agrid=agrid,
#                          alphamax=alphamax,
#                          TR=truncation)$Kplus
#       qq <- quantile(samp, 1-tail.prob)
#       if(qq < U){
#         lam_tmp <- lam_tmp - lam_tun
#       } else {
#         lam_tmp <- lam_tmp + lam_tun
#       }
#       #      cat("qq = ", qq, "\n")
#       #      cat("lam_tmp = ", lam_tmp, "\n")
#       if(lam_tmp == lam_tun) break
#     }
#   }
#   lam_tmp
# }



### using optimize()
locate.previous <- function(x, U){
  res <- x-U
  max(x[res<0])
}

locate.next <- function(x, U){
  res <- x-U
  min(x[res>0])
}

refine.lambda.int <- function(y, x, U, K, basemodel){
  n <- length(x)
  if (basemodel==1){
    # locating U for basemodel=1
    if (all(x>U)) stop("initialize lambda.interval to a SMALLER lower bound")
    if (all(x[2:n] == K)) stop("initialize lambda.interval to a SMALLER upper bound")
    U.previous <- miscPack:::locate.previous(x, U)
    U.next <- miscPack:::locate.next(x, U)
    if(sum(x >=U.previous & x <= U.next)>0)
      refined.interval <- range(y[sum(x >=U.previous & x <= U.next)>0])
    else  stop("initialize lambda.interval to a BIGGER upper bound")
  } else {
    # locating U for basemodel=1
    if (locate.previous(x, U) < 0)  stop("initialize lambda.interval to a LARGER upper bound")
    if (locate.next(x, U) >= K)  stop("initialize lambda.interval to other values")
    if (all(x>U)) stop("initialize lambda.interval using a BIGGER upper bound")
    if (all(x[2:n] == K)) stop("initialize lambda.interval to a BIGGER upper bound")
    U.previous <- miscPack:::locate.previous(x, U)
    U.next <- miscPack:::locate.next(x, U)
    if(sum(x >=U.previous & x <= U.next)>0)
      refined.interval <- range(y[sum(x >=U.previous & x <= U.next)>0])
    else  stop("initialize lambda.interval using a BIGGER upper bound")
  }
}

objective <- function(lam, U, tail.prob,
                      K=K,
                      basemodel=basemodel,
                      alpha0=alpha0,
                      n.obs=n.obs, n.samples=n.samples,
                      agrid = agrid,
                      alphamax=alphamax){
  if (basemodel==0){
    (1 - miscPack:::cdf.impliedprior.Kplus(U, lam, K=K,
                                           basemodel = basemodel,
                                           alpha0=alpha0,
                                           n.obs=n.obs, n.samples=n.samples,
                                           agrid = agrid,
                                           alphamax=alphamax,
                                           TR=TRUE) - tail.prob)^2
  } else {
    (miscPack:::cdf.impliedprior.Kplus(U, lam, K=K,
                                       basemodel = basemodel,
                                       alpha0=alpha0,
                                       n.obs=n.obs, n.samples=n.samples,
                                       agrid = agrid,
                                       alphamax=alphamax,
                                       TR=TRUE) - tail.prob)^2
  }
}

cdf.impliedprior.Kplus <- function(U,
                                   lam,
                                   K=K,
                                   basemodel,
                                   alpha0=alpha0,
                                   n.obs=n.obs,
                                   n.samples=n.samples,
                                   agrid = agrid,
                                   alphamax=alphamax,
                                   TR=TRUE){
  sample.Kplus <- miscPack:::draw.Kplus(lam=lam,
                                        alpha0=alpha0,
                                        K=K,
                                        basemodel = basemodel,
                                        n.obs=n.obs, n.samples=n.samples,
                                        agrid = agrid,
                                        alphamax=alphamax,
                                        TR=TRUE)
  tab.Kplus <- table(sample.Kplus$Kplus)/n.samples
  sum(tab.Kplus[as.numeric(names(tab.Kplus))<=U])
}


lambda_finder <- function(U, tail.prob,
                          alpha0,
                          agrid,
                          K,
                          basemodel=0,
                          n.obs,
                          n.samples,
                          alphamax,
                          truncation,
                          lambda.interval = c(0.0001, 0.2),
                          presearch.step=FALSE){
  if(presearch.step){
    # pre-search
    print("Start pre-search step   ")
    ini.lamba.lwr <- lambda.interval[1]
    ini.lamba.upr <- lambda.interval[2]
    lambda.grid <- seq(ini.lamba.lwr,
                       ini.lamba.upr,
                       (ini.lamba.upr-ini.lamba.lwr)/20)
    quantile.grid <- rep(NA, length(lambda.grid))
    for (i in 1:length(lambda.grid)){
      sample.Kplus <- miscPack:::draw.Kplus(lam=lambda.grid[i],
                                            K,
                                            basemodel=basemodel,
                                            alpha0=alpha0,
                                            n.obs=n.obs,
                                            n.samples=n.samples,
                                            agrid=agrid,
                                            alphamax=alphamax,
                                            TR=truncation)

      if (basemodel==1) {
        quantile.grid[i] <- quantile(sample.Kplus$Kplus, tail.prob)
      } else {
        quantile.grid[i] <- quantile(sample.Kplus$Kplus, 1-tail.prob)
      }
    }

    print("Quantiles are: ")
    print(quantile.grid)

    lambda.ini.range <- miscPack:::refine.lambda.int(lambda.grid,
                                                     round(quantile.grid, 0),
                                                     U, K,
                                                     basemodel = basemodel)
  } else {
    lambda.ini.range <- lambda.interval
  }

  result <- optimize(f=miscPack:::objective,
                     interval=lambda.ini.range,
                     maximum = FALSE,
                     U=U,
                     K=K,
                     tail.prob=tail.prob,
                     basemodel = basemodel,
                     alpha0=alpha0,
                     n.obs=n.obs, n.samples=n.samples,
                     agrid = agrid,
                     alphamax=alphamax)
  result$minimum
}




################
#### select lambda such that Pr(Kplus > U)=0.99
#   OBSOLETE
################

#
#
# lambda_givenU <- function(U,
#                           tail.prob,
#                           K,
#                           basemodel,
#                           alpha0,
#                           n.obs,n.samples,agrid,
#                           alphamax=alphamax,
#                           interval=c(1e-3,2.5))
# {
#   uniroot(objective,interval=interval,
#           #extendInt = "yes",
#           U=U,          ### this is U such that Pr(Kplus>U)=tail.prob
#           tail.prob=tail.prob,  ### tail.prob can be either fixed by us, or we may let the user handle it
#           K=K, alpha0=alpha0,
#           n.obs=n.obs, n.samples=n.samples,
#           agrid = agrid,
#           alphamax=alphamax)$root
# }





#### #### #### #### #### Asymmetric Dir(a1,a2) ####

# EVALUATE pc prior density for a1 \in (0,U]
# with a2 fixed at alpha2.fixed=1-e5

pcprior.asym.a1 <- function(a.eval,
                            lam=1,
                            U=5,
                            K=25,
                            alpha1.base,
                            alpha2.base,
                            alpha2.fixed=alpha2.base,
                            agrid,
                            TR=T){

  a <- agrid
  n.agrid <- length(a)
  # alpha.vec.base <- c(rep(alpha1.base, U-1),
  #                     rep(alpha2.base, K-U+1))
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha1 <- a[i]
    #alpha <- c(rep(alpha1, U-1), rep(alpha2.fixed, K-U+1))
    alpha <- c(rep(alpha1, U), rep(alpha2.fixed, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)

  # expontential density (PC prior is based on exponential of KL distance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR){
    dmax <- b.texp
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(a, fexp(lambda = lam,
                                   d=a.to.d(a)) * abs(a.to.d(a, deriv=1)))

  tmp <- pcprior(a.eval)
  return(ifelse(tmp < 0, 0, tmp))
}


# EVALUATE pc prior density for log(a1) \in (-inf,0]
# with a2 fixed at alpha2.fixed=1-e5
pcprior.asym.log_a1 <- function(a.eval,
                                lam=1,
                                U=5,
                                K=25,
                                alpha1.base,
                                alpha2.base,
                                alpha2.fixed=alpha2.base,
                                agrid,
                                TR=T){

  a <- agrid
  n.agrid <- length(a)
  # alpha.vec.base <- c(rep(alpha1.base, U-1),
  #                     rep(alpha2.base, K-U+1))
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha1 <- a[i]
    #alpha <- c(rep(alpha1, U-1), rep(alpha2.fixed, K-U+1))
    alpha <- c(rep(alpha1, U), rep(alpha2.fixed, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR){
    dmax <- b.texp
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(a, fexp(lambda = lam,
                                   d=a.to.d(a)) * abs(a.to.d(a, deriv=1)))
  pcprior.spline <- splinefun(a, pcprior(a))
  pcprior.spline.logscale <- splinefun(log(a), pcprior.spline(a)*a)
  tmp <- pcprior.spline.logscale(log_a.eval)
  return(ifelse(tmp < 0, 0, tmp))
}


## find lambda
# (decay rate par for the pcprior on a1 with fixed a2, asymm dir)
objective.pcprior.asym.a1 <- function(lam,
                                      U,
                                      U.lwr = U-1,
                                      tail.prob,
                                      K,
                                      alpha1.base,
                                      alpha2.base,
                                      alpha2.fixed=alpha2.base,
                                      agrid,
                                      n.obs,
                                      n.samples,
                                      TR=TRUE){
    (miscPack:::cdf.impliedpcprior.asym.a1.Kplus(lam,
                                               U = U,
                                               U.lwr = U.lwr,
                                               K = K,
                                               alpha1.base=alpha1.base,
                                               alpha2.base=alpha2.base,
                                               alpha2.fixed=alpha2.fixed,
                                               agrid=agrid,
                                               n.obs=n.obs,
                                               n.samples=n.samples,
                                               TR=TRUE) - tail.prob)^2
  }

cdf.impliedpcprior.asym.a1.Kplus <- function(lam,
                                   U,
                                   U.lwr = U-1,
                                   K,
                                   alpha1.base,
                                   alpha2.base,
                                   alpha2.fixed=alpha2.base,
                                   agrid,
                                   n.obs,
                                   n.samples,
                                   TR=TRUE){
  a <- agrid
  n.agrid <- length(a)
  # alpha.vec.base <- c(rep(alpha1.base, U-1),
  #                     rep(alpha2.base, K-U+1))
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha1 <- a[i]
    #alpha <- c(rep(alpha1, U-1), rep(alpha2.fixed, K-U+1))
    alpha <- c(rep(alpha1, U), rep(alpha2.fixed, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)
  d.to.a <- splinefun(d_val, a)

  # exponential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR) {
    sample.d <- rtexp(n.samples, rate = lam, a=0, b=b.texp)
  } else {
    sample.d <- rexp(n.samples, rate = lam)
  }
  sample.a <- d.to.a(sample.d)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(a.eval, U),
                                    rep(alpha2.fixed, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  tab.Kplus <- table(Kplus)/n.samples
  sum(tab.Kplus[as.numeric(names(tab.Kplus))<=U.lwr])
}

lambda_finder_pcprior.asym.a1 <- function(lam,
                                          U,
                                          U.lwr = U-1,
                                          tail.prob,
                                          K,
                                          alpha1.base,
                                          alpha2.base,
                                          alpha2.fixed=alpha2.base,
                                          agrid,
                                          n.obs,
                                          n.samples,
                                          TR=TRUE,
                                          lambda.interval, ...){
  result <- optimize(f=miscPack:::objective.pcprior.asym.a1,
                     interval=lambda.interval,
                     maximum = FALSE,
                     U = U,
                     U.lwr = U.lwr,
                     tail.prob = tail.prob,
                     K = K,
                     alpha1.base=alpha1.base,
                     alpha2.base=alpha2.base,
                     alpha2.fixed=alpha2.fixed,
                     agrid=agrid,
                     n.obs=n.obs,
                     n.samples=n.samples,
                     TR=TR, ...)
  result$minimum
}


### this would be 'draw.Kplus' (old name used for the symm dirichlet case)
impliedpcprior.asym.a1.Kplus <- function(lam,
                                       U,
                                       K,
                                       alpha1.base,
                                       alpha2.base,
                                       alpha2.fixed=alpha2.base,
                                       agrid,
                                       n.obs,
                                       n.samples,
                                       TR=TRUE){
  a <- agrid
  n.agrid <- length(a)
  # alpha.vec.base <- c(rep(alpha1.base, U-1),
  #                     rep(alpha2.base, K-U+1))
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha1 <- a[i]
    #alpha <- c(rep(alpha1, U-1), rep(alpha2.fixed, K-U+1))
    alpha <- c(rep(alpha1, U), rep(alpha2.fixed, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)
  d.to.a <- splinefun(d_val, a)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR) {
    sample.d <- rtexp(n.samples, rate = lam, a=0, b=b.texp)
  } else {
    sample.d <- rexp(n.samples, rate = lam)
  }
  sample.a <- d.to.a(sample.d)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(a.eval, U),
                                    rep(alpha2.fixed, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus, a=sample.a)
}




#### #### #### #### #### Asymmetric Dir(a1,a2) ####

# EVALUATE pc prior density for a2 \in (0,1]
# with a1 fixed at alpha1.fixed=1

pcprior.asym.a2 <- function(a.eval,
                            lam=1,
                            U=5,
                            K=25,
                            alpha1.base,
                            alpha2.base,
                            alpha1.fixed=alpha1.base,
                            agrid,
                            TR=T){

  a <- agrid
  n.agrid <- length(a)
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha2 <- a[i]
    alpha <- c(rep(alpha1.fixed, U), rep(alpha2, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR){
    dmax <- b.texp
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(agrid, fexp(lambda = lam,
                                   d=a.to.d(agrid)) * abs(a.to.d(agrid, deriv=1)))

  tmp <- pcprior(a.eval)
  return(ifelse(tmp < 0, 0, tmp))
}


# EVALUATE pc prior density for log(a2) \in (-inf,0]
# with a1 fixed at alpha1.fixed=1
pcprior.asym.log_a2 <- function(a.eval,
                                lam=1,
                                U=5,
                                K=25,
                                alpha1.base,
                                alpha2.base,
                                alpha1.fixed=alpha1.base,
                                agrid,
                                TR=T){

  a <- agrid
  n.agrid <- length(a)
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha2 <- a[i]
    alpha <- c(rep(alpha1.fixed, U), rep(alpha2, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR){
    dmax <- b.texp
    fexp <- function(lambda, d) (lambda * exp(-lambda * d)/(1-exp(-lambda*dmax)))
  } else {
    fexp <- function(lambda, d) (lambda * exp(-lambda * d))
  }

  # This produces a spline representation of the density of "alpha"
  pcprior <- splinefun(a, fexp(lambda = lam,
                               d=a.to.d(a)) * abs(a.to.d(a, deriv=1)))
  pcprior.spline <- splinefun(a, pcprior(a))
  pcprior.spline.logscale <- splinefun(log(a), pcprior.spline(a)*a)
  tmp <- pcprior.spline.logscale(log_a.eval)
  return(ifelse(tmp < 0, 0, tmp))
}


## find lambda
# (decay rate par for the pcprior on a1 with fixed a2, asymm dir)
objective.pcprior.asym.a2 <- function(lam,
                                      U,
                                      U.upr = U+1,
                                      tail.prob,
                                      K,
                                      alpha1.base=1,
                                      alpha2.base=1e-5,
                                      alpha1.fixed=alpha1.base,
                                      agrid,
                                      n.obs,
                                      n.samples,
                                      TR=TRUE){
  (miscPack:::cdf.impliedpcprior.asym.a2.Kplus(lam,
                                             U = U,
                                             U.upr = U.upr,
                                             K = K,
                                             alpha1.base=alpha1.base,
                                             alpha2.base=alpha2.base,
                                             alpha1.fixed=alpha1.fixed,
                                             agrid=agrid,
                                             n.obs=n.obs,
                                             n.samples=n.samples,
                                             TR=TRUE) - tail.prob)^2
}

cdf.impliedpcprior.asym.a2.Kplus <- function(lam,
                                           U,
                                           U.upr = U+1,
                                           K,
                                           alpha1.base,
                                           alpha2.base,
                                           alpha1.fixed=alpha1.base,
                                           agrid,
                                           n.obs,
                                           n.samples,
                                           TR=TRUE){
  a <- agrid
  n.agrid <- length(a)
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha2 <- a[i]
    alpha <- c(rep(alpha1.fixed, U), rep(alpha2, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)
  d.to.a <- splinefun(d_val, a)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR) {
    sample.d <- rtexp(n.samples, rate = lam, a=0, b=b.texp)
  } else {
    sample.d <- rexp(n.samples, rate = lam)
  }
  sample.a <- d.to.a(sample.d)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(alpha1.fixed, U),
                                    rep(a.eval, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  tab.Kplus <- table(Kplus)/n.samples
  1 - sum(tab.Kplus[as.numeric(names(tab.Kplus)) <= U.upr])
}

lambda_finder_pcprior.asym.a2 <- function(lam,
                                          U,
                                          U.upr = U-1,
                                          tail.prob,
                                          K,
                                          alpha1.base=1,
                                          alpha2.base=1e-5,
                                          alpha1.fixed=alpha1.base,
                                          agrid,
                                          n.obs,
                                          n.samples,
                                          TR=TRUE,
                                          lambda.interval,
                                          ...){
  result <- optimize(f=miscPack:::objective.pcprior.asym.a2,
                     interval=lambda.interval,
                     maximum = FALSE,
                     U = U,
                     U.upr = U.upr,
                     tail.prob = tail.prob,
                     K = K,
                     alpha1.base=alpha1.base,
                     alpha2.base=alpha2.base,
                     alpha1.fixed=alpha1.fixed,
                     agrid=agrid,
                     n.obs=n.obs,
                     n.samples=n.samples,
                     TR=TR, ...)
  result$minimum
}


### this would be 'draw.Kplus' (old name used for the symm dirichlet case)
impliedpcprior.asym.a2.Kplus <- function(lam,
                                       U,
                                       K,
                                       alpha1.base,
                                       alpha2.base,
                                       alpha1.fixed=alpha1.base,
                                       agrid,
                                       n.obs,
                                       n.samples,
                                       TR=TRUE){
  a <- agrid
  n.agrid <- length(a)
  alpha.vec.base <- c(rep(alpha1.base, U),
                      rep(alpha2.base, K-U))

  # agrid: seq values for which kullback-Leibler distance will be evaluated
  kld_val <- d_val <- rep(NA, n.agrid)
  for(i in 1:n.agrid){
    alpha2 <- a[i]
    alpha <- c(rep(alpha1.fixed, U), rep(alpha2, K-U))
    alpha0 <- alpha.vec.base
    c1 <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - (
      lgamma(sum(alpha0)) - sum(lgamma(alpha0)))
    c2 <- sum((alpha - alpha0)*(digamma(alpha) - digamma(sum(alpha))))
    kld_val[i] <- c1 + c2
    d_val[i] <- sqrt(2*kld_val[i])
  }
  a.to.d <- splinefun(a, d_val)
  d.to.a <- splinefun(d_val, a)

  # expontential density (PC prior is based on exponential of KL dist(ance)
  ## compute the pc prior on a1 (not available in closed form)
  b.texp <- max(d_val)
  if(TR) {
    sample.d <- rtexp(n.samples, rate = lam, a=0, b=b.texp)
  } else {
    sample.d <- rexp(n.samples, rate = lam)
  }
  sample.a <- d.to.a(sample.d)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(alpha1.fixed, U),
                                    rep(a.eval, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus, a=sample.a)
}






#### find a1~Gamma(a,a) parameter,
# we want this centred at a1=1;
# find 'a' such that Pr(Kplus<=U.lwr) = b,
# with b small, e.g. 10% or 20%

objective.gammaprior.asym.a1 <- function(gamma.par,
                                         U,
                                         U.lwr = U-1,
                                         tail.prob,
                                         K,
                                         alpha2.fixed=1e-5,
                                         n.obs,
                                         n.samples,
                                         sc){
  (miscPack:::cdf.implied.gammaprior.asym.a1.Kplus(gamma.par,
                                                   U = U,
                                                   U.lwr = U.lwr,
                                                   K = K,
                                                   alpha2.fixed=alpha2.fixed,
                                                   n.obs=n.obs,
                                                   n.samples=n.samples,
                                                   sc=sc) - tail.prob)^2
}

cdf.implied.gammaprior.asym.a1.Kplus <- function(gamma.par,
                                                 U,
                                                 U.lwr = U-1,
                                                 K,
                                                 alpha2.fixed=1e-5,
                                                 n.obs,
                                                 n.samples,
                                                 sc){

  sample.a <- rgamma(n.samples,
                     shape = gamma.par,
                     rate = gamma.par*sc)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(a.eval, U),
                                    rep(alpha2.fixed, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  tab.Kplus <- table(Kplus)/n.samples
  sum(tab.Kplus[as.numeric(names(tab.Kplus))<=U.lwr])
}


gammapar_finder_gammaprior.asym.a1 <- function(gamma.par,
                                               U,
                                               U.lwr = U-1,
                                               tail.prob,
                                               K,
                                               alpha2.fixed=1e-5,
                                               n.obs,
                                               n.samples,
                                               sc,
                                               gamma.par.interval,
                                               ...){
  result <- optimize(f=miscPack:::objective.gammaprior.asym.a1,
                     interval=gamma.par.interval,
                     maximum = FALSE,
                     U = U,
                     U.lwr = U.lwr,
                     tail.prob = tail.prob,
                     K = K,
                     alpha2.fixed=alpha2.fixed,
                     n.obs=n.obs,
                     n.samples=n.samples,
                     sc= sc,
                     ...)
  result$minimum
}


### this would be 'draw.Kplus' (old name used for the symm dirichlet case)
implied.gammaprior.asym.a1.Kplus <- function(gamma.par,
                                             U,
                                             K,
                                             alpha2.fixed=1e-5,
                                             n.obs,
                                             n.samples,
                                             sc){

  sample.a <- rgamma(n.samples,
                     shape = gamma.par,
                     rate = gamma.par*sc)
  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    a.eval <- sample.a[i]
    sample.w[i,] <- rdirichlet(1, c(rep(a.eval, U),
                                    rep(alpha2.fixed, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
    #print(i)
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus, a=sample.a)
}





########## implied prior on K+ for an Asym Dir with fixed a1 and a2
### this would be 'draw.Kplus' (old name used for the symm dirichlet case)
implied.static.asym.Kplus <- function(U,
                                      K,
                                      alpha1.fixed=U,
                                      alpha2.fixed=1e-5,
                                      n.obs,
                                      n.samples){

  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    sample.w[i,] <- rdirichlet(1, c(rep(alpha1.fixed, U),
                                    rep(alpha2.fixed, K-U)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus)
}




########## implied prior on K+ for an Sym Dir with fixed a
### this would be 'draw.Kplus' (old name used for the symm dirichlet case)
implied.static.sym.Kplus <- function(K,
                                     alpha.fixed=1e-5,
                                     n.obs,
                                     n.samples){

  sample.z <- sample.w <- array(NA, c(n.samples, K))
  for(i in 1:n.samples){
    sample.w[i,] <- rdirichlet(1, c(rep(alpha.fixed, K)))
    if(is.na(sample.w[i,1])){
      sample.z[i,] <- c(1, rep(0, K-1))
    } else {
      sample.z[i,] <- rmultinom(n=1,
                                size=n.obs,
                                prob=sample.w[i,])
    }
  }
  Kplus <- K-apply(sample.z, 1, function(x) sum(x==0))
  list(Kplus=Kplus)
}


