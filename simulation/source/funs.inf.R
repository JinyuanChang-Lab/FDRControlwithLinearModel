
# Assuming that grid is in sorted order from smallest to largest,
# and vals are monotonically increasing function values over the
# grid, returns the grid end points such that the corresponding
# vals are approximately equal to {val1, val2}

grid.search <- function(grid, fun, val1, val2, gridpts=100, griddepth=2) {
  n = length(grid)
  vals = fun(grid)
    
  ii = which(vals >= val1)
  jj = which(vals <= val2)
  if (length(ii)==0) return(c(grid[n],Inf))   # All vals < val1
  if (length(jj)==0) return(c(-Inf,grid[1]))  # All vals > val2
  # RJT: the above logic is correct ... but for simplicity, instead,
  # we could just return c(-Inf,Inf) 

  i1 = min(ii); i2 = max(jj)
  if (i1==1) lo = -Inf
  else lo = grid.bsearch(grid[i1-1],grid[i1],fun,val1,gridpts,
         griddepth-1,below=TRUE)
  if (i2==n) hi = Inf
  else hi = grid.bsearch(grid[i2],grid[i2+1],fun,val2,gridpts,
         griddepth-1,below=FALSE)
  return(c(lo,hi))
}

# Repeated bin search to find the point x in the interval [left, right]
# that satisfies f(x) approx equal to val. If below=TRUE, then we seek
# x such that the above holds and f(x) <= val; else we seek f(x) >= val.

grid.bsearch <- function(left, right, fun, val, gridpts=100, griddepth=1, below=TRUE) {
  n = gridpts
  depth = 1

  while (depth <= griddepth) {
    grid = seq(left,right,length=n)
    vals = fun(grid)
    
    if (below) {
      ii = which(vals >= val)
      if (length(ii)==0) return(grid[n])   # All vals < val (shouldn't happen)
      if ((i0=min(ii))==1) return(grid[1]) # All vals > val (shouldn't happen)
      left = grid[i0-1]
      right = grid[i0]
    }
    
    else {
      ii = which(vals <= val)
      if (length(ii)==0) return(grid[1])   # All vals > val (shouldn't happen)
      if ((i0=max(ii))==n) return(grid[n]) # All vals < val (shouldn't happen)
      left = grid[i0]
      right = grid[i0+1]
    }

    depth = depth+1
  }

  return(ifelse(below, left, right))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

active_constraint <- function(X, active, signs, lambda) {
  XE <- X[, active, drop = FALSE]
  Ginv <- solve(crossprod(XE))
  XE_dag <- Ginv %*% t(XE)
  z <- as.numeric(signs)
  wE <- as.vector(t(XE_dag) %*% (lambda * z))
  denom <- 1 - sum(wE * wE)
  
  if (denom <= 0) {
    stop("invalid square-root Lasso KKT denominator")
  }
  
  D <- diag(z, nrow = length(z), ncol = length(z))
  A <- -D %*% XE_dag
  b <- -z * as.vector(Ginv %*% (lambda * z)) / sqrt(denom)
  
  P <- XE %*% XE_dag
  
  list(A = A, b = b, P = P, XE_dag = XE_dag)
}

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 

  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

# Returns Prob(Z>z | Z in [a,b]), where mean cane be a vector, using
# multi precision floating point calculations thanks to the Rmpfr package

mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                      (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)

  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)

  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
         (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

# Return Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# Riemann approximation tricks, by Max G'Sell

gsell.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  return(max.approx.frac(a/sd,b/sd,z/sd,mean/sd))
}


estimate_sigma <- function(observed,
                           truncated_df,
                           lower_bound,
                           upper_bound,
                           untruncated_df = 0,
                           factor = 3,
                           npts = 50,
                           nsample = 2000) {
  
  if (!is.finite(observed) || observed <= 0) {
    stop("observed must be a positive finite residual standard deviation")
  }
  if (!is.finite(truncated_df) || truncated_df <= 0) {
    stop("truncated_df must be positive")
  }
  if (!is.finite(lower_bound) || lower_bound < 0) {
    stop("lower_bound must be finite and nonnegative")
  }
  if (is.finite(upper_bound) && upper_bound <= lower_bound) {
    stop("upper_bound must be greater than lower_bound")
  }
  if (factor <= 1 || npts < 2 || nsample < 1) {
    stop("invalid factor, npts, or nsample in estimate_sigma")
  }
  
  # Same regularization rule as the Python implementation.
  if (untruncated_df < 50) {
    linear_term <- sqrt(truncated_df)
  } else {
    linear_term <- 0
  }
  
  total_df <- untruncated_df + truncated_df
  
  # Candidate sigma values: [observed / factor, observed * factor].
  values <- seq(
    from = 1 / factor,
    to = factor,
    length.out = npts
  ) * observed
  
  expected <- numeric(npts)
  
  # For each candidate sigma, approximate the conditional expected RSS
  # by sampling from the truncated chi distribution, exactly as in Python.
  for (i in seq_along(values)) {
    value <- values[i]
    
    # Python uses scipy.stats.chi.cdf(). Since chi_d^2 ~ chisq_d,
    # F_chi(a) = F_chisq(a^2).
    P_upper <- pchisq(
      (upper_bound * sqrt(truncated_df) / value)^2,
      df = truncated_df
    )
    
    P_lower <- pchisq(
      (lower_bound * sqrt(truncated_df) / value)^2,
      df = truncated_df
    )
    
    U <- runif(nsample)
    prob <- (P_upper - P_lower) * U + P_lower
    
    # Python: chidist.ppf(prob, df)^2
    truncated_chisq <- qchisq(prob, df = truncated_df)
    
    if (untruncated_df > 0) {
      # Python: chidist.rvs(untruncated_df)^2
      untruncated_chisq <- rchisq(nsample, df = untruncated_df)
      sample_value <- (
        truncated_chisq + untruncated_chisq
      ) * value^2
    } else {
      sample_value <- truncated_chisq * value^2
    }
    
    expected[i] <- mean(sample_value)
    
    if (expected[i] >= 1.5 * (
      observed^2 * total_df + observed^2 * linear_term
    )) {
      break
    }
  }
  
  # Python interp1d(values, expected + values**2 * linear_term).
  interpolant_y <- expected + values^2 * linear_term
  
  V <- seq(
    from = 1 / factor,
    to = factor,
    length.out = 10 * npts
  ) * observed
  
  interp_V <- approx(
    x = values,
    y = interpolant_y,
    xout = V,
    rule = 1
  )$y
  
  # This is the Python RPL root-search criterion. With
  # untruncated_df = 0 and d = truncated_df, dividing by d gives
  #
  #   H_d(L,U,sigma^2) + sigma^2 / sqrt(d)
  #     = (1 + 1/sqrt(d)) * observed^2,
  #
  # i.e. Tian et al.'s regularized pseudolikelihood equation (21).
  target <- observed^2 * total_df + observed^2 * linear_term
  
  candidate <- V[
    !is.na(interp_V) & interp_V >= target
  ]
  
  # Same fallback as the Python try/except block.
  if (length(candidate) == 0) {
    sigma_hat <- observed
  } else {
    sigma_hat <- min(candidate)
  }
  
  sigma_hat
}

# ---------------------------------------------------------------------------
# Build the RPL truncation interval from active_constraint(), then call the
# Python-translated estimate_sigma() above.
#
# active_constraint() uses
#
#   A y <= ||(I-P)y||_2 * b.
#
# Therefore r = ||(I-P)y||_2 is truncated to [lower_r, upper_r]. The Python
# estimate_sigma() expects bounds on sigma_P = r / sqrt(nu), so we pass
#
#   lower_bound = lower_r / sqrt(nu),
#   upper_bound = upper_r / sqrt(nu),
#   observed    = r_obs   / sqrt(nu).
# ---------------------------------------------------------------------------
estimate_sigma2_rpl <- function(y, cons, nu,
                                constraint_tol = 1e-8,
                                factor = 3,
                                npts = 50,
                                nsample = 2000) {
  if (nu <= 0) {
    stop("RPL variance estimator requires positive degrees of freedom")
  }
  
  resid <- as.vector(y - cons$P %*% y)
  rss_obs <- sum(resid^2)
  
  if (!is.finite(rss_obs) || rss_obs <= 0) {
    stop("invalid residual sum of squares for RPL variance estimation")
  }
  
  r_obs <- sqrt(rss_obs)
  observed <- r_obs / sqrt(nu)
  sigma2_ols <- observed^2
  
  Ay <- as.vector(cons$A %*% y)
  b0 <- as.vector(cons$b)
  
  if (length(Ay) != length(b0)) {
    stop("incompatible dimensions in RPL selection constraints")
  }
  
  scale_tol <- constraint_tol * max(1, max(abs(Ay)), r_obs)
  
  pos <- which(b0 > constraint_tol)
  neg <- which(b0 < -constraint_tol)
  zer <- which(abs(b0) <= constraint_tol)
  
  # Constraints with b_i = 0 must already satisfy (Ay)_i <= 0.
  if (length(zer) && any(Ay[zer] > scale_tol)) {
    stop("observed data do not satisfy zero-slope RPL selection constraints")
  }
  
  lower_r <- 0
  if (length(pos)) {
    lower_r <- max(0, max(Ay[pos] / b0[pos]))
  }
  
  upper_r <- Inf
  if (length(neg)) {
    upper_r <- min(Ay[neg] / b0[neg])
  }
  
  if (!is.finite(lower_r) || lower_r < 0) {
    stop("invalid lower truncation boundary for RPL variance estimation")
  }
  
  if (is.finite(upper_r) && upper_r <= 0) {
    stop("invalid upper truncation boundary for RPL variance estimation")
  }
  
  if (is.finite(upper_r) && lower_r >= upper_r) {
    stop("empty truncation interval for RPL variance estimation")
  }
  
  # The observed residual norm must lie in its own selection interval.
  bound_tol <- constraint_tol * max(
    1,
    r_obs,
    lower_r,
    if (is.finite(upper_r)) upper_r else 0
  )
  
  if (r_obs < lower_r - bound_tol ||
      (is.finite(upper_r) && r_obs > upper_r + bound_tol)) {
    stop("observed residual norm is outside the RPL truncation interval")
  }
  
  # Convert residual-norm bounds to the scale expected by estimate_sigma().
  lower_bound <- lower_r / sqrt(nu)
  upper_bound <- if (is.infinite(upper_r)) Inf else upper_r / sqrt(nu)
  
  sigma_hat <- estimate_sigma(
    observed = observed,
    truncated_df = nu,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    untruncated_df = 0,
    factor = factor,
    npts = npts,
    nsample = nsample
  )
  
  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    stop("invalid RPL sigma estimate")
  }
  
  list(
    sigma2 = sigma_hat^2,
    sigma = sigma_hat,
    sigma2_ols = sigma2_ols,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    lower_rss = lower_r^2,
    upper_rss = if (is.infinite(upper_r)) Inf else upper_r^2
  )
}

##############################

forwardStop <- function(pv, alpha=.10){
  if (alpha<0 || alpha>1) stop("alpha must be in [0,1]")
  if (min(pv,na.rm=T)<0 || max(pv,na.rm=T)>1) stop("pvalues must be in [0,1]")
  val=-(1/(1:length(pv)))*cumsum(log(1-pv))
  oo = which(val <= alpha)
  if (length(oo)==0) out=0
  else out = oo[length(oo)]
  return(out)
}

##############################

aicStop <- function(x, y, action, df, sigma, mult=2, ntimes=2) {
  n = length(y)
  k = length(action)
  aic = numeric(k)
  G = matrix(0,nrow=0,ncol=n)
  u = numeric(0)
  count = 0
  
  for (i in 1:k) {
    A = action[1:i]
    aic[i] = sum(lsfit(x[,A],y,intercept=F)$res^2) + mult*sigma^2*df[i]

    j = action[i]
    if (i==1) xtil = x[,j]
    else xtil = lsfit(x[,action[1:(i-1)]],x[,j],intercept=F)$res
    s = sign(sum(xtil*y))
    
    if (i==1 || aic[i] <= aic[i-1]) {
      G = rbind(G,s*xtil/sqrt(sum(xtil^2)))
      u = c(u,sqrt(mult)*sigma)
      count = 0
    }

    else {
      G = rbind(G,-s*xtil/sqrt(sum(xtil^2)))
      u = c(u,-sqrt(mult)*sigma)
      count = count+1
      if (count == ntimes) break
    }
  }

  if (i < k) {
    khat = i - ntimes
    aic = aic[1:i]
  }
  else khat = k
  
  return(list(khat=khat,G=G,u=u,aic=aic,stopped=(i<k)))
}

#these next two functions are used by the binomial and Cox options of fixedLassoInf

# Compute the truncation interval and SD of the corresponding Gaussian

TG.limits = function(Z, A, b, eta, Sigma=NULL) {

    target_estimate = sum(as.numeric(eta) * as.numeric(Z))

    if (max(A %*% as.numeric(Z) - b) > 0) {
        warning('Constraint not satisfied. A %*% Z should be elementwise less than or equal to b')
    }

    if (is.null(Sigma)) {
        Sigma = diag(rep(1, n))
    }

    # compute pvalues from poly lemma:  full version from Lee et al for full matrix Sigma

    n = length(Z)
    eta = matrix(eta, ncol=1, nrow=n)
    b = as.vector(b)
    var_estimate = sum(matrix(eta, nrow=1, ncol=n) %*% (Sigma %*% matrix(eta, ncol=1, nrow=n)))
    cross_cov = Sigma %*% matrix(eta, ncol=1, nrow=n)
   
    resid = (diag(n) - matrix(cross_cov / var_estimate, ncol=1, nrow=n) %*% matrix(eta, nrow=1, ncol=n)) %*% Z
    rho = A %*% cross_cov / var_estimate
    vec = (b - as.numeric(A %*% resid)) / rho

    vlo = suppressWarnings(max(vec[rho < 0]))
    vup = suppressWarnings(min(vec[rho > 0]))

    sd = sqrt(var_estimate)
    return(list(vlo=vlo, vup=vup, sd=sd, estimate=target_estimate))
}

TG.pvalue = function(Z, A, b, eta, Sigma=NULL, null_value=0, bits=NULL) {

    limits.info = TG.limits(Z, A, b, eta, Sigma)

    return(TG.pvalue.base(limits.info, null_value=null_value, bits=bits))
}

TG.interval = function(Z, A, b, eta, Sigma=NULL, alpha=0.1, 
                       gridrange=c(-100,100),
                       gridpts=100, 
                       griddepth=2, 
                       flip=FALSE, 
                       bits=NULL) {

    limits.info = TG.limits(Z, A, b, eta, Sigma)

    return(TG.interval.base(limits.info, 
                            alpha=alpha, 
                            gridrange=gridrange,
                            griddepth=griddepth,
			    flip=flip,
			    bits=bits))
}

TG.interval.base = function(limits.info, alpha=0.1, 
                            gridrange=c(-100,100),
                            gridpts=100, 
                            griddepth=2, 
                            flip=FALSE, 
                            bits=NULL) {

    # compute sel intervals from poly lemmma, full version from Lee et al for full matrix Sigma

    param_grid = seq(gridrange[1] * limits.info$sd, gridrange[2] * limits.info$sd, length=gridpts)

    pivot = function(param) {
        pv = tnorm.surv(limits.info$estimate, param, limits.info$sd, limits.info$vlo, limits.info$vup, bits) 
        return(pv)
    }

    interval = grid.search(param_grid, pivot, alpha/2, 1-alpha/2, gridpts, griddepth)
    tailarea = c(pivot(interval[1]), 1- pivot(interval[2]))

    if (flip) {
        interval = -interval[2:1]
        tailarea = tailarea[2:1]
     }
 
     # int is not a good variable name, synonymous with integer...
     return(list(int=interval,
                 tailarea=tailarea))
}

TG.pvalue.base = function(limits.info, null_value=0, bits=NULL) {
    pv = tnorm.surv(limits.info$estimate, null_value, limits.info$sd, limits.info$vlo, limits.info$vup, bits)
    return(list(pv=pv, vlo=limits.info$vlo, vup=limits.info$vup, sd=limits.info$sd))
}


mydiag=function(x){
    if(length(x)==1) out=x
    if(length(x)>1) out=diag(x)
       return(out)
   }

