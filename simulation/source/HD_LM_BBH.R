#' High-dimensional Bonferroni-Benjamini-Hochberg method
#'
#' @param X input design matrix of dimension n x p; each row is an observation vector.
#' @param y response variable.
#' @param alpha target FDR level.
#' @param nlambda number of lambda values in cross-validation for Lasso, default 100.
#' @param resol resolution parameter for InverseLinfty, default 1.5.
#' @param mu tuning parameter for InverseLinfty; if NULL, it is chosen automatically inside InverseLinfty.
#' @param maxiter maximum iterations for InverseLinfty, default 20.
#' @param threshold convergence threshold for InverseLinfty, default 1e-2.
#' @return algorithm1_selected: the index of selected variable of algorithm 1 .
#' @return algorithm2_selected: the index of selected variable of algorithm 2 .
#' @import scalreg glmnet knockoff

HD_LM_BBH = function(X, y, alpha,
                     nlambda = 100,
                     resol = 1.5,
                     mu = NULL,
                     maxiter = 20,
                     threshold = 1e-2){
  
  # Load packages
  library(scalreg);  
  library(glmnet);
  library(knockoff);
  
  n = nrow(X)
  p = ncol(X)
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  # Transformation matrix
  ATrans = rbind(diag(p),diag(p))
  BTrans = rbind(diag(p),-diag(p))
  Trans = cbind(ATrans, BTrans)
  
  # Run the knockoff filter to achieve model-X knockoff
  second_order_MX = function(X) create.second_order(X, method='asdp')
  result = knockoff.filter(X, y, knockoffs = second_order_MX, fdr = alpha)
  KX = result$Xk
  
  # Augmented linear model and Lasso
  AX = cbind(X,KX)
  fit = cv.glmnet(AX, y, nlambda = nlambda, alpha = 1, standardize = TRUE)
  Ab = as.vector(coef(fit, s = fit$lambda.1se))[-1]
  
  # Estimate the precision matrix
  G = t(AX) %*% AX / n
  Inv = InverseLinfty(G, n, resol = resol, mu = mu,
                      maxiter = maxiter, threshold = threshold)
  
  # Debiased Lasso
  Adb = as.numeric(Ab + (Inv %*% t(AX) %*% (y - AX %*% Ab)) / n)
  w2 = Adb[1:p] - Adb[(p+1):(2*p)]
  w1 = Adb[1:p] + Adb[(p+1):(2*p)]
  
  # Scaled Lasso
  Scaledlasso = scalreg(AX, y)
  hsigma = Scaledlasso$hsigma
  
  # Paired test statistics
  Cov = Trans %*% Inv %*% G %*% t(Inv) %*% Trans
  t1 = sqrt(n) * w1 / (hsigma * sqrt(diag(Cov[1:p, 1:p])))
  t2 = sqrt(n) * w2 / (hsigma * sqrt(diag(Cov[(p+1):(2*p), (p+1):(2*p)])))
  
  # Obtain paired p-values
  p1 = 2 * (1 - pnorm(abs(t1)))
  p2 = 2 * (1 - pnorm(abs(t2)))
  p21 = p2
  
  # Knockoff-assisted BH
  p2a = p.adjust(p2, method = "BH")
  index2 = seq(along = p2a)[p2a <= alpha]
  
  algorithm1_selected = index2
  
  # Knockoff-assisted BBH
  p21[p1 > sqrt(alpha)] = 1
  p21a = p.adjust(p21, method = "BH")
  index21 = seq(along = p21a)[p21a <= sqrt(alpha)]
  
  algorithm2_selected = index21
  
  # Bind the result
  return(results = list(algorithm1_selected = algorithm1_selected,
                        algorithm2_selected = algorithm2_selected))
}



InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=20, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=20, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}
