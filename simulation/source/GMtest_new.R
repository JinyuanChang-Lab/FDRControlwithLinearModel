GMtest_new = function(n,p,k,amplitude,alpha,type){
  
  # Load packages
  library(doParallel);
  library(foreach);
  library(glmnet);
  
  # Generate covariates from a multivariate normal distribution
  mu = rep(0,p)
  
  # Choose the dependence structure: Independent/AR1/Block
  if (type == 1){
    Sigma = toeplitz(0^(0:(p-1)))
  } else if (type == 2){
    Sigma = toeplitz(0.4^(0:(p-1)))
  } else {
    Sigma = matrix(data = 0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p){
        Sigma[i,j] = 0.2*(ceiling(i/20)==ceiling(j/20))
      }
    }
    diag(Sigma) = 1
  }
  
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Generate responses from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  # Calculate fdp
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  # Calculate power
  power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)
  
  # Fit GM
  result = gm_new(y, X, q = alpha, ncores=1)
  
  fdp_G = fdp(result$gm_selected)
  pwr_G = power(result$gm_selected)
  
  rst = cbind(fdp_G, pwr_G)
  return (rst)
}

