OBBH = function(n,p,k,amplitude,alpha,rho){
  
  # Load packages
  library(scalreg);  
  library(glmnet);
  
  # Generate covariates from a multivariate normal distribution
  mu = rep(0,p)
  Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Generate responses from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  # Debiased Lasso
  fit = cv.glmnet(X, y, nlambda=100, alpha=1, standardize=T)
  b = as.vector(coef(fit, s=fit$lambda.1se))[-1]
  
  G = t(X)%*%X/n
  Inv = InverseLinfty(G, n)
  
  db = as.numeric(b + (Inv%*%t(X)%*%(y - X %*% b))/n)
  
  # Scaled Lasso
  Scaledlasso = scalreg(X, y)
  hsigma = Scaledlasso$hsigma
  
  # Test statistic
  Cov = Inv%*%G%*%t(Inv)
  t = sqrt(n)*db/(hsigma*sqrt(diag(Cov[1:p, 1:p])))
  
  # p-value
  p = 2*(1 - pnorm(abs(t)))
  pp = p
  
  # adjusted p-value
  p1 = p.adjust(p, method = "BH")
  
  # BH
  index1 = seq(along = p1)[p1 <= alpha]
  
  # BBH
  pp[p > sqrt(alpha)] = 1
  p2 = p.adjust(pp, method = "BH")
  index2 = seq(along = p2)[p2 <= sqrt(alpha)]
  
  # Calculate fdp
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  # Calculate power
  power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)
  
  fdp_OBH = fdp(index1)
  pwr_OBH = power(index1)
  fdp_OBBH = fdp(index2)
  pwr_OBBH = power(index2)
  
  rst = cbind(fdp_OBH, pwr_OBH, fdp_OBBH, pwr_OBBH)
  return (rst)
}

