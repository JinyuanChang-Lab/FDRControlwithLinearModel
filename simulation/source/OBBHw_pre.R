OBBHw_pre = function(n,p,k,amplitude,alpha,rho){
  
  library(scalreg);  
  library(glmnet);
  library(knockoff);
  
  # Generate the variables from a multivariate normal distribution
  mu = rep(0,p)
  Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Generate the response from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  ########################################
  #pre-selection with lasso
  
  lambda_max = 1*max(abs(t(X) %*% y))/(n^(3/2))
  lambda_min = lambda_max / 1e2
  nlambda=100
  kk = (0:(nlambda-1)) / nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^kk
  fit0.lasso = cv.glmnet(X, y, family = "gaussian", lambda=lambda, nfold=10,alpha=1)
  lambda_ast = fit0.lasso$lambda.min
  
  coef_ast = coef(fit0.lasso, s='lambda.min')[-1]
  if(sum(coef_ast!=0)>0.9*n){
    ord_idx = order(abs(coef_ast[which(coef_ast!=0)]))
    coef_ast[which(coef_ast!=0)[1:(sum(coef_ast!=0)-0.9*n)]]=0
  }
  
  idx =  intersect(which(coef_ast!=0), c(1:p))
  
  ######################################## 
  #Run original procedure on the subindex
  sX = X[ , idx]
  sp = ncol(sX)
  
  fit = cv.glmnet(sX, y, nlambda=100, alpha=1, standardize=T)
  b0 = as.vector(coef(fit, s=fit$lambda.1se))[-1]
  
  G = t(sX)%*%sX/n
  Inv = InverseLinfty(G, n)
  
  db = as.numeric(b0 + (Inv%*%t(sX)%*%(y - sX %*% b0))/n)
  
  # Scaled lasso
  Scaledlasso = scalreg(sX, y)
  hsigma = Scaledlasso$hsigma
  
  # Test statistic
  Cov = Inv%*%G%*%t(Inv)
  t = sqrt(n)*db/(hsigma*sqrt(diag(Cov[1:sp, 1:sp])))
  
  # p-value
  p = 2*(1 - pnorm(abs(t)))
  
  pp = p
  
  # adjusted p-value
  p1 = p.adjust(p, method = "BH")
  
  # BH
  index1 = seq(along = p1)[p1 <= alpha]
  
  # Trans into original index
  index1 = idx[index1]
  
  # BBH
  pp[p > sqrt(alpha)] = 1
  p2 = p.adjust(pp, method = "BH")
  index2 = seq(along = p2)[p2 <= sqrt(alpha)]
  
  # Trans into original index
  index2 = idx[index2]
  
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

