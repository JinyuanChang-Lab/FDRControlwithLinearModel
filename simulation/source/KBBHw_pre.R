KBBHw_pre = function(n,p,k,amplitude,alpha,rho){
  
  # Load packages
  library(scalreg);  
  library(glmnet);
  library(knockoff);
  
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
  
  ########################################
  # Pre-selection with Lasso
  
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
  #Run debiased Lasso on the subindex
  sX = X[ , idx]
  sp = ncol(sX)
  
  # Run the knockoff filter to achieve MX knockoff X
  KX = create.second_order(X, method='asdp')
  sKX = KX[ , idx]
  
  # Augmented linear model and debiased Lasso
  AX = cbind(sX, sKX)
  
  AX0 = cbind(X, KX)
  fit = cv.glmnet(AX0, y, nlambda=100, alpha=1, standardize=T)
  Ab0 = as.vector(coef(fit, s=fit$lambda.1se))[-1]
  
  kidx = idx + p
  aidx = c(idx, kidx)
  
  Ab = Ab0[aidx]
  
  G = t(AX)%*%AX/n
  Inv = InverseLinfty(G, n)
  
  Adb = as.numeric(Ab + (Inv%*%t(AX)%*%(y - AX %*% Ab))/n)
  w2 = Adb[1:sp] - Adb[(sp+1):(2*sp)]
  w1 = Adb[1:sp] + Adb[(sp+1):(2*sp)]
  
  # Scaled Lasso
  Scaledlasso = scalreg(AX, y)
  hsigma = Scaledlasso$hsigma
  
  # Transformation matrix
  ATrans = rbind(diag(sp),diag(sp))
  BTrans = rbind(diag(sp),-diag(sp))
  Trans = cbind(ATrans, BTrans)
  
  # Test statistic
  Cov = Trans%*%Inv%*%G%*%t(Inv)%*%Trans
  t1 = sqrt(n)*w1/(hsigma*sqrt(diag(Cov[1:sp, 1:sp])))
  t2 = sqrt(n)*w2/(hsigma*sqrt(diag(Cov[(sp+1):(2*sp), (sp+1):(2*sp)])))
  
  # p-value
  p1 = 2*(1 - pnorm(abs(t1)))
  p2 = 2*(1 - pnorm(abs(t2)))
  p21 = p2
  
  # adjusted p-value
  p2a = p.adjust(p2, method = "BH")
  
  # Knockoff-assisted BH
  index2 = seq(along = p2a)[p2a <= alpha]
  
  # Trans into original index
  index2 = idx[index2]
  
  # Knockoff-assisted BBH
  p21[p1 >= sqrt(alpha)] = 1
  p21a = p.adjust(p21, method = "BH")
  index21 = seq(along = p21a)[p21a <= sqrt(alpha)]
  
  # Trans into original index
  index21 = idx[index21]
  
  # Calculate fdp
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  # Calculate power
  power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)
  
  fdp_KBH = fdp(index2)
  pwr_KBH = power(index2)
  fdp_KBBH = fdp(index21)
  pwr_KBBH = power(index21)
  
  rst = cbind(fdp_KBH, pwr_KBH, fdp_KBBH, pwr_KBBH)
  return (rst)
}

