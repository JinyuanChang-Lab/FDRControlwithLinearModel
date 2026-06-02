DSBBH = function(n,p,k,amplitude,alpha,rho){
  
  library(scalreg);  
  library(glmnet);
  library(knockoff);
  library(expm);
                                          
  # Generate covariates from a multivariate normal distribution
  mu = rep(0,p)
  Sigma = toeplitz(rho^(0:(p-1)))                                                  
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Generate responses from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  # Transformation matrix
  ATrans = rbind(diag(p),diag(p))
  BTrans = rbind(diag(p),-diag(p))
  Trans = cbind(ATrans, BTrans)
  
  # Data splitting
  half <- n %/% 2
  n1 = half
  n2 = n - n1
  y1 = y[1:n1] 
  X1 = X[1:n1,] 
  X2 = X[(n1+1):n,]
  
  # Centering X and y
  X1 = scale(X1, center = TRUE, scale = FALSE)
  X2 = scale(X2, center = TRUE, scale = FALSE)
  y1 = scale(y1, center = TRUE, scale = FALSE)
  
  # Generate knockoff
  Sigma2 = t(X2)%*%X2/n2
  Omega = InverseLinfty(Sigma2, n2)
  lambdamax = max(Re(eigen(Omega)$values))
  D = 1.7/lambdamax*diag(p)
  xi = matrix(rnorm(n1 * p), nrow = n1, ncol = p)
  KX = X1%*%(diag(p) - Omega%*%D) + xi%*%sqrtm(2*D - D%*%Omega%*%D)
  
  # Augmented linear model and Lasso
  AX = cbind(X1,KX)
  fit = cv.glmnet(AX, y1, nlambda=100, alpha=1, standardize=T)
  Ab = as.vector(coef(fit, s=fit$lambda.1se))[-1]
  
  # Estimate the precision matrix
  G = t(AX)%*%AX/n1
  Inv = InverseLinfty(G, n1)
  
  # Debiased Lasso
  Adb = as.numeric(Ab + (Inv%*%t(AX)%*%(y1 - AX %*% Ab))/n1)
  w2 = Adb[1:p] - Adb[(p+1):(2*p)]
  w1 = Adb[1:p] + Adb[(p+1):(2*p)]
  
  # Scaled Lasso
  Scaledlasso = scalreg(AX, y1)
  hsigma = Scaledlasso$hsigma
  
  # Paired test statistics
  Cov = Trans%*%Inv%*%G%*%t(Inv)%*%Trans
  t1 = sqrt(n1)*w1/(hsigma*sqrt(diag(Cov[1:p, 1:p])))
  t2 = sqrt(n1)*w2/(hsigma*sqrt(diag(Cov[(p+1):(2*p), (p+1):(2*p)])))
  
  # Obtain paired p-values
  p1 = 2*(1 - pnorm(abs(t1)))
  p2 = 2*(1 - pnorm(abs(t2)))
  p21 = p2
  
  # Knockoff-assisted BH
  p2a = p.adjust(p2, method = "BH")
  index2 = seq(along = p2a)[p2a <= alpha]
  
  # Knockoff-assisted BBH
  p21[p1 > sqrt(alpha)] = 1
  p21a = p.adjust(p21, method = "BH")
  index21 = seq(along = p21a)[p21a <= sqrt(alpha)]
  
  # Calculate fdp
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  # Calculate power
  power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)
  
  fdp_KBH = fdp(index2) # fdp of Algorithm 1
  pwr_KBH = power(index2) # power of Algorithm 1
  fdp_KBBH = fdp(index21) # fdp of Algorithm 2
  pwr_KBBH = power(index21) # power of Algorithm 2
  
  rst = cbind(fdp_KBH, pwr_KBH, fdp_KBBH, pwr_KBBH)
  return (rst)
}

