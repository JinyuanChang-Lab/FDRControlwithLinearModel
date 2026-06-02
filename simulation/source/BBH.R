BBH = function(n,p,k,amplitude,alpha,rho){
  
  # Load packages
  library(knockoff);
  
  # Generate the variables from a multivariate normal distribution
  mu = rep(0,p)
  Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Transformation matrix
  ATrans = rbind(diag(p),diag(p))
  BTrans = rbind(diag(p),-diag(p))
  Trans = cbind(ATrans, BTrans)
  
  # Generate the response from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n)
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  # Generate knockoffs
  KX = create.fixed(X, method = "asdp")$Xk
  AX = cbind(X, KX)
  
  # Covariance matrix
  G = t(X)%*%X
  GD = t(X)%*%KX
  
  # 2Sigma-D
  IGMD = solve(G+GD)
  # D
  ID = solve(G-GD)
  
  # Estimation
  w1 = IGMD%*%t(X+KX)%*%y
  w2 = ID%*%t(X-KX)%*%y
  
  # hsigma
  ############################################
  lm.fit = lm(y ~ AX - 1)
  hsigma = sigma(lm.fit)
  
  # test statistics
  t1 = w1/(hsigma*sqrt(2)*sqrt(diag(IGMD)))
  t2 = w2/(hsigma*sqrt(2)*sqrt(diag(ID)))
  
  # p-value
  p1 = 2*(1 - pt(abs(t1), df = n - 2*p))
  p2 = 2*(1 - pt(abs(t2), df = n - 2*p))
  p21 = p2
  
  # adjusted p-value
  p2a = p.adjust(p2, method = "BH")
  
  # Knockoff-assisted BH
  index2 = seq(along = p2a)[p2a <= alpha]
  
  # Knockoff-assisted BBH
  p21[p1 >= sqrt(alpha)] = 1
  p21a = p.adjust(p21, method = "BH")
  index21 = seq(along = p21a)[p21a <= alpha]
  
  # Calculate fdp
  fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  # Calculate power
  power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)
  
  fdp_BBH = fdp(index21)
  pwr_BBH = power(index21)
  
  rst = cbind(fdp_BBH, pwr_BBH)
  return (rst)
}
