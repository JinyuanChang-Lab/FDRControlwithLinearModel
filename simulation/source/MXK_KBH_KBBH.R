K_BH_BBH = function(n,p,k,amplitude,alpha,rho){

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
fit = cv.glmnet(AX, y, nlambda=100, alpha=1, standardize=T)
Ab = as.vector(coef(fit, s=fit$lambda.1se))[-1]

# Estimate the precision matrix
G = t(AX)%*%AX/n
Inv = InverseLinfty(G, n)

# Debiased Lasso
Adb = as.numeric(Ab + (Inv%*%t(AX)%*%(y - AX %*% Ab))/n)
w2 = Adb[1:p] - Adb[(p+1):(2*p)]
w1 = Adb[1:p] + Adb[(p+1):(2*p)]

# Scaled Lasso
Scaledlasso = scalreg(AX, y)
hsigma = Scaledlasso$hsigma

# Paired test statistics
Cov = Trans%*%Inv%*%G%*%t(Inv)%*%Trans
t1 = sqrt(n)*w1/(hsigma*sqrt(diag(Cov[1:p, 1:p])))
t2 = sqrt(n)*w2/(hsigma*sqrt(diag(Cov[(p+1):(2*p), (p+1):(2*p)])))

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
fdp_K = fdp(result$selected) # fdp of model-X knockoff
pwr_K = power(result$selected) # power of model-X knockoff

# Bind the result
rst = cbind(fdp_KBH, pwr_KBH, fdp_KBBH, pwr_KBBH, fdp_K, pwr_K)
return (rst)
}

