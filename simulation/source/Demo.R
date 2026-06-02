source("./HD_LM_BBH.R")

# Parameters
n = 200
p = 400
k = 16
amplitude = 0.5
rho = 0.4

# Generate covariates from a multivariate normal distribution
mu = rep(0,p)
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate responses from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)

# Run HD_LM_BBH
Results = HD_LM_BBH(X, y, alpha = 0.1)

# Calculate fdp
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))

# Calculate power
power = function(selected) sum(beta[selected] != 0) / sum(beta != 0)

# Check fdp and power
fdp(Results$algorithm1_selected)
fdp(Results$algorithm2_selected)
power(Results$algorithm1_selected)
power(Results$algorithm2_selected)