BBH = function(n, p, k, amplitude, alpha, rho, safety = 0.99) {
  stopifnot(
    n > 2*p, p >= 1, k >= 0, k <= p,
    is.finite(amplitude), is.finite(rho), abs(rho) < 1,
    alpha > 0, alpha < 1,
    length(safety) == 1L, is.finite(safety), safety > 0, safety < 1
  )
  
  # Generate data
  Sigma = toeplitz(rho^(0:(p - 1)))
  X = matrix(rnorm(n*p), n, p) %*% chol(Sigma)
  beta = numeric(p)
  beta[sample.int(p, k)] = amplitude
  y = as.vector(X %*% beta + rnorm(n))
  
  # Normalize column lengths
  column_norms = sqrt(colSums(X^2))
  Xs = sweep(X, 2, column_norms, "/")
  X.svd = svd(Xs)
  if (any(X.svd$d <= 1e-5 * max(X.svd$d)))
    stop("X must have full column rank for BBH")
  
  # The complement needs to be orthogonal to X
  Q = qr.Q(qr(cbind(X.svd$u, matrix(0, n, p))))
  U_perp = Q[, p + seq_len(p), drop = FALSE]
  
  s = safety * min(2 * min(X.svd$d)^2, 1)
  s_diff = pmax(0, 2*s - (s/X.svd$d)^2)
  KX = (sweep(X.svd$u, 2, X.svd$d - s/X.svd$d, "*") +
          sweep(U_perp, 2, sqrt(s_diff), "*")) %*% t(X.svd$v)
  KX = sweep(KX, 2, column_norms, "*")
  AX = cbind(X, KX)
  
  # Sum and difference estimators
  G = crossprod(X)
  GD = crossprod(X, KX)
  IGMD = solve(G + GD)
  ID = solve(G - GD)
  w1 = as.vector(IGMD %*% crossprod(X + KX, y))
  w2 = as.vector(ID %*% crossprod(X - KX, y))
  
  # OLS
  fit_ols = lm.fit(AX, y)
  if (fit_ols$rank != 2*p)
    stop("The augmented design is numerically rank deficient")
  nu = n - 2*p
  hsigma = sqrt(sum(fit_ols$residuals^2) / nu)
  
  t1 = w1 / (hsigma * sqrt(2 * diag(IGMD)))
  t2 = w2 / (hsigma * sqrt(2 * diag(ID)))
  p1 = 2 * pt(-abs(t1), df = nu)
  p2 = 2 * pt(-abs(t2), df = nu)
  
  # Bonferroni-BH
  p21 = p2
  p21[p1 >= sqrt(alpha)] = 1
  selected = which(p.adjust(p21, method = "BH") <= sqrt(alpha))
  
  fdp_BBH = sum(beta[selected] == 0) / max(1, length(selected))
  n_signal = sum(beta != 0)
  pwr_BBH = if (n_signal > 0) sum(beta[selected] != 0) / n_signal else NA_real_
  cbind(fdp_BBH = fdp_BBH, pwr_BBH = pwr_BBH)
}
