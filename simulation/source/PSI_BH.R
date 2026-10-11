PSI_BH <- function(n, p, k, amplitude, alpha, rho, estim_tol = 1e-3) {
  library(selectiveInference)
  library(glmnet)
  library(MASS)
  
  # Generate design matrix
  mu <- rep(0, p)
  if (rho == 0.2) {
    Sigma <- matrix(data = 0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p){
        Sigma[i,j] <- rho*(ceiling(i/20)==ceiling(j/20))
      }
    }
    diag(Sigma) <- 1
  } else {
    Sigma <- toeplitz(rho^(0:(p-1)))
  }
  X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
  
  # Generate true coefficients and response
  nonzero <- sample(p, k)
  beta_true <- amplitude * (1:p %in% nonzero)
  y <- X %*% beta_true + rnorm(n)
  
  # Center data
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- scale(y, center = TRUE, scale = FALSE)
  
  # Fit Lasso via cross-validation
  cvfit <- cv.glmnet(X, y, nlambda = 100, alpha = 1, standardize = FALSE)
  lambda_glmnet <- cvfit$lambda.1se
  lambda_psi <- lambda_glmnet * n
  
  betahat <- as.vector(coef(cvfit, x = X, y = y, s = lambda_glmnet, exact = TRUE))[-1]
  
  # Use true noise standard deviation
  hsigma <- 1
  
  # Run selective inference with error handling
  has_constraint_warning <- FALSE
  
  out <- tryCatch({
    result <- withCallingHandlers(
      expr = fixedLassoInf(X, y, betahat, lambda_psi, sigma = hsigma),
      warning = function(w) {
        if (grepl("Constraint not satisfied", w$message, fixed = TRUE)) {
          has_constraint_warning <<- TRUE
        }
      }
    )
    if (has_constraint_warning) NULL else result
  }, error = function(e) {
    NULL
  })
  
  if (is.null(out)) {
    selected_idx <- integer(0)
  } else {
    vars_active <- out$vars
    pv_active <- out$pv
    k_active <- length(vars_active)
    
    if (k_active == 0) {
      selected_idx <- integer(0)
    } else {
      # Near-estimability check
      XM <- X[, vars_active, drop = FALSE]
      XtX <- t(XM) %*% XM
      XtX_inv <- ginv(XtX)
      
      residual_vec <- numeric(k_active)
      for (j in 1:k_active) {
        ej <- rep(0, k_active)
        ej[j] <- 1
        cj <- XM %*% XtX_inv %*% ej
        residual_vec[j] <- sqrt(sum((t(XM) %*% cj - ej)^2))
      }
      
      keep_idx <- which(residual_vec <= estim_tol)
      valid_vars <- vars_active[keep_idx]
      valid_pv <- pv_active[keep_idx]
      
      if (length(valid_pv) == 0) {
        selected_idx <- integer(0)
      } else {
        # BH correction
        p_adj <- p.adjust(valid_pv, method = "BH")
        selected_idx <- valid_vars[p_adj <= alpha]
      }
    }
  }
  
  # Calculate FDP and power
  fdp <- function(selected) {
    sum(beta_true[selected] == 0) / max(1, length(selected))
  }
  
  power <- function(selected) {
    n_true <- sum(beta_true != 0)
    if (n_true == 0) return(0)
    sum(beta_true[selected] != 0) / n_true
  }
  
  fdp_PSI <- fdp(selected_idx)
  pwr_PSI <- power(selected_idx)
  
  rst <- cbind(fdp_PSI, pwr_PSI)
  return(rst)
}