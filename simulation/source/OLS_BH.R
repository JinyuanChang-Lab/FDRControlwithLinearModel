OLS_BH <- function(n, p, k, amplitude, alpha, rho, 
                   n0 = floor(n * 0.4),
                   kmax_ratio = 0.25,
                   estim_tol = 1e-3) {
  library(glmnet)
  library(MASS)
  
  # Generate design matrix
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
  y <- as.vector(X %*% beta_true + rnorm(n))
  
  # Split data into screening set and inference set
  idx_screen <- sample(n, n0)
  X0 <- X[idx_screen, , drop = FALSE]
  y0 <- y[idx_screen]
  X1 <- X[-idx_screen, , drop = FALSE]
  y1 <- y[-idx_screen]
  n1 <- n - n0
  
  # Center screening and inference data separately after splitting
  X0 <- scale(X0, center = TRUE, scale = FALSE)
  y0 <- as.vector(scale(y0, center = TRUE, scale = FALSE))
  X1 <- scale(X1, center = TRUE, scale = FALSE)
  y1 <- as.vector(scale(y1, center = TRUE, scale = FALSE))
  
  # ========== Screening: Lasso path, top kmax by entry order ==========
  fit_path <- glmnet(X0, y0, alpha = 1, standardize = FALSE, nlambda = 200)
  beta_path <- as.matrix(coef(fit_path))[-1, ]
  
  # Record entry lambda and sign at first entry for each feature
  entry_lambda <- rep(0, p)
  entry_sign <- rep(0, p)
  for (j in 1:p) {
    nonzero_idx <- which(abs(beta_path[j, ]) > 1e-8)
    if (length(nonzero_idx) > 0) {
      first_idx <- min(nonzero_idx)
      entry_lambda[j] <- fit_path$lambda[first_idx]
      entry_sign[j] <- sign(beta_path[j, first_idx])
    }
  }
  
  # Select top kmax features by entry order
  kmax <- min(floor(n1 * kmax_ratio), floor(n1 / 2) - 1)
  active_idx <- order(entry_lambda, decreasing = TRUE)[1:kmax]
  active_idx <- active_idx[entry_lambda[active_idx] > 0]
  k_active <- length(active_idx)
  
  if (k_active == 0 || k_active >= n1) {
    selected_idx <- integer(0)
  } else {
    s_sign <- entry_sign[active_idx]
    
    # ========== Inference stage: OLS estimation on inference set ==========
    XM <- X1[, active_idx, drop = FALSE]
    XtX <- t(XM) %*% XM
    XtX_inv <- tryCatch({
      solve(XtX)
    }, error = function(e) {
      ginv(XtX)
    })
    
    beta_ls <- as.vector(XtX_inv %*% t(XM) %*% y1)
    resid <- as.vector(y1 - XM %*% beta_ls)
    df <- n1 - qr(XM)$rank - 1
    sigma_hat <- sqrt(sum(resid^2)/df)
    
    # Calculate t-statistics
    se <- sigma_hat * sqrt(diag(XtX_inv))
    t_vec <- beta_ls / se
    
    # Calculate one-sided p-values 
    pv_active <- numeric(k_active)
    for (j in 1:k_active) {
      if (s_sign[j] == 1) {
        pv_active[j] <- pt(t_vec[j], df = df, lower.tail = FALSE)
      } else {
        pv_active[j] <- pt(t_vec[j], df = df, lower.tail = TRUE)
      }
    }
    
    # Near-estimability check
    residual_vec <- numeric(k_active)
    for (j in 1:k_active) {
      ej <- rep(0, k_active)
      ej[j] <- 1
      cj <- XM %*% XtX_inv %*% ej
      residual_vec[j] <- sqrt(sum((t(XM) %*% cj - ej)^2))
    }
    
    keep_idx <- which(residual_vec <= estim_tol)
    valid_vars <- active_idx[keep_idx]
    valid_pv <- pv_active[keep_idx]
    
    if (length(valid_pv) == 0) {
      selected_idx <- integer(0)
    } else {
      # BH correction
      p_adj <- p.adjust(valid_pv, method = "BH")
      selected_idx <- valid_vars[p_adj <= alpha]
    }
  }
  
  # ========== Evaluate FDP and power ==========
  fdp <- function(selected) {
    sum(beta_true[selected] == 0) / max(1, length(selected))
  }
  
  power <- function(selected) {
    n_true <- sum(beta_true != 0)
    if (n_true == 0) return(0)
    sum(beta_true[selected] != 0) / n_true
  }
  
  fdp_OLS <- fdp(selected_idx)
  pwr_OLS <- power(selected_idx)
  
  rst <- cbind(fdp_OLS, pwr_OLS)
  return(rst)
}