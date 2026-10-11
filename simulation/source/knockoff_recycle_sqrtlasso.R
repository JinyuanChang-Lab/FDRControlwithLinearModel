knockoff_recycle_sqrtlasso <- function(n, p, k, amplitude, alpha, rho,
                                       n0 = floor(n * 0.4),
                                       kmax_ratio = 0.25,
                                       kappa = 0.9) {
  library(glmnet)
  library(knockoff)
  library(MASS)
  
  # -------------------------------------------------------------------------
  # Baseline lambda: quantile of ||X^T e||_inf / ||e||_2,
  # e ~ N(0, I), normalized so that ||e||_2 = 1
  # -------------------------------------------------------------------------
  lambda_baseline <- function(X, quantile = 0.95, ndraw = 10000) {
    n_obs <- nrow(X)
    E <- matrix(rnorm(n_obs * ndraw), n_obs, ndraw)
    E <- scale(E, center = FALSE, scale = sqrt(colSums(E^2)))
    XtE_abs <- abs(crossprod(X, E))
    max_vals <- apply(XtE_abs, 2, max)
    return(as.numeric(quantile(max_vals, probs = quantile)))
  }
  
  # -------------------------------------------------------------------------
  # Sign-restricted square-root Lasso
  # -------------------------------------------------------------------------
  solve_sign_restricted_sqrt_lasso <- function(X, y, lambda, sign_vec,
                                               maxit = 2000,
                                               kkt_tol = 1e-3) {
    q <- ncol(X)
    if (length(sign_vec) != q) {
      stop("length(sign_vec) must equal ncol(X)")
    }
    if (any(!sign_vec %in% c(-1, 1))) {
      stop("sign_vec must contain only -1 or +1")
    }
    
    Z <- sweep(X, 2, sign_vec, FUN = "*")
    
    f <- function(theta) {
      r <- as.vector(y - Z %*% theta)
      sqrt(sum(r^2)) + lambda * sum(theta)
    }
    
    gr <- function(theta) {
      r <- as.vector(y - Z %*% theta)
      nr <- sqrt(sum(r^2))
      g <- if (nr < 1e-12) {
        rep(0, q)
      } else {
        -as.vector(crossprod(Z, r)) / nr
      }
      g + lambda
    }
    
    opt <- optim(
      par = rep(0, q),
      fn = f,
      gr = gr,
      method = "L-BFGS-B",
      lower = rep(0, q),
      control = list(maxit = maxit, factr = 1e5, pgtol = 1e-10)
    )
    
    if (opt$convergence != 0) {
      stop(sprintf(
        "sign-restricted square-root Lasso optimizer failed: %s",
        opt$message
      ))
    }
    
    theta <- opt$par
    theta[abs(theta) < 2e-7] <- 0
    beta <- sign_vec * theta
    
    r_kkt <- as.vector(y - Z %*% theta)
    nr_kkt <- sqrt(sum(r_kkt^2))
    
    if (!is.finite(nr_kkt) || nr_kkt <= 0) {
      stop("invalid residual norm in sign-restricted square-root Lasso KKT check")
    }
    
    score_kkt <- as.vector(crossprod(Z, r_kkt)) / nr_kkt
    active_theta <- which(theta > 0)
    inactive_theta <- which(theta == 0)
    
    active_kkt_error <- if (length(active_theta)) {
      max(abs(score_kkt[active_theta] - lambda))
    } else {
      0
    }
    
    inactive_kkt_excess <- if (length(inactive_theta)) {
      max(0, max(score_kkt[inactive_theta] - lambda))
    } else {
      0
    }
    
    if (active_kkt_error > kkt_tol) {
      stop(sprintf(
        "active sign-restricted square-root Lasso KKT check failed: max error = %.6g",
        active_kkt_error
      ))
    }
    
    if (inactive_kkt_excess > kkt_tol) {
      stop(sprintf(
        "inactive sign-restricted square-root Lasso KKT check failed: max excess = %.6g",
        inactive_kkt_excess
      ))
    }
    
    return(beta)
  }
  
  # ========== 1. Generate design and response ==========
  if (rho == 0.2) {
    Sigma <- matrix(data = 0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p){
        Sigma[i,j] <- rho*(ceiling(i/20)==ceiling(j/20))
      }
    }
    diag(Sigma) <- 1
  } else {
    Sigma <- toeplitz(rho^(0:(p - 1)))
  }
  X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
  
  nonzero <- sample(p, k)
  beta_true <- amplitude * (1:p %in% nonzero)
  y <- as.vector(X %*% beta_true + rnorm(n))
  
  # ========== 2. Split data; first part is used only for screening ==========
  idx_screen <- sample(n, n0)
  idx_remaining <- setdiff(seq_len(n), idx_screen)
  
  X0 <- X[idx_screen, , drop = FALSE]
  y0 <- y[idx_screen]
  X1 <- X[idx_remaining, , drop = FALSE]
  y1 <- y[idx_remaining]
  n1 <- length(idx_remaining)
  
  fit_path <- glmnet(
    X0,
    y0,
    alpha = 1,
    standardize = FALSE,
    intercept = TRUE,
    nlambda = 200
  )
  beta_path <- as.matrix(coef(fit_path))[-1, , drop = FALSE]
  
  entry_lambda <- rep(0, p)
  entry_sign <- rep(0, p)
  for (j in seq_len(p)) {
    nonzero_idx <- which(abs(beta_path[j, ]) > 1e-8)
    if (length(nonzero_idx) > 0) {
      first_idx <- min(nonzero_idx)
      entry_lambda[j] <- fit_path$lambda[first_idx]
      entry_sign[j] <- sign(beta_path[j, first_idx])
    }
  }
  
  kmax <- min(floor(n1 * kmax_ratio), floor(n1 / 2) - 1)
  active_idx <- order(entry_lambda, decreasing = TRUE)[seq_len(kmax)]
  active_idx <- active_idx[entry_lambda[active_idx] > 0]
  k_active <- length(active_idx)
  
  if (k_active == 0 || k_active >= n1 / 2) {
    selected <- integer(0)
  } else {
    s_sign <- entry_sign[active_idx]
    X0_active <- X0[, active_idx, drop = FALSE]
    X1_active <- X1[, active_idx, drop = FALSE]
    
    # ========== 3. Construct knockoffs ONLY on the remaining data ==========
    Xk1 <- tryCatch({
      res <- suppressWarnings(create.fixed(X1_active, method = "equi"))
      sweep(res$Xk, 2, sqrt(colSums(X1_active^2)), FUN = "*")
    }, warning = function(w) NULL, error = function(e) NULL)
    
    if (is.null(Xk1)) {
      selected <- integer(0)
    } else {
      # ========== 4. Data recycling construction ==========
      X_recycle <- rbind(X0_active, X1_active)
      Xk_recycle <- rbind(X0_active, Xk1)
      y_recycle <- c(y0, y1)
      
      # ========== 5. Sign-restricted square-root Lasso on recycled data ==========
      X_aug <- cbind(X_recycle, Xk_recycle)
      base_lam <- lambda_baseline(X_aug)
      lambda_val <- kappa * base_lam
      
      sign_aug <- c(s_sign, s_sign)
      
      beta_hat <- tryCatch({
        solve_sign_restricted_sqrt_lasso(
          X = X_aug,
          y = y_recycle,
          lambda = lambda_val,
          sign_vec = sign_aug
        )
      }, error = function(e) NULL)
      
      if (is.null(beta_hat)) {
        selected <- integer(0)
      } else {
        beta_orig <- beta_hat[seq_len(k_active)]
        beta_k <- beta_hat[k_active + seq_len(k_active)]
        W <- abs(beta_orig) - abs(beta_k)
        
        T_val <- knockoff.threshold(W, fdr = alpha, offset = 1)
        selected <- active_idx[W >= T_val]
      }
    }
  }
  
  # ========== 6. Evaluate FDP and power ==========
  fdp_val <- sum(beta_true[selected] == 0) / max(1, length(selected))
  power_denom <- sum(beta_true != 0)
  power_val <- if (power_denom > 0) {
    sum(beta_true[selected] != 0) / power_denom
  } else {
    0
  }
  
  rst <- cbind(fdp = fdp_val, power = power_val)
  return(rst)
}