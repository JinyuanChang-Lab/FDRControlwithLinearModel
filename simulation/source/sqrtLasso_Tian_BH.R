sqrtLasso_Tian_BH <- function(n, p, k, amplitude, alpha, rho, kappa = 0.8, estim_tol = 1e-3, bits = NULL) {
  library(glmnet)
  library(MASS)
  library(Rmpfr)
  
  if (!exists("TG.pvalue", mode = "function")) {
    stop(
      paste0(
      )
    )
  }
  
  if (!is.null(bits) && !requireNamespace("Rmpfr", quietly = TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits <- NULL
  }
  
  choose_lambda <- function(X, kappa, seed = 2024, draws = 10000) {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    
    set.seed(seed)
    
    n_obs <- nrow(X)
    G <- matrix(rnorm(draws * n_obs), nrow = draws, ncol = n_obs)
    ratios <- apply(G %*% X, 1, function(v) max(abs(v))) / sqrt(rowSums(G^2))
    as.numeric(kappa * unname(quantile(ratios, probs = 0.95)))
  }
  
  solve_sqrt_lasso <- function(X, y, lambda, maxit = 2000) {
    p <- ncol(X)
    f <- function(u) {
      beta <- u[1:p] - u[(p+1):(2*p)]
      r <- as.vector(y - X %*% beta)
      sqrt(sum(r^2)) + lambda * sum(u)
    }
    gr <- function(u) {
      beta <- u[1:p] - u[(p+1):(2*p)]
      r <- as.vector(y - X %*% beta)
      nr <- sqrt(sum(r^2))
      g <- if (nr < 1e-12) rep(0, p) else -as.vector(crossprod(X, r)) / nr
      c(g + lambda, -g + lambda)
    }
    opt <- optim(
      par = rep(0, 2 * p), fn = f, gr = gr,
      method = "L-BFGS-B", lower = rep(0, 2 * p),
      control = list(maxit = maxit, factr = 1e5, pgtol = 1e-10)
    )
    
    if (opt$convergence != 0) {
      stop(sprintf("square-root Lasso optimizer failed: %s", opt$message))
    }
    
    beta <- opt$par[1:p] - opt$par[(p+1):(2*p)]
    beta[abs(beta) < 2e-7] <- 0
    active <- which(beta != 0)
    signs <- sign(beta[active])
    
    r <- as.vector(y - X %*% beta)
    score <- as.vector(crossprod(X, r)) / sqrt(sum(r^2))
    inactive <- setdiff(seq_len(p), active)
    
    if (length(active) && max(abs(score[active] - lambda * signs)) > 1e-3) {
      stop("active square-root Lasso KKT check failed")
    }
    if (length(inactive) && max(abs(score[inactive]) - lambda) > 1e-3) {
      stop("inactive square-root Lasso KKT check failed")
    }
    
    list(beta = beta, active = active, signs = signs)
  }
  
  # Generate data
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
  
  nonzero <- sample(p, k)
  beta_true <- numeric(p)
  beta_true[nonzero] <- amplitude
  y <- as.vector(X %*% beta_true + rnorm(n))
  
  # Center data
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- scale(y, center = TRUE, scale = FALSE)
  
  # Choose lambda
  lambda_val <- choose_lambda(X, kappa = kappa)
  fit <- tryCatch({
    solve_sqrt_lasso(X, y, lambda = lambda_val)
  }, error = function(e) NULL)
  
  if (is.null(fit) || length(fit$active) == 0 || length(fit$active) >= n - 1) {
    selected <- integer(0)
  } else {
    active_idx <- fit$active
    s_sign <- fit$signs
    k_active <- length(active_idx)
    
    # Constraint  
    cons <- tryCatch({
      active_constraint(X, active_idx, s_sign, lambda_val)
    }, error = function(e) NULL)
    
    if (is.null(cons)) {
      selected <- integer(0)
    } else {
      df <- n - k_active - 1
      resid_ols <- as.vector(y - cons$P %*% y)
      sigma_hat <- sqrt(sum(resid_ols^2) / df)
      
      p_values <- rep(1, k_active)
      
      b_gaussian <- sigma_hat * sqrt(df) * cons$b
      Sigma_gaussian <- diag(rep(sigma_hat^2, n))
      
      for (j in seq_len(k_active)) {
        eta_j <- as.vector(cons$XE_dag[j, ])
        eta_norm <- sqrt(sum(eta_j^2))
        
        if (!is.finite(eta_norm) || eta_norm <= 0) {
          next
        }
        
        vj <- s_sign[j] * eta_j / eta_norm
        
        # p-values
        tg_out <- tryCatch({
          TG.pvalue(
            Z = y,
            A = cons$A,
            b = b_gaussian,
            eta = vj,
            Sigma = Sigma_gaussian,
            null_value = 0,
            bits = bits
          )
        }, error = function(e) NULL)
        
        if (!is.null(tg_out) &&
            length(tg_out$pv) == 1 &&
            is.finite(tg_out$pv)) {
          p_values[j] <- min(1, max(0, as.numeric(tg_out$pv)))
        }
      }
      
      #Near-estimability check
      XM <- X[, active_idx, drop = FALSE]
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
      valid_vars <- active_idx[keep_idx]
      valid_pv <- p_values[keep_idx]
      
      if (length(valid_pv) == 0) {
        selected <- integer(0)
      } else {
        # BH
        p_adj <- p.adjust(valid_pv, method = "BH")
        selected <- valid_vars[p_adj <= alpha]
      }
    }
  }
  
  # FDP and power
  fdp_val <- sum(beta_true[selected] == 0) / max(1, length(selected))
  power_val <- sum(beta_true[selected] != 0) / sum(beta_true != 0)
  
  rst <- cbind(fdp = fdp_val, power = power_val)
  return(rst)
}