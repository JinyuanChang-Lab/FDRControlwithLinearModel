library(ggplot2)

env_ar1           <- new.env()
env_indep         <- new.env()
env_tian          <- new.env()
env_knock         <- new.env()
env_knock_recycle <- new.env()
env_ols           <- new.env()
env_psi           <- new.env()

load("./SI_plot_S2/S2_simulation.RData", envir = env_ar1)
load("./SI_plot_S2/S2_simulation_indep.RData", envir = env_indep)
load("./SI_plot_S2/S2_simulation_Tian_BH.RData", envir = env_tian)
load("./SI_plot_S2/S2_simulation_knockoff_split.RData", envir = env_knock)
load("./SI_plot_S2/S2_simulation_knockoff_recycle.RData", envir = env_knock_recycle)
load("./SI_plot_S2/S2_simulation_OLS_BH.RData", envir = env_ols)
load("./SI_plot_S2/S2_simulation_PSI_BH.RData", envir = env_psi)

load("./SI_plot_S2/S2_simulation_005.RData", envir = env_ar1)
load("./SI_plot_S2/S2_simulation_indep_005.RData", envir = env_indep)
load("./SI_plot_S2/S2_simulation_Tian_BH_005.RData", envir = env_tian)
load("./SI_plot_S2/S2_simulation_knockoff_split_005.RData", envir = env_knock)
load("./SI_plot_S2/S2_simulation_knockoff_recycle_005.RData", envir = env_knock_recycle)
load("./SI_plot_S2/S2_simulation_OLS_BH_005.RData", envir = env_ols)
load("./SI_plot_S2/S2_simulation_PSI_BH_005.RData", envir = env_psi)

replace_na_with_0 <- function(env) {
  for (var_name in ls(envir = env)) {
    var <- get(var_name, envir = env)
    if (is.matrix(var) || is.data.frame(var)) {
      var[is.na(var)] <- 0
      assign(var_name, var, envir = env)
    }
  }
}

replace_na_with_0(env_ar1)
replace_na_with_0(env_indep)
replace_na_with_0(env_tian)
replace_na_with_0(env_knock)
replace_na_with_0(env_knock_recycle)
replace_na_with_0(env_ols)
replace_na_with_0(env_psi)

get_plot_matrix <- function(metric, structure, n, p, amp_seq, fdr_level) {
  suffix <- if (fdr_level == 0.05) "_A005" else ""
  
  mat <- matrix(0, nrow = length(amp_seq), ncol = 7)
  colnames(mat) <- c("Alg1", "Alg2", "Knockoff", "Knockoff_Recycle", "OLS_BH", "PSI_BH", "Tian_BH")
  
  env_user <- if (structure == "ar1") env_ar1 else env_indep
  
  for (i in seq_along(amp_seq)) {
    amp <- amp_seq[i]
    
    r_var             <- paste0("r_n", n, "_p", p, "_amp", amp, suffix) 
    knock_var         <- paste0("r_knock_", structure, "_n", n, "_p", p, "_amp", amp, suffix)
    knock_recycle_var <- paste0("r_knock_recycle_", structure, "_n", n, "_p", p, "_amp", amp, suffix)
    ols_var           <- paste0("r_ols_", structure, "_n", n, "_p", p, "_amp", amp, suffix)
    psi_var           <- paste0("r_psi_", structure, "_n", n, "_p", p, "_amp", amp, suffix)
    tian_var          <- paste0("r_tian_", structure, "_n", n, "_p", p, "_amp", amp, suffix)
    
    idx_alg1 <- if (metric == "FDR") 1 else 2
    idx_alg2 <- if (metric == "FDR") 3 else 4
    idx_ext  <- if (metric == "FDR") 1 else 2
    
    if (exists(r_var, envir = env_user)) {
      r_val <- get(r_var, envir = env_user)
      mat[i, 1] <- colMeans(r_val)[idx_alg1]
      mat[i, 2] <- colMeans(r_val)[idx_alg2]
    }
    
    if (exists(knock_var, envir = env_knock)) {
      knock_val <- get(knock_var, envir = env_knock)
      mat[i, 3] <- colMeans(knock_val)[idx_ext]
    }
    
    if (exists(knock_recycle_var, envir = env_knock_recycle)) {
      knock_rec_val <- get(knock_recycle_var, envir = env_knock_recycle)
      mat[i, 4] <- colMeans(knock_rec_val)[idx_ext]
    }
    
    if (exists(ols_var, envir = env_ols)) {
      ols_val <- get(ols_var, envir = env_ols)
      mat[i, 5] <- colMeans(ols_val)[idx_ext] 
    }
    
    if (exists(psi_var, envir = env_psi)) {
      psi_val <- get(psi_var, envir = env_psi)
      mat[i, 6] <- colMeans(psi_val)[idx_ext]
    }
    
    if (exists(tian_var, envir = env_tian)) {
      tian_val <- get(tian_var, envir = env_tian)
      mat[i, 7] <- colMeans(tian_val)[idx_ext]
    }
  }
  return(mat)
}

amp_seq <- seq(0.1, 0.5, by = 0.05)
method_levels <- c("Alg1", "Alg2", "Knockoff", "Knockoff_Recycle", "OLS_BH", "PSI_BH", "Tian_BH")

points <- c(
  "Alg1"             = 15,
  "Alg2"             = 16,
  "Knockoff"         = 3, 
  "Knockoff_Recycle" = 8, 
  "OLS_BH"           = 17,
  "PSI_BH"           = 18,
  "Tian_BH"          = 25 
)

line_types <- c(
  "Alg1"             = "solid",
  "Alg2"             = "longdash",
  "Knockoff"         = "twodash",
  "Knockoff_Recycle" = "dotdash",
  "OLS_BH"           = "dashed",
  "PSI_BH"           = "dotted",
  "Tian_BH"          = "F1"
)

plot_metrics <- function(mat, amp_seq, title, save_name, metric) {
  type_vec <- rep(method_levels, each = length(amp_seq))
  x_vals   <- rep(amp_seq, times = length(method_levels))
  y_vals   <- as.vector(mat)
  df <- data.frame(
    x = x_vals, 
    Method = factor(type_vec, levels = method_levels), 
    value = y_vals
  )
  
  if (metric == "FDR") {
    y_limits <- c(0, 0.5)
    y_breaks <- seq(0, 0.5, by = 0.1)
    y_labels <- c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50")
    y_name   <- "Simulated FDR"
  } else {
    y_limits <- c(0, 1)
    y_breaks <- seq(0, 1, by = 0.2)
    y_labels <- c("0.00", "0.20", "0.40", "0.60", "0.80", "1.00")
    y_name   <- "Simulated Power"
  }
  
  p <- ggplot(data = df, mapping = aes(x = x, y = value, colour = Method, shape = Method, linetype = Method)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_x_continuous(name = "Amplitude", breaks = amp_seq, labels = as.character(amp_seq)) +
    scale_y_continuous(name = y_name, breaks = y_breaks, labels = y_labels, limits = y_limits) +
    ggtitle(title) +
    theme(
      legend.position = 'none', 
      plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
      axis.title = element_text(size = 26),
      axis.text = element_text(size = 22)
    )
  
  ggsave(paste0(save_name, ".pdf"), p, device = "pdf", path = "./SI_plot_S2", width = 8, height = 6)
  return(p)
}

n_values <- c(200, 500)
structures <- c("indep", "ar1")
metrics <- c("FDR", "Power")
fdr_levels <- c(0.1, 0.05)

for (n in n_values) {
  p_values <- if (n == 200) c(200, 300, 400) else c(500, 750, 1000)
  
  for (p in p_values) {
    for (struct in structures) {
      for (metric in metrics) {
        for (fdr in fdr_levels) {
          
          title <- paste0("n = ", n, ", d = ", p)
          
          suffix_name <- if (fdr == 0.05) "_005" else ""
          save_name <- paste0(struct, "_n", n, "_p", p, "_", tolower(metric), suffix_name)
          
          mat <- get_plot_matrix(metric, struct, n, p, amp_seq, fdr)
          plot_metrics(mat, amp_seq, title, save_name, metric)
        }
      }
    }
  }
}