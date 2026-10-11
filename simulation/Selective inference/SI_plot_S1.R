library(ggplot2)

env_sim  <- new.env()
env_ols  <- new.env()
env_psi  <- new.env()
env_tian <- new.env()

load("./SI_plot_S1/S1_simulation.RData", envir = env_sim)
load("./SI_plot_S1/S1_simulation_OLS_BH.RData", envir = env_ols)
load("./SI_plot_S1/S1_simulation_PSI_BH.RData", envir = env_psi)
load("./SI_plot_S1/S1_simulation_Tian_BH.RData", envir = env_tian)

load("./SI_plot_S1/S1_simulation_005.RData", envir = env_sim)
load("./SI_plot_S1/S1_simulation_OLS_BH_005.RData", envir = env_ols)
load("./SI_plot_S1/S1_simulation_PSI_BH_005.RData", envir = env_psi)
load("./SI_plot_S1/S1_simulation_Tian_BH_005.RData", envir = env_tian)

replace_na_with_0 <- function(env) {
  for (var_name in ls(envir = env)) {
    var <- get(var_name, envir = env)
    if (is.matrix(var) || is.data.frame(var)) {
      var[is.na(var)] <- 0
      assign(var_name, var, envir = env)
    }
  }
}

replace_na_with_0(env_sim)
replace_na_with_0(env_ols)
replace_na_with_0(env_psi)
replace_na_with_0(env_tian)

get_fdr_matrix <- function(fdr_level, structure, n, p_seq) {
  suffix <- switch(as.character(fdr_level),
                   "0.05" = "_A005",
                   "0.1"  = "",
                   "0.2"  = "_A02")
  
  struct_prefix <- switch(structure,
                          ind   = "0",
                          AR1   = "",
                          block = "E")
  
  mat <- matrix(0, nrow = length(p_seq), ncol = 5)
  colnames(mat) <- c("Alg1", "Alg2", "OLS_BH", "PSI_BH", "Tian_BH")
  
  for (i in seq_along(p_seq)) {
    p <- p_seq[i]
    
    r_var    <- paste0("r_", struct_prefix, "n", n, "_p", p, suffix)
    ols_var  <- paste0("r_ols_", struct_prefix, "n", n, "_p", p, suffix)
    psi_var  <- paste0("r_psi_", struct_prefix, "n", n, "_p", p, suffix)
    tian_var <- paste0("r_tian_", struct_prefix, "n", n, "_p", p, suffix)
    
    if (exists(r_var, envir = env_sim)) {
      r_val <- get(r_var, envir = env_sim)
      mat[i, 1] <- colMeans(r_val)[1]
      mat[i, 2] <- colMeans(r_val)[3]
    }
    
    if (exists(ols_var, envir = env_ols)) {
      ols_val <- get(ols_var, envir = env_ols)
      mat[i, 3] <- colMeans(ols_val)[1] 
    }
    
    if (exists(psi_var, envir = env_psi)) {
      psi_val <- get(psi_var, envir = env_psi)
      mat[i, 4] <- colMeans(psi_val)[1]
    }
    
    if (exists(tian_var, envir = env_tian)) {
      tian_val <- get(tian_var, envir = env_tian)
      mat[i, 5] <- colMeans(tian_val)[1]
    }
  }
  return(mat)
}

p_seq_200 <- c(200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400)
p_seq_500 <- c(500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

method_levels <- c("Alg1", "Alg2", "OLS_BH", "PSI_BH", "Tian_BH")
points <- c("Alg1" = 15, "Alg2" = 16, "OLS_BH" = 17, "PSI_BH" = 18, "Tian_BH" = 25)
line_types <- c("Alg1" = "solid", "Alg2" = "longdash", "OLS_BH" = "dashed", "PSI_BH" = "dotted", "Tian_BH" = "dotdash")

plot_fdr <- function(mat, p_seq, title, save_name) {
  type_vec <- rep(method_levels, each = length(p_seq))
  x_vals   <- rep(p_seq, times = 5)
  fdr_vals <- as.vector(mat)
  df <- data.frame(x = x_vals, Method = factor(type_vec, levels = method_levels), fdr = fdr_vals)
  
  p <- ggplot(data = df, mapping = aes(x = x, y = fdr, colour = Method, shape = Method, linetype = Method)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_x_continuous(name = "d", breaks = p_seq, labels = as.character(p_seq)) +
    scale_y_continuous(name = "Simulated FDR",
                       breaks = seq(0, 0.5, by = 0.1),
                       labels = c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50"),
                       limits = c(0, 0.5)) +
    ggtitle(title) +
    theme(legend.position = 'none', 
          plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
          axis.title = element_text(size = 26),
          axis.text = element_text(size = 22))
  
  ggsave(paste0(save_name, ".pdf"), p, device = "pdf", path = "./SI_plot_S1", width = 8, height = 6)
  return(p)
}

plot_fdr(get_fdr_matrix(0.1, "ind", 200, p_seq_200), p_seq_200, "Setting 1, n = 200", "ind_01")
plot_fdr(get_fdr_matrix(0.1, "AR1", 200, p_seq_200), p_seq_200, "Setting 2, n = 200", "AR1_01")
plot_fdr(get_fdr_matrix(0.1, "block", 200, p_seq_200), p_seq_200, "Setting 3, n = 200", "block_01")

plot_fdr(get_fdr_matrix(0.1, "ind", 500, p_seq_500), p_seq_500, "Setting 1, n = 500", "indL_01")
plot_fdr(get_fdr_matrix(0.1, "AR1", 500, p_seq_500), p_seq_500, "Setting 2, n = 500", "AR1L_01")
plot_fdr(get_fdr_matrix(0.1, "block", 500, p_seq_500), p_seq_500, "Setting 3, n = 500", "blockL_01")

plot_fdr(get_fdr_matrix(0.05, "ind", 200, p_seq_200), p_seq_200, "Setting 1, n = 200", "ind_005")
plot_fdr(get_fdr_matrix(0.05, "AR1", 200, p_seq_200), p_seq_200, "Setting 2, n = 200", "AR1_005")
plot_fdr(get_fdr_matrix(0.05, "block", 200, p_seq_200), p_seq_200, "Setting 3, n = 200", "block_005")

plot_fdr(get_fdr_matrix(0.05, "ind", 500, p_seq_500), p_seq_500, "Setting 1, n = 500", "indL_005")
plot_fdr(get_fdr_matrix(0.05, "AR1", 500, p_seq_500), p_seq_500, "Setting 2, n = 500", "AR1L_005")
plot_fdr(get_fdr_matrix(0.05, "block", 500, p_seq_500), p_seq_500, "Setting 3, n = 500", "blockL_005")