library(ggplot2)

env_sim <- new.env()
env_ds <- new.env()
load("./S2_plot/S2_simulation.RData", envir = env_sim)
load("./S2_plot_ds/S2_simulation_ds.RData", envir = env_ds)

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
replace_na_with_0(env_ds)

get_matrices <- function(n, p) {
  amp_seq <- seq(0.1, 0.5, by = 0.05)
  
  fdr_mat <- matrix(nrow = length(amp_seq), ncol = 4)
  power_mat <- matrix(nrow = length(amp_seq), ncol = 4)
  
  for (i in seq_along(amp_seq)) {
    amp <- amp_seq[i]
    r_var <- paste0("r_n", n, "_p", p, "_amp", amp)
    
    val_ds <- get(r_var, envir = env_ds)
    val_sim <- get(r_var, envir = env_sim)
    
    fdr_mat[i, ] <- c(colMeans(val_sim)[c(1, 3)], colMeans(val_ds)[c(1, 3)])
    power_mat[i, ] <- c(colMeans(val_sim)[c(2, 4)], colMeans(val_ds)[c(2, 4)])
  }
  return(list(fdr = fdr_mat, power = power_mat))
}

amp_seq <- seq(0.1, 0.5, by = 0.05)
type <- rep(c('A', 'B', 'C', 'D'), each = length(amp_seq))
points <- c('A' = 15, 'B' = 16, 'C' = 17, 'D' = 18)
line_types <- c("solid", "solid", "dotted", "dashed")
x_vals <- rep(amp_seq, times = 4)

plot_result <- function(fdr_mat, power_mat, n, p, save_path) {
  # FDR Plot
  fdr_vals <- c(fdr_mat[, 1], fdr_mat[, 2], fdr_mat[, 3], fdr_mat[, 4])
  df_fdr <- data.frame(x = x_vals, type = type, y = fdr_vals)
  
  p_fdr <- ggplot(data = df_fdr, mapping = aes(x = x, y = y, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "Amplitude", breaks = seq(0.1, 0.5, by = 0.05),
                       labels = c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"),
                       limits = c(0.1, 0.5)) +
    scale_y_continuous(name = "Simulated FDR", breaks = seq(0, 0.8, by = 0.1),
                       labels = c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80"),
                       limits = c(0, 0.8)) +
    ggtitle(sprintf("n = %d, d = %d", n, p)) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(sprintf("n%d_d%d_fdr.pdf", n, p), p_fdr, device = "pdf", path = save_path, width = 8, height = 6)
  
  # Power Plot
  power_vals <- c(power_mat[, 1], power_mat[, 2], power_mat[, 3], power_mat[, 4])
  df_power <- data.frame(x = x_vals, type = type, y = power_vals)
  
  p_power <- ggplot(data = df_power, mapping = aes(x = x, y = y, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "Amplitude", breaks = seq(0.1, 0.5, by = 0.05),
                       labels = c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"),
                       limits = c(0.1, 0.5)) +
    scale_y_continuous(name = "Simulated power", breaks = seq(0, 1, by = 0.25),
                       labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                       limits = c(0, 1)) +
    ggtitle(sprintf("n = %d, d = %d", n, p)) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(sprintf("n%d_d%d_power.pdf", n, p), p_power, device = "pdf", path = save_path, width = 8, height = 6)
}

n_list <- list("200" = c(200, 300, 400), "500" = c(500, 750, 1000))
save_path <- "./S2_plot_ds"

for (n_str in names(n_list)) {
  n <- as.integer(n_str)
  p_values <- n_list[[n_str]]
  for (p in p_values) {
    res <- get_matrices(n, p)
    plot_result(res$fdr, res$power, n, p, save_path)
  }
}