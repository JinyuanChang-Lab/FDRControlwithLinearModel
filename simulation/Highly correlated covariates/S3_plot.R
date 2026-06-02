library(ggplot2)

env_sim <- new.env()
load("./S3_plot/S3_simulation.RData", envir = env_sim)

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

get_matrices <- function(n, p, fdr_level) {
  suffix <- switch(as.character(fdr_level),
                   "0.05" = "_A005",
                   "0.1" = "",
                   "0.2" = "_A02")
  amp_seq <- seq(0.1, 0.5, by = 0.05)
  
  fdr_mat <- matrix(nrow = length(amp_seq), ncol = 5)
  power_mat <- matrix(nrow = length(amp_seq), ncol = 5)
  
  for (i in seq_along(amp_seq)) {
    amp <- amp_seq[i]
    r_var <- paste0("r_n", n, "_p", p, "_amp", amp, suffix)
    rgm_var <- paste0("rgm_n", n, "_p", p, "_amp", amp, suffix)
    rNgm_var <- paste0("rNgm_n", n, "_p", p, "_amp", amp, suffix)
    
    r_val <- get(r_var, envir = env_sim)
    rgm_val <- get(rgm_var, envir = env_sim)
    rNgm_val <- get(rNgm_var, envir = env_sim)
    
    fdr_mat[i, ] <- c(colMeans(r_val)[c(1, 3, 5)], colMeans(rgm_val)[1], colMeans(rNgm_val)[1])
    power_mat[i, ] <- c(colMeans(r_val)[c(2, 4, 6)], colMeans(rgm_val)[2], colMeans(rNgm_val)[2])
  }
  return(list(fdr = fdr_mat, power = power_mat))
}

amp_seq <- seq(0.1, 0.5, by = 0.05)
type <- rep(c('A', 'B', 'C', 'D', 'E'), each = length(amp_seq))
points <- c('A' = 15, 'B' = 16, 'C' = 17, 'D' = 18, 'E' = 15)
line_types <- c("solid", "solid", "dotted", "dashed", "dashed")
x_vals <- rep(amp_seq, times = 5)

plot_result <- function(fdr_mat, power_mat, n, p, fdr_str, save_path = "./S3_plot") {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # FDR Plot
  fdr_vals <- c(fdr_mat[, 1], fdr_mat[, 2], fdr_mat[, 3], fdr_mat[, 4], fdr_mat[, 5])
  df_fdr <- data.frame(x = x_vals, type = type, y = fdr_vals)
  
  p_fdr <- ggplot(data = df_fdr, mapping = aes(x = x, y = y, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "Amplitude", breaks = seq(0.1, 0.5, by = 0.05),
                       labels = c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"),
                       limits = c(0.1, 0.5)) +
    scale_y_continuous(name = "Simulated FDR", breaks = seq(0, 0.6, by = 0.1),
                       labels = c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60"),
                       limits = c(0, 0.6)) +
    ggtitle(sprintf("n = %d, d = %d", n, p)) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(sprintf("n%dp%dfdr_%s.pdf", n, p, fdr_str), p_fdr, device = "pdf", path = save_path, width = 8, height = 6)
  
  # Power Plot
  power_vals <- c(power_mat[, 1], power_mat[, 2], power_mat[, 3], power_mat[, 4], power_mat[, 5])
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
  
  ggsave(sprintf("n%dp%dpower_%s.pdf", n, p, fdr_str), p_power, device = "pdf", path = save_path, width = 8, height = 6)
}

n_list <- list("600" = 600, "800" = 800, "1000" = 1000)
fdr_list <- list(list(level = 0.05, str = "005"), list(level = 0.1, str = "01"), list(level = 0.2, str = "02"))

for (n_str in names(n_list)) {
  n <- as.integer(n_str)
  p <- n_list[[n_str]]
  for (fdr_cfg in fdr_list) {
    res <- get_matrices(n, p, fdr_cfg$level)
    plot_result(res$fdr, res$power, n, p, fdr_cfg$str)
  }
}