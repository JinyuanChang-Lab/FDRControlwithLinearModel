library(ggplot2)

env <- new.env()
load("./S5_plot/S5_simulation.RData", envir = env)

replace_na_with_0 <- function(env) {
  for (var_name in ls(envir = env)) {
    var <- get(var_name, envir = env)
    if (is.matrix(var) || is.data.frame(var)) {
      var[is.na(var)] <- 0
      assign(var_name, var, envir = env)
    }
  }
}
replace_na_with_0(env)

get_matrices <- function(n, p, k, fdr_level) {
  suffix <- switch(as.character(fdr_level),
                   "0.05" = "_A005",
                   "0.1" = "",
                   "0.2" = "_A02")
  amp_seq <- seq(0.1, 0.5, by = 0.05)
  
  fdr_mat <- matrix(nrow = length(amp_seq), ncol = 3)
  power_mat <- matrix(nrow = length(amp_seq), ncol = 3)
  
  fdr_cols <- c(1, 3, 1)
  power_cols <- c(2, 4, 2)
  
  for (i in seq_along(amp_seq)) {
    amp <- amp_seq[i]
    amp_str <- gsub("\\.", "", sprintf("%.2f", amp))
    r_var <- paste0("r_n", n, "_p", p, "_k", k, "_amp", amp_str, suffix)
    rbh_var <- paste0("r_bh_n", n, "_p", p, "_k", k, "_amp", amp_str, suffix)
    
    r_val <- get(r_var, envir = env)
    rbh_val <- get(rbh_var, envir = env)
    
    fdr_mat[i, ] <- c(colMeans(r_val)[fdr_cols[1:2]], colMeans(rbh_val)[fdr_cols[3]])
    power_mat[i, ] <- c(colMeans(r_val)[power_cols[1:2]], colMeans(rbh_val)[power_cols[3]])
  }
  
  fdr_mat[is.na(fdr_mat)] <- 0
  power_mat[is.na(power_mat)] <- 0
  return(list(fdr = fdr_mat, power = power_mat))
}

amp_seq <- seq(0.1, 0.5, by = 0.05)
type <- rep(c('A', 'B', 'C'), each = length(amp_seq))
points <- c('A' = 15, 'B' = 16, 'C' = 17)
line_types <- c("solid", "solid", "dotted")
x_vals <- rep(amp_seq, times = 3)

plot_result <- function(fdr_mat, power_mat, n, p, fdr_str, save_path = "./S5_plot") {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # FDR Plot
  fdr_vals <- c(fdr_mat[, 1], fdr_mat[, 2], fdr_mat[, 3])
  df_fdr <- data.frame(x = x_vals, type = type, y = fdr_vals)
  
  p_fdr <- ggplot(data = df_fdr, mapping = aes(x = x, y = y, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "Amplitude", breaks = amp_seq,
                       labels = sprintf("%.2f", amp_seq),
                       limits = range(amp_seq)) +
    scale_y_continuous(name = "Simulated FDR", breaks = seq(0, 0.5, by = 0.1),
                       labels = sprintf("%.2f", seq(0, 0.5, by = 0.1)),
                       limits = c(0, 0.5)) +
    ggtitle(sprintf("n = %d, d = %d", n, p)) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(sprintf("n%dp%dfdr_%s.pdf", n, p, fdr_str), p_fdr, device = "pdf", path = save_path, width = 8, height = 6)
  
  # Power Plot
  power_vals <- c(power_mat[, 1], power_mat[, 2], power_mat[, 3])
  df_power <- data.frame(x = x_vals, type = type, y = power_vals)
  
  p_power <- ggplot(data = df_power, mapping = aes(x = x, y = y, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = line_types) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "Amplitude", breaks = amp_seq,
                       labels = sprintf("%.2f", amp_seq),
                       limits = range(amp_seq)) +
    scale_y_continuous(name = "Simulated power", breaks = seq(0, 1, by = 0.25),
                       labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
                       limits = c(0, 1)) +
    ggtitle(sprintf("n = %d, d = %d", n, p)) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(sprintf("n%dp%dpower_%s.pdf", n, p, fdr_str), p_power, device = "pdf", path = save_path, width = 8, height = 6)
}

npk_list <- list(list(n=150, p=30, k=3), list(n=200, p=40, k=4), list(n=250, p=50, k=5))
fdr_list <- list(list(level = 0.05, str = "005"), list(level = 0.1, str = "01"), list(level = 0.2, str = "02"))

for (npk in npk_list) {
  n <- npk$n
  p <- npk$p
  k <- npk$k
  for (fdr_cfg in fdr_list) {
    res <- get_matrices(n, p, k, fdr_cfg$level)
    plot_result(res$fdr, res$power, n, p, fdr_cfg$str)
  }
}