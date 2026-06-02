library(ggplot2)

env_ds <- new.env()
env_sim <- new.env()
load("./S1_plot_ds/S1_simulation_ds.RData", envir = env_ds)
load("./S1_plot/S1_simulation.RData", envir = env_sim)

replace_na_with_0 <- function(env) {
  for (var_name in ls(envir = env)) {
    var <- get(var_name, envir = env)
    if (is.matrix(var) || is.data.frame(var)) {
      var[is.na(var)] <- 0
      assign(var_name, var, envir = env)
    }
  }
}
replace_na_with_0(env_ds)
replace_na_with_0(env_sim)

get_fdr_matrix <- function(structure, n, p_seq) {
  mat <- matrix(nrow = length(p_seq), ncol = 4)
  for (i in seq_along(p_seq)) {
    p <- p_seq[i]
    prefix <- switch(structure,
                     ind = paste0("r_0n", n, "_p", p),
                     AR1 = paste0("r_n", n, "_p", p),
                     block = paste0("r_En", n, "_p", p))
    val_sim <- get(prefix, envir = env_sim)
    val_ds <- get(prefix, envir = env_ds)
    mat[i, ] <- c(colMeans(val_sim)[c(1, 3)], colMeans(val_ds)[c(1, 3)])
  }
  return(mat)
}

p_seq_200 <- c(200, 240, 280, 320, 360, 400)
p_seq_500 <- c(500, 600, 700, 800, 900, 1000)
points <- c('A' = 15, 'B' = 16, 'C' = 17, 'D' = 18)

plot_fdr <- function(mat, p_seq, title, save_name) {
  type <- rep(c('A', 'B', 'C', 'D'), each = length(p_seq))
  x_vals <- rep(p_seq, times = 4)
  fdr_vals <- c(mat[, 1], mat[, 2], mat[, 3], mat[, 4])
  df <- data.frame(x = x_vals, type = type, fdr = fdr_vals)
  
  p <- ggplot(data = df, mapping = aes(x = x, y = fdr, colour = type, shape = type, linetype = type)) +
    scale_linetype_manual(values = c("solid", "solid", "dotted", "dashed")) +
    scale_shape_manual(values = points) +
    geom_line(linewidth = 2) +
    geom_point(size = 4.5) +
    scale_x_continuous(name = "d", breaks = p_seq, labels = as.character(p_seq)) +
    scale_y_continuous(name = "Simulated FDR",
                       breaks = seq(0, 0.5, by = 0.1),
                       labels = c("0.00", "0.10", "0.20", "0.30", "0.40", "0.50"),
                       limits = c(0, 0.5)) +
    ggtitle(title) +
    theme(legend.position = 'none', text = element_text(size = 27), plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(save_name, ".pdf"), p, device = "pdf", path = "./S1_plot_ds", width = 8, height = 6)
  return(p)
}

plot_fdr(get_fdr_matrix("ind", 200, p_seq_200), p_seq_200, "Setting 1, n = 200", "ind_01")
plot_fdr(get_fdr_matrix("AR1", 200, p_seq_200), p_seq_200, "Setting 2, n = 200", "AR1_01")
plot_fdr(get_fdr_matrix("block", 200, p_seq_200), p_seq_200, "Setting 3, n = 200", "block_01")

plot_fdr(get_fdr_matrix("ind", 500, p_seq_500), p_seq_500, "Setting 1, n = 500", "indL_01")
plot_fdr(get_fdr_matrix("AR1", 500, p_seq_500), p_seq_500, "Setting 2, n = 500", "AR1L_01")
plot_fdr(get_fdr_matrix("block", 500, p_seq_500), p_seq_500, "Setting 3, n = 500", "blockL_01")