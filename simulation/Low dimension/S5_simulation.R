library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/MXK_KBH_KBBH.R",
  "../source/BBH.R",
  "../source/Inverse.R"
)
for (f in source_files) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

rho <- 0.4
npk_list <- list(
  c(150, 30, 3),
  c(200, 40, 4),
  c(250, 50, 5)
)
amp_values <- seq(0.1, 0.5, by = 0.05)
fdr_configs <- list(
  list(fdr = 0.05, suffix = "_A005"),
  list(fdr = 0.1,  suffix = ""),
  list(fdr = 0.2,  suffix = "_A02")
)

run_simulation <- function(n, p, k, amp, fdr, suffix, func_str, var_prefix) {
  amp_str <- gsub("\\.", "", sprintf("%.2f", amp))
  var_name <- paste0(var_prefix, "n", n, "_p", p, "_k", k, "_amp", amp_str, suffix)
  
  loc_n <- n
  loc_p <- p
  loc_k <- k
  loc_amp <- amp
  loc_fdr <- fdr
  loc_func_str <- func_str
  loc_src <- source_files
  
  assign(var_name, foreach(i = 1:100, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                           .options.RNG = 2024,
                           .export = c("source_files")) %dorng% {
                             for (f in loc_src) source(f)
                             f <- get(loc_func_str)
                             f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, rho)
                           }, envir = .GlobalEnv)
}

for (npk in npk_list) {
  n <- npk[1]
  p <- npk[2]
  k <- npk[3]
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "K_BH_BBH", "r_")
    }
  }
}

for (npk in npk_list) {
  n <- npk[1]
  p <- npk[2]
  k <- npk[3]
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "BBH", "r_bh_")
    }
  }
}

save.image(file = "./S5_simulation.RData")
stopCluster(cl)