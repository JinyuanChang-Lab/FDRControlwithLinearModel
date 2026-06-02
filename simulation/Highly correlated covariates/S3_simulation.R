library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
  "../source/MXK_KBH_KBBH.R",
  "../source/gm.R",
  "../source/GMtest067.R",
  "../source/gm_new.R",
  "../source/GMtest067_new.R"
)
for (f in source_files) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(1000, 800, 600)
k <- 15
amp_values <- seq(0.1, 0.5, by = 0.05)
rho <- 0.67
fdr_configs <- list(
  list(fdr = 0.05, suffix = "_A005"),
  list(fdr = 0.1,  suffix = ""),
  list(fdr = 0.2,  suffix = "_A02")
)

run_simulation <- function(n, p, k, amp, fdr, suffix, func_str, var_prefix, extra_param) {
  var_name <- paste0(var_prefix, "n", n, "_p", p, "_amp", amp, suffix)
  
  loc_n <- n
  loc_p <- p
  loc_k <- k
  loc_amp <- amp
  loc_fdr <- fdr
  loc_extra <- extra_param
  loc_func_str <- func_str
  loc_src <- source_files
  
  assign(var_name, foreach(i = 1:100, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                           .options.RNG = 2024,
                           .export = c("source_files")) %dorng% {
                             for (f in loc_src) source(f)
                             f <- get(loc_func_str)
                             f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, loc_extra)
                           }, envir = .GlobalEnv)
}

for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "K_BH_BBH", "r_", rho)
    }
  }
}

for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "GMtest067", "rgm_", 2)
    }
  }
}

for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "GMtest067_new", "rNgm_", 2)
    }
  }
}

save.image(file = "./S3_simulation.RData")
stopCluster(cl)