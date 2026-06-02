library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
  "../source/TBBHw_pre.R",
  "../source/TBBH.R",
  "../source/gm.R",
  "../source/GMtestT.R",
  "../source/gm_new.R",
  "../source/GMtestT_new.R"
)
for (f in source_files) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(400, 600, 800)
k <- 15
amp_values <- seq(0.1, 0.5, by = 0.05)
rho <- 0.4
alpha <- 0.1

fdr_configs <- list(
  list(fdr = alpha, suffix = "")
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

# Loop for TBBHw_pre
for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "TBBHw_pre", "r_", rho)
    }
  }
}

# Loop for TBBH
for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "TBBH", "r_O", rho)
    }
  }
}

# Loop for GMtestT
for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "GMtestT", "rgm_", 2)
    }
  }
}

# Loop for GMtestT_new
for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix, "GMtestT_new", "rNgm_", 2)
    }
  }
}

save.image(file = "./S4_simulation_nongaussian.RData")
stopCluster(cl)