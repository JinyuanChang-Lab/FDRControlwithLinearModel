library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
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

n_values <- c(200)
p_list <- list(
  "200" = c(200, 300, 400)
)
amp_values <- seq(0.1, 0.5, by = 0.05)
alpha <- 0.1
rho <- 0.4

run_simulation <- function(n, p, amp, func_str, var_prefix, extra_param) {
  var_name <- paste0(var_prefix, "n", n, "_p", p, "_amp", amp)
  k <- 0.04 * p
  
  loc_n <- n
  loc_p <- p
  loc_k <- k
  loc_amp <- amp
  loc_alpha <- alpha
  loc_extra <- extra_param
  loc_func_str <- func_str
  loc_src <- source_files
  
  assign(var_name, foreach(i = 1:100, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                           .options.RNG = 2024,
                           .export = c("source_files")) %dorng% {
                             for (f in loc_src) source(f)
                             f <- get(loc_func_str)
                             f(loc_n, loc_p, loc_k, loc_amp, loc_alpha, loc_extra)
                           }, envir = .GlobalEnv)
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp, "TBBH", "r_", rho)
      run_simulation(n, p, amp, "GMtestT", "rgm_", 2)
      run_simulation(n, p, amp, "GMtestT_new", "rNgm_", 2)
    }
  }
}

save.image(file = "./S2_simulation_nongaussian.RData")
stopCluster(cl)