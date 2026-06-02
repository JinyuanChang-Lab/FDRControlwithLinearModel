library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
  "../source/OBBH.R"
)
for (f in source_files) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(200, 500)
p_list <- list(
  "200" = c(200, 300, 400),
  "500" = c(500, 750, 1000)
)
amp_values <- seq(0.1, 0.5, by = 0.05)
alpha_base <- 0.05
rho <- 0.4
fdr_configs <- list(
  list(fdr = alpha_base, suffix = "_A005"),
  list(fdr = 2 * alpha_base, suffix = ""),
  list(fdr = 4 * alpha_base, suffix = "_A02")
)

run_simulation <- function(n, p, amp, fdr, suffix) {
  var_name <- paste0("r_n", n, "_p", p, "_amp", amp, suffix)
  k <- 0.04 * p
  
  loc_n <- n
  loc_p <- p
  loc_k <- k
  loc_amp <- amp
  loc_fdr <- fdr
  loc_rho <- rho
  loc_src <- source_files
  
  assign(var_name, foreach(i = 1:100, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                           .options.RNG = 2024,
                           .export = c("source_files")) %dorng% {
                             for (f in loc_src) source(f)
                             f <- get("OBBH")
                             f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, loc_rho)
                           }, envir = .GlobalEnv)
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      for (fdr_cfg in fdr_configs) {
        run_simulation(n, p, amp, fdr_cfg$fdr, fdr_cfg$suffix)
      }
    }
  }
}

save.image(file = "./S2_simulation_origin.RData")
stopCluster(cl)