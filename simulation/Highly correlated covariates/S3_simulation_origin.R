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

n_values <- c(1000, 800, 600)
k <- 15
amp_values <- seq(0.1, 0.5, by = 0.05)
rho <- 0.67
fdr_configs <- list(
  list(fdr = 0.05, suffix = "_A005"),
  list(fdr = 0.1,  suffix = ""),
  list(fdr = 0.2,  suffix = "_A02")
)

run_simulation <- function(n, p, k, amp, fdr, suffix) {
  var_name <- paste0("r_n", n, "_p", p, "_amp", amp, suffix)
  
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

# Loop for OBBH
for (n in n_values) {
  p <- n
  for (amp in amp_values) {
    for (fdr_cfg in fdr_configs) {
      run_simulation(n, p, k, amp, fdr_cfg$fdr, fdr_cfg$suffix)
    }
  }
}

save.image(file = "./S3_simulation_origin.RData")
stopCluster(cl)