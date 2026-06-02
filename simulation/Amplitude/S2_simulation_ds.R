library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)
library(expm)

source_files <- c(
  "../source/Inverse.R",
  "../source/DSBBH.R"
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
alpha <- 0.1
rho <- 0.4

run_simulation <- function(n, p, amp) {
  var_name <- paste0("r_n", n, "_p", p, "_amp", amp)
  n_ds <- 2 * n
  k <- 0.04 * p
  
  loc_n_ds <- n_ds
  loc_p <- p
  loc_k <- k
  loc_amp <- amp
  loc_alpha <- alpha
  loc_rho <- rho
  loc_src <- source_files
  
  assign(var_name, foreach(i = 1:100, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'expm'),
                           .options.RNG = 2024,
                           .export = c("source_files")) %dorng% {
                             for (f in loc_src) source(f)
                             f <- get("DSBBH")
                             f(loc_n_ds, loc_p, loc_k, loc_amp, loc_alpha, loc_rho)
                           }, envir = .GlobalEnv)
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp)
    }
  }
}

save.image(file = "./S2_simulation_ds.RData")
stopCluster(cl)