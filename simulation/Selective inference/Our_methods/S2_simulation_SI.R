library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../../source/Inverse.R",
  "../../source/MXK_KBH_KBBH.R"
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

run_simulation <- function(n, p, amp, fdr, suffix, rho, var_prefix) {
  var_name <- paste0(var_prefix, "n", n, "_p", p, "_amp", amp, suffix)
  k <- 0.04 * p
  
  assign(var_name, foreach(i = 1:300, .combine = rbind,
                           .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                           .options.RNG = 2024,
                           .export = c("source_files", "n", "p", "k", "amp", "fdr", "rho")) %dorng% {
                             
                             for (f in source_files) source(f)
                             
                             f_func <- get("K_BH_BBH")
                             f_func(n, p, k, amp, fdr, rho)
                             
                           }, envir = .GlobalEnv)
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp, 0.1, "", 0, "r_")
    }
  }
}
save.image(file = "./S2_simulation_indep.RData")

rm(list = ls(pattern = "^r_"))


for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp, 0.1, "", 0.4, "r_")
    }
  }
}
save.image(file = "./S2_simulation.RData")

rm(list = ls(pattern = "^r_"))


for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp, 0.05, "_A005", 0, "r_")
    }
  }
}
save.image(file = "./S2_simulation_indep_005.RData")

rm(list = ls(pattern = "^r_"))


for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_simulation(n, p, amp, 0.05, "_A005", 0.4, "r_")
    }
  }
}
save.image(file = "./S2_simulation_005.RData")

stopCluster(cl)