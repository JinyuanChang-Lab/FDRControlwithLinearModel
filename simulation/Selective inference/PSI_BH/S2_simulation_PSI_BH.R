library(glmnet)
library(selectiveInference)
library(doParallel)
library(doRNG)
library(MASS)

source("../../source/PSI_BH.R")

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

psi_cov <- list(
  list(rho = 0.4, label = "ar1"),
  list(rho = 0,   label = "indep")
)

run_psi <- function(n, p, amp, fdr, cov_list, var_prefix, suffix) {
  for (cov in cov_list) {
    var_name <- paste0(var_prefix, cov$label, "_n", n, "_p", p, "_amp", amp, suffix)
    k <- 0.04 * p
    loc_rho <- cov$rho
    
    assign(var_name, foreach(i = 1:300, .combine = rbind,
                             .packages = c("glmnet", "scalreg", "selectiveInference", "MASS"),
                             .options.RNG = 2024,
                             .export = c("PSI_BH")) %dorng% {
                               PSI_BH(n, p, k, amp, fdr, loc_rho)
                             }, envir = .GlobalEnv)
  }
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_psi(n, p, amp, 0.1, psi_cov, "r_psi_", "")
    }
  }
}

save.image(file = "./S2_simulation_PSI_BH.RData")

rm(list = ls(pattern = "^r_psi_"))

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      run_psi(n, p, amp, 0.05, psi_cov, "r_psi_", "_A005")
    }
  }
}

save.image(file = "./S2_simulation_PSI_BH_005.RData")

stopCluster(cl)