library(glmnet)
library(doParallel)
library(doRNG)
library(MASS)

source("../../source/OLS_BH.R")

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

amplitude <- 0

cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(200, 500)
p_list <- list(
  "200" = c(200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400),
  "500" = c(500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)
)

ols_cov <- list(
  list(prefix = "",  rho = 0.4),
  list(prefix = "0", rho = 0),
  list(prefix = "E", rho = 0.2)
)

run_ols <- function(n, p, fdr, suffix, cov_list, var_prefix) {
  for (cov in cov_list) {
    var_name <- paste0(var_prefix, cov$prefix, "n", n, "_p", p, suffix)
    
    assign(var_name, 
           foreach(i = 1:300, .combine = rbind,
                   .packages = c("glmnet", "MASS"),
                   .options.RNG = 2024,
                   .export = c("amplitude", "OLS_BH")) %dorng% {
                     
                     OLS_BH(n, p, 0, amplitude, fdr, cov$rho)
                     
                   }, envir = .GlobalEnv)
  }
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    run_ols(n, p, 0.1, "", ols_cov, "r_ols_")
  }
}

save.image(file = "./S1_simulation_OLS_BH.RData")

rm(list = ls(pattern = "^r_ols_"))

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    run_ols(n, p, 0.05, "_A005", ols_cov, "r_ols_")
  }
}

save.image(file = "./S1_simulation_OLS_BH_005.RData")

stopCluster(cl)