library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../../source/Inverse.R",
  "../../source/MXK_KBH_KBBH.R",
  "../../source/MXK_KBH_KBBH_E.R"
)

for (f in source_files) source(f)

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

kbh_cov <- list(
  list(prefix = "", rho = 0.4, func_str = "K_BH_BBH"),
  list(prefix = "0", rho = 0, func_str = "K_BH_BBH"),
  list(prefix = "E", rho = 0.2, func_str = "K_BH_BBH_E")
)

run_simulation <- function(n, p, fdr, suffix, cov_list, var_prefix) {
  for (cov in cov_list) {
    var_name <- paste0(var_prefix, cov$prefix, "n", n, "_p", p, suffix)
    
    assign(var_name, foreach(i = 1:300, .combine = rbind,
                             .packages = c('knockoff', 'scalreg', 'glmnet'),
                             .options.RNG = 2024,
                             .export = c("amplitude", "source_files")) %dorng% {
                               
                               for (f in source_files) source(f)
                               
                               f <- get(cov$func_str)
                               f(n, p, 0, amplitude, fdr, cov$rho)
                               
                             }, envir = .GlobalEnv)
  }
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    run_simulation(n, p, 0.05, "_A005", kbh_cov, "r_")
  }
}

save.image(file = "./S1_simulation_005.RData")

rm(list = ls(pattern = "_A005$"))

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    run_simulation(n, p, 0.1, "", kbh_cov, "r_")
  }
}

save.image(file = "./S1_simulation.RData")

stopCluster(cl)