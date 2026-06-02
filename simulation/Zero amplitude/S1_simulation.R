library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
  "../source/MXK_KBH_KBBH.R",
  "../source/MXK_KBH_KBBH_E.R",
  "../source/gm.R",
  "../source/GMtest.R",
  "../source/gm_new.R",
  "../source/GMtest_new.R"
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
gm_cov <- list(
  list(prefix = "0", code = 1),
  list(prefix = "", code = 2),
  list(prefix = "E", code = 3)
)

run_simulation <- function(n, p, fdr, suffix, cov_list, func_name, var_prefix) {
  for (cov in cov_list) {
    var_name <- paste0(var_prefix, cov$prefix, "n", n, "_p", p, suffix)
    
    assign(var_name, foreach(i = 1:100, .combine = rbind,
                             .packages = c('knockoff', 'scalreg', 'glmnet'),
                             .options.RNG = 2024,
                             .export = c("amplitude", "source_files")) %dorng% {
                               
                               for (f in source_files) source(f)
                               
                               if (func_name == "KBH") {
                                 f <- get(cov$func_str)
                                 f(n, p, 0, amplitude, fdr, cov$rho)
                               } else {
                                 f <- get(func_name)
                                 f(n, p, 0, amplitude, fdr, cov$code)
                               }
                             }, envir = .GlobalEnv)
  }
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    run_simulation(n, p, 0.05, "_A005", kbh_cov, "KBH", "r_")
    run_simulation(n, p, 0.1, "", kbh_cov, "KBH", "r_")
    run_simulation(n, p, 0.2, "_A02", kbh_cov, "KBH", "r_")
    
    run_simulation(n, p, 0.05, "_A005", gm_cov, "GMtest", "rgm_")
    run_simulation(n, p, 0.1, "", gm_cov, "GMtest", "rgm_")
    run_simulation(n, p, 0.2, "_A02", gm_cov, "GMtest", "rgm_")
    
    run_simulation(n, p, 0.05, "_A005", gm_cov, "GMtest_new", "rNgm_")
    run_simulation(n, p, 0.1, "", gm_cov, "GMtest_new", "rNgm_")
    run_simulation(n, p, 0.2, "_A02", gm_cov, "GMtest_new", "rNgm_")
  }
}

save.image(file = "./S1_simulation.RData")
stopCluster(cl)