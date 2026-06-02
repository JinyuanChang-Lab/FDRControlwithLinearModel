library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files_orig <- c(
  "../source/Inverse.R",
  "../source/OBBH.R",
  "../source/OBBH_E.R"
)
for (f in source_files_orig) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

amplitude <- 0
base_alpha <- 0.05
cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(200, 500)
p_list <- list(
  "200" = c(200, 240, 280, 320, 360, 400),
  "500" = c(500, 600, 700, 800, 900, 1000)
)

fdr_list <- list(
  list(mult = 1, suffix = "_A005", extra_prefix = ""),
  list(mult = 2, suffix = "",      extra_prefix = ""),
  list(mult = 4, suffix = "_A02",  extra_prefix = "")
)

cov_list <- list(
  list(prefix = "",  rho = 0.4, func_str = "OBBH"),
  list(prefix = "0", rho = 0,   func_str = "OBBH"),
  list(prefix = "E", rho = 0.2, func_str = "OBBH_E")
)

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  
  for (p in p_values) {
    for (fdr in fdr_list) {
      for (cov in cov_list) {
        
        var_prefix_full <- paste0(fdr$extra_prefix, cov$prefix)
        var_name <- paste0("r_", var_prefix_full, "n", n, "_p", p, fdr$suffix)
        
        loc_n <- n
        loc_p <- p
        loc_amp <- amplitude
        loc_alpha <- base_alpha * fdr$mult
        loc_rho <- cov$rho
        loc_func_str <- cov$func_str
        loc_src <- source_files_orig
        
        assign(var_name, foreach(i = 1:100, .combine = rbind,
                                 .packages = c('knockoff', 'scalreg', 'glmnet'),
                                 .options.RNG = 2024) %dorng% {
                                   
                                   for (f in loc_src) source(f)
                                   f <- get(loc_func_str)
                                   f(loc_n, loc_p, 0, loc_amp, loc_alpha, loc_rho)
                                   
                                 }, envir = .GlobalEnv)
      }
    }
  }
}

save.image(file = "./S1_simulation_origin.RData")
stopCluster(cl)