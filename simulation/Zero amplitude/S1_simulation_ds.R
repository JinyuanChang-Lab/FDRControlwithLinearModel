library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)
library(expm)

source_files_ds <- c(
  "../source/Inverse.R",
  "../source/DSBBH.R",
  "../source/DSBBH_E.R"
)
for (f in source_files_ds) source(f)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2024)

amplitude <- 0
alpha <- 0.1
cl <- makeCluster(32, type = "SOCK")
registerDoParallel(cl)

n_values <- c(200, 500)
p_list <- list(
  "200" = c(200, 240, 280, 320, 360, 400),
  "500" = c(500, 600, 700, 800, 900, 1000)
)

ds_cov <- list(
  list(prefix = "", rho = 0.4, func_str = "DSBBH"),
  list(prefix = "0", rho = 0, func_str = "DSBBH"),
  list(prefix = "E", rho = 0.2, func_str = "DSBBH_E")
)

for (n_base in n_values) {
  n_ds <- 2 * n_base  
  p_values <- p_list[[as.character(n_base)]]
  
  for (p in p_values) {
    for (cov in ds_cov) {
      var_name <- paste0("r_", cov$prefix, "n", n_base, "_p", p)
      
      loc_n <- n_ds
      loc_p <- p
      loc_amp <- amplitude
      loc_alpha <- alpha
      loc_rho <- cov$rho
      loc_func_str <- cov$func_str
      loc_src <- source_files_ds
      
      assign(var_name, foreach(i = 1:100, .combine = rbind,
                               .packages = c('knockoff', 'scalreg', 'glmnet', 'expm'),
                               .options.RNG = 2024) %dorng% {
                                 
                                 for (f in loc_src) source(f)
                                 f <- get(loc_func_str)
                                 f(loc_n, loc_p, 0, loc_amp, loc_alpha, loc_rho)
                                 
                               }, envir = .GlobalEnv)
    }
  }
}

save.image(file = "./S1_simulation_ds.RData")
stopCluster(cl)