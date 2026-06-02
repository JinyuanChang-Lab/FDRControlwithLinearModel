library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files_ng <- c(
  "../source/Inverse.R",
  "../source/TBBH.R",
  "../source/TBBH_E.R",
  "../source/gm.R",
  "../source/GMtestT.R",
  "../source/gm_new.R",
  "../source/GMtestT_new.R"
)
for (f in source_files_ng) source(f)

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

method_list <- list(
  list(prefix = "r_", func_str = "TBBH", 
       covs = list(
         list(cov_prefix = "", rho = 0.4),
         list(cov_prefix = "0", rho = 0),
         list(cov_prefix = "E", rho = 0.2, func_suffix = "_E")
       )),
  list(prefix = "rgm_", func_str = "GMtestT",
       covs = list(
         list(cov_prefix = "0", code = 1),
         list(cov_prefix = "", code = 2),
         list(cov_prefix = "E", code = 3)
       )),
  list(prefix = "rNgm_", func_str = "GMtestT_new",
       covs = list(
         list(cov_prefix = "0", code = 1),
         list(cov_prefix = "", code = 2),
         list(cov_prefix = "E", code = 3)
       ))
)

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (meth in method_list) {
      for (cov in meth$covs) {
        var_name <- paste0(meth$prefix, cov$cov_prefix, "n", n, "_p", p)
        
        current_func <- meth$func_str
        if (!is.null(cov$func_suffix)) {
          current_func <- paste0(meth$func_str, cov$func_suffix)
        }
        
        loc_n <- n
        loc_p <- p
        loc_amp <- amplitude
        loc_alpha <- alpha
        loc_cov <- cov
        loc_func_str <- current_func
        loc_src <- source_files_ng
        
        assign(var_name, foreach(i = 1:100, .combine = rbind,
                                 .packages = c('knockoff', 'scalreg', 'glmnet'),
                                 .options.RNG = 2024) %dorng% {
                                   
                                   for (f in loc_src) source(f)
                                   f <- get(loc_func_str)
                                   
                                   if (!is.null(loc_cov$code)) {
                                     f(loc_n, loc_p, 0, loc_amp, loc_alpha, loc_cov$code)
                                   } else {
                                     f(loc_n, loc_p, 0, loc_amp, loc_alpha, loc_cov$rho)
                                   }
                                   
                                 }, envir = .GlobalEnv)
      }
    }
  }
}

save.image(file = "./S1_simulation_nongaussian.RData")
stopCluster(cl)