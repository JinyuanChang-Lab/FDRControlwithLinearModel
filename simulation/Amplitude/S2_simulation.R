library(glmnet)
library(knockoff)
library(doParallel)
library(scalreg)
library(doRNG)

source_files <- c(
  "../source/Inverse.R",
  "../source/MXK_KBH_KBBH.R",
  "../source/gm.R",
  "../source/GMtest.R",
  "../source/gm_new.R",
  "../source/GMtest_new.R"
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
fdr_configs <- list(
  list(fdr = 0.05, suffix = "_A005"),
  list(fdr = 0.1,  suffix = ""),
  list(fdr = 0.2,  suffix = "_A02")
)

kbh_config <- list(list(rho = 0.4, func_str = "K_BH_BBH"))
gm_config <- list(list(code = 2))

run_simulation <- function(n, p, amp, fdr, suffix, cov_list, func_name, var_prefix) {
  for (cov in cov_list) {
    var_name <- paste0(var_prefix, "n", n, "_p", p, "_amp", amp, suffix)
    k <- 0.04 * p
    
    loc_n <- n
    loc_p <- p
    loc_k <- k
    loc_amp <- amp
    loc_fdr <- fdr
    loc_cov <- cov
    loc_func_name <- func_name
    loc_src <- source_files
    
    assign(var_name, foreach(i = 1:100, .combine = rbind,
                             .packages = c('knockoff', 'scalreg', 'glmnet', 'doParallel', 'foreach'),
                             .options.RNG = 2024,
                             .export = c("source_files")) %dorng% {
                               for (f in loc_src) source(f)
                               if (loc_func_name == "KBH") {
                                 f <- get(loc_cov$func_str)
                                 f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, loc_cov$rho)
                               } else if (loc_func_name == "GMtest") {
                                 f <- get("GMtest")
                                 f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, loc_cov$code)
                               } else if (loc_func_name == "GMtest_new") {
                                 f <- get("GMtest_new")
                                 f(loc_n, loc_p, loc_k, loc_amp, loc_fdr, loc_cov$code)
                               }
                             }, envir = .GlobalEnv)
  }
}

for (n in n_values) {
  p_values <- p_list[[as.character(n)]]
  for (p in p_values) {
    for (amp in amp_values) {
      for (fdr_cfg in fdr_configs) {
        run_simulation(n, p, amp, fdr_cfg$fdr, fdr_cfg$suffix, kbh_config, "KBH", "r_")
        run_simulation(n, p, amp, fdr_cfg$fdr, fdr_cfg$suffix, gm_config, "GMtest", "rgm_")
        run_simulation(n, p, amp, fdr_cfg$fdr, fdr_cfg$suffix, gm_config, "GMtest_new", "rNgm_")
      }
    }
  }
}

save.image(file = "./S2_simulation.RData")
stopCluster(cl)