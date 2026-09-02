# Standalone script — run manually to regenerate tests/testthat/fixtures/init-regression.rds
# Do NOT source() from test files.
setwd("/Users/sophiehuebler/Documents/bhCRR")

devtools::load_all(quiet = TRUE)

# -- data --
d <- make_ssl_fixture_data()
x <- d$x
y <- d$y

# -- init_ref: replicate fit_ssl_psdh.r lines 124-130 exactly --
init_lam_path <- 10^seq(log10(0.1), log10(0.001), length = 10)
full_model_cv <- cv_fastCrrp_cpp(x, y[, 1], y[, 2], k = 5, penalty = "LASSO",
                                 lambda_path = init_lam_path,
                                 tuning = "wolbers", eval_quantile = 0.5)
lambda_min <- full_model_cv$lambda_min
init_ref <- full_model_cv$full_model$coef[, full_model_cv$lambda == lambda_min]

# -- fit_ref --
fit <- fit_ssl_psdh(x, y,
  ss = c(0.04, 0.5), initial_sparsity = 0.05,
  theta_a = 1, theta_b = ncol(x),
  maxit = 50, epsilon = 1e-04,
  init = NULL,
  init_lam_path = init_lam_path,
  inner_maxit_start = 1000)

fit_ref <- list(
  Estimate     = fit$coefficients$Estimate,
  pips         = fit$pips,
  penalty.factor = fit$penalty.factor,
  iterations   = fit$iterations,
  conv         = fit$conv,
  ss           = fit$ss
)

# -- provenance --
provenance <- list(
  R_version        = R.version.string,
  fastcmprsk_version = as.character(packageVersion("fastcmprsk")),
  created_at       = Sys.time()
)

# -- save --
out_path <- "tests/testthat/fixtures/init-regression.rds"
saveRDS(list(x = x, y = y, init_ref = init_ref, fit_ref = fit_ref,
             provenance = provenance),
        out_path)
cat("Saved:", out_path, "\n")
cat("File size:", file.size(out_path), "bytes\n")
cat("Cause-1 events:", sum(y[, 2] == 1), "\n")
cat("conv:", fit$conv, "\n")
cat("iterations:", fit$iterations, "\n")
