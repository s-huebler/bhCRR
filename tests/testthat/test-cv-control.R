# Tests for R/cv_control.r

test_that("bhcrr_cv_control() defaults construct without error", {
  ctrl <- bhcrr_cv_control()
  expect_s3_class(ctrl, "bhcrr_cv_control")
  expect_equal(ctrl$nfolds, 10L)
  expect_equal(ctrl$ncv,    2L)
  expect_equal(ctrl$strata, "cause1")
  expect_equal(ctrl$init_method, "LASSO_cv")
  expect_equal(ctrl$init_label,  "LASSO_cv")
  expect_false(ctrl$parallel)
})

test_that("print.bhcrr_cv_control() runs without error and returns invisibly", {
  ctrl <- bhcrr_cv_control()
  out  <- capture.output(ret <- print(ctrl))
  expect_identical(ret, ctrl)
  expect_true(any(grepl("bhcrr_cv_control", out)))
  expect_true(any(grepl("resampling", out)))
  expect_true(any(grepl("init", out)))
})

test_that("strata is matched by match.arg()", {
  expect_equal(bhcrr_cv_control(strata = "cause1")$strata, "cause1")
  expect_equal(bhcrr_cv_control(strata = "none")$strata,   "none")
  expect_equal(bhcrr_cv_control(strata = "status")$strata, "status")
  # partial match
  expect_equal(bhcrr_cv_control(strata = "cau")$strata, "cause1")
})

test_that("unknown strata value errors", {
  expect_error(bhcrr_cv_control(strata = "random"), regexp = "arg")
})

test_that("a bad init_method errors at construction, not later", {
  expect_error(
    bhcrr_cv_control(init_method = "not_a_real_method"),
    regexp = "LASSO_cv"   # .resolve_init_method lists the built-in names
  )
})

test_that("fit_args containing 'init' errors with a message naming init_method", {
  expect_error(
    bhcrr_cv_control(fit_args = list(init = rep(0, 5))),
    regexp = "init_method"
  )
})

test_that("fit_args containing 'init_method' errors", {
  expect_error(
    bhcrr_cv_control(fit_args = list(init_method = "zero")),
    regexp = "init_method"
  )
})

test_that("fit_args containing 'x' or 'y' errors", {
  expect_error(bhcrr_cv_control(fit_args = list(x = matrix(1))), regexp = "'x'|x")
  expect_error(bhcrr_cv_control(fit_args = list(y = matrix(1))), regexp = "'y'|y")
})

test_that("fit_args containing a non-formal of fit_ssl_psdh errors and names it", {
  expect_error(
    bhcrr_cv_control(fit_args = list(not_a_real_arg = 99)),
    regexp = "not_a_real_arg"
  )
})

test_that("valid fit_args names are accepted", {
  ctrl <- bhcrr_cv_control(fit_args = list(maxit = 100, epsilon = 1e-5,
                                           theta_a = 1, theta_b = 20,
                                           initial_sparsity = 0.1,
                                           inner_maxit_start = 500))
  expect_equal(ctrl$fit_args$maxit, 100)
})

test_that("pool = 0.1 with strata = 'cause1' warns about rsample pooling", {
  expect_warning(
    bhcrr_cv_control(pool = 0.1, strata = "cause1"),
    regexp = "pool|stratum|10"
  )
})

test_that("pool = 0.1 with strata = 'none' does NOT warn", {
  expect_no_warning(bhcrr_cv_control(pool = 0.1, strata = "none"))
})

test_that("nfolds < 2 errors", {
  expect_error(bhcrr_cv_control(nfolds = 1L), regexp = "nfolds")
})

test_that("ncv < 1 errors", {
  expect_error(bhcrr_cv_control(ncv = 0L), regexp = "ncv")
})

test_that("eval_quantile outside (0,1) errors", {
  expect_error(bhcrr_cv_control(eval_quantile = 0),  regexp = "eval_quantile")
  expect_error(bhcrr_cv_control(eval_quantile = 1),  regexp = "eval_quantile")
  expect_error(bhcrr_cv_control(eval_quantile = -1), regexp = "eval_quantile")
})

test_that("pool outside (0, 0.5) errors", {
  expect_error(bhcrr_cv_control(pool = 0),   regexp = "pool")
  expect_error(bhcrr_cv_control(pool = 0.5), regexp = "pool")
})

test_that("foldid is validated as positive-integer matrix", {
  good <- matrix(c(1L, 2L, 1L, 2L), ncol = 1)
  ctrl <- bhcrr_cv_control(foldid = good)
  expect_true(!is.null(ctrl$foldid))

  expect_error(bhcrr_cv_control(foldid = matrix(c(0L, 1L), ncol = 1)),
               regexp = "positive")
})

test_that("a function passed as init_method resolves with label 'custom'", {
  my_fn <- function(x, y, ...) list(init = rep(0, ncol(x)), meta = NULL)
  ctrl  <- bhcrr_cv_control(init_method = my_fn)
  expect_equal(ctrl$init_label, "custom")
})

test_that("seed and workers are coerced to integer", {
  ctrl <- bhcrr_cv_control(seed = 42, parallel = TRUE, workers = 4)
  expect_identical(ctrl$seed,    42L)
  expect_identical(ctrl$workers, 4L)
})

# ---------------------------------------------------------------------------
# eval_time field
# ---------------------------------------------------------------------------

test_that("eval_time accepts a positive scalar and is stored", {
  ctrl <- bhcrr_cv_control(eval_time = 12.5)
  expect_equal(ctrl$eval_time, 12.5)
})

test_that("eval_time = NULL (default) is stored as NULL", {
  ctrl <- bhcrr_cv_control()
  expect_null(ctrl$eval_time)
})

test_that("eval_time = 0 errors", {
  expect_error(bhcrr_cv_control(eval_time = 0), regexp = "eval_time")
})

test_that("eval_time negative errors", {
  expect_error(bhcrr_cv_control(eval_time = -1), regexp = "eval_time")
})

test_that("eval_time = NA errors", {
  expect_error(bhcrr_cv_control(eval_time = NA_real_), regexp = "eval_time")
})

test_that("eval_time = Inf errors", {
  expect_error(bhcrr_cv_control(eval_time = Inf), regexp = "eval_time")
})

test_that("eval_time as a length-2 vector errors", {
  expect_error(bhcrr_cv_control(eval_time = c(1, 2)), regexp = "eval_time")
})

test_that("eval_time as character errors", {
  expect_error(bhcrr_cv_control(eval_time = "12"), regexp = "eval_time")
})

# ---------------------------------------------------------------------------
# keep_coefs field
# ---------------------------------------------------------------------------

test_that("keep_coefs defaults to FALSE", {
  ctrl <- bhcrr_cv_control()
  expect_false(ctrl$keep_coefs)
})

test_that("keep_coefs = TRUE is stored", {
  ctrl <- bhcrr_cv_control(keep_coefs = TRUE)
  expect_true(ctrl$keep_coefs)
})

test_that("keep_coefs = 1L (non-logical) errors", {
  expect_error(bhcrr_cv_control(keep_coefs = 1L), regexp = "keep_coefs")
})

test_that("keep_coefs = NA errors", {
  expect_error(bhcrr_cv_control(keep_coefs = NA), regexp = "keep_coefs")
})

# ---------------------------------------------------------------------------
# print shows horizon correctly
# ---------------------------------------------------------------------------

test_that("print shows eval_time when supplied", {
  ctrl <- bhcrr_cv_control(eval_time = 365)
  out  <- capture.output(print(ctrl))
  expect_true(any(grepl("365", out)))
  expect_true(any(grepl("fixed|eval_time", out)))
})

test_that("print shows 'derived' when eval_time is NULL", {
  ctrl <- bhcrr_cv_control()
  out  <- capture.output(print(ctrl))
  expect_true(any(grepl("derived", out)))
})
