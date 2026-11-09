test_that("null model (all-zero init, tight spike) converges with all-zero coefficients", {
  dat <- make_ssl_fixture_data()
  x <- dat$x
  y <- dat$y

  expect_no_warning(
    fit <- fit_ssl_psdh(x, y,
      ss = c(1e-4, 0.5),
      init = rep(0, ncol(x)))
  )

  expect_true(fit$conv)
  expect_true(all(fit$coefficients$Estimate == 0))
})

test_that("non-converging fit warns and warning mentions support changes", {
  dat <- make_ssl_fixture_data()
  # maxit = 3 guarantees non-convergence: the deviance check requires iter > 5
  expect_warning(
    fit_ssl_psdh(dat$x, dat$y,
      ss = c(0.04, 0.5),
      maxit = 3,
      init_method = "zero"),
    regexp = "support changes"
  )
})
