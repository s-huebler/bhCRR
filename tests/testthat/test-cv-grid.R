# Tests for R/cv_grid.r — .cv_grid()
# Internal function; visible here because tests run inside the package namespace.

test_that("a 5x4 grid with all valid pairs returns 20 rows", {
  # s0_seq and s1_seq chosen so every s1 > every s0 (all 20 pairs valid)
  s0 <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  s1 <- c(0.5,  0.6,  0.7,  0.8)
  g  <- .cv_grid(s0, s1)
  expect_equal(nrow(g), 20L)
  expect_equal(ncol(g), 3L)
  expect_named(g, c("s0", "s1", "pair"))
  expect_equal(g$pair, seq_len(20L))
})

test_that("pairs with s1 <= s0 are dropped", {
  s0 <- c(0.1, 0.3, 0.5)
  s1 <- c(0.2, 0.4)
  # Valid pairs: (0.1,0.2),(0.1,0.4),(0.3,0.4) — dropped: (0.3,0.2),(0.5,0.2),(0.5,0.4)
  g <- .cv_grid(s0, s1)
  expect_true(all(g$s1 > g$s0))
  expect_equal(nrow(g), 3L)
})

test_that("traversal order matches tune_ssl_psdh logic exactly", {
  s0_seq <- c(0.02, 0.05, 0.01)   # deliberate non-sorted s0
  s1_seq <- c(0.5,  0.3,  0.4)    # deliberate non-sorted s1

  g <- .cv_grid(s0_seq, s1_seq)

  # Replicate tune_ssl_psdh's order derivation directly:
  #   expand.grid -> filter -> unique(s1) in grid order -> sort s0 per group
  ref_grid  <- expand.grid(s0 = s0_seq, s1 = s1_seq, stringsAsFactors = FALSE)
  ref_valid <- ref_grid[ref_grid$s1 > ref_grid$s0, , drop = FALSE]
  ref_s1    <- unique(ref_valid$s1)   # s1 in traversal order
  ref_rows  <- do.call(rbind, lapply(ref_s1, function(s1v) {
    data.frame(s0 = sort(ref_valid$s0[ref_valid$s1 == s1v]), s1 = s1v)
  }))

  expect_equal(g$s0, ref_rows$s0)
  expect_equal(g$s1, ref_rows$s1)
})

test_that("no valid pairs stops with an informative message", {
  # All s1 <= s0
  expect_error(
    .cv_grid(s0_seq = c(0.5, 0.8), s1_seq = c(0.1, 0.3)),
    regexp = "No valid"
  )
  # Error names the ranges
  expect_error(
    .cv_grid(s0_seq = c(0.5, 0.8), s1_seq = c(0.1, 0.3)),
    regexp = "0\\.5.*0\\.8|0\\.1.*0\\.3"
  )
})

test_that("pair column is a contiguous 1..npairs integer sequence", {
  g <- .cv_grid(c(0.01, 0.02), c(0.1, 0.2, 0.3))
  expect_equal(g$pair, seq_len(nrow(g)))
})
