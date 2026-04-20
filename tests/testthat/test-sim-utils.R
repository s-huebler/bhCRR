test_that("%||% returns left side when non-NULL", {
  expect_equal("a" %||% "b", "a")
  expect_equal(0L  %||% 99,  0L)
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that("%||% returns right side when left is NULL", {
  expect_equal(NULL %||% "b", "b")
  expect_equal(NULL %||% 42,  42)
})
