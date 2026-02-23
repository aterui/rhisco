# tests/testthat/test-set_r.R

test_that("set_r returns vector of length s", {
  s <- sample(2:10, 1)
  r <- par.control(dist = "constant", mu = 1, seed = 123)

  v <- set_r(s, r = r)

  expect_true(is.numeric(v))
  expect_equal(length(v), s)
})

test_that("set_r errors if r is not a list", {
  s <- sample(2:10, 1)
  expect_error(set_r(s, r = "wrong"), "'r' must be a list")
})

test_that("set_r supports scalar and vector mu", {
  s <- sample(2:5, 1)

  # scalar mu
  r <- par.control(dist = "constant", mu = 5, seed = 123)
  v <- set_r(s, r = r)
  expect_equal(v, rep(5, s))

  # vector mu
  r <- par.control(dist = "constant", mu = seq_len(s), seed = 123)
  v <- set_r(s, r = r)
  expect_equal(v, seq_len(s))
})

test_that("set_r errors if vector mu length mismatch", {
  s <- sample(3:5, 1)
  r <- list(dist = "constant", mu = seq_len(s - 1), seed = 123)
  expect_error(set_r(s, r = r), "'mu' in 'r' must be scalar or vector of length")
})

test_that("set_r supports different dist types", {
  s <- sample(2:5, 1)
  dist_types <- c("constant", "exp", "unif", "normal")

  for (d in dist_types) {
    r <- par.control(dist = d, mu = 1, min = 0, max = 2, sd = 0.1, seed = 123)
    v <- set_r(s, r = r)
    expect_equal(length(v), s)
  }
})

test_that("set_r errors for unsupported dist", {
  s <- sample(2:5, 1)
  r <- list(dist = "unsupported", mu = 1, seed = 123)
  expect_error(set_r(s, r = r), "Unsupported dist in 'r'")
})

test_that("set_r is reproducible with seed", {
  s <- sample(2:5, 1)
  r <- par.control(dist = "unif", min = 0, max = 1, seed = 123)

  v1 <- set_r(s, r = r)
  v2 <- set_r(s, r = r)

  expect_equal(v1, v2)
})
