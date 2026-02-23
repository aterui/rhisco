# tests/testthat/test-set_coef.R

test_that("set_coef returns s x s matrix", {
  s <- sample(2:10, 1)
  alpha <- par.control(dist = "constant", mu = 1, seed = 123)
  beta <- par.control(dist = "constant", mu = 2, seed = 123)

  m <- set_coef(s, alpha = alpha, beta = beta)

  expect_true(is.matrix(m))
  expect_equal(dim(m), c(s, s))
})

test_that("set_coef errors if alpha or beta not list", {
  s <- sample(2:10, 1)
  expect_error(set_coef(s, alpha = "wrong"), "'alpha' must be a list")
  expect_error(set_coef(s, beta = "wrong"), "'beta' must be a list")
})

test_that("set_coef negative argument works", {
  s <- sample(2:5, 1)
  alpha <- par.control(dist = "constant", mu = 1, seed = 123)
  beta <- par.control(dist = "constant", mu = 2, seed = 123)

  m_neg <- set_coef(s, alpha = alpha, beta = beta, negative = TRUE)
  m_pos <- set_coef(s, alpha = alpha, beta = beta, negative = FALSE)

  expect_equal(m_neg, -m_pos)
})

test_that("set_coef supports scalar and vector mu for alpha", {
  s <- sample(2:5, 1)
  # scalar mu
  alpha <- par.control(dist = "constant", mu = 5, seed = 123)
  beta <- par.control(dist = "constant", mu = 2, seed = 123)
  m <- set_coef(s, alpha = alpha, beta = beta)
  expect_equal(diag(m), -rep(5, s))

  # vector mu
  alpha <- par.control(dist = "constant", mu = seq_len(s), seed = 123)
  m <- set_coef(s, alpha = alpha, beta = beta)
  expect_equal(diag(m), -seq_len(s))
})

test_that("set_coef errors if alpha vector length mismatch", {
  s <- sample(3:5, 1)
  alpha <- par.control(dist = "constant", mu = seq_len(s - 1), seed = 123)
  beta <- par.control(dist = "constant", mu = 2, seed = 123)
  expect_error(set_coef(s, alpha = alpha, beta = beta),
               "'mu' in 'alpha' must be scalar or vector of length")
})

test_that("set_coef supports different dist types", {
  s <- sample(2:5, 1)
  dist_types <- c("constant", "exp", "unif", "normal")
  for (d in dist_types) {
    alpha <- par.control(dist = d, mu = 1, min = 0, max = 2, sd = 0.1, seed = 123)
    beta <- par.control(dist = d, mu = 2, min = 0, max = 2, sd = 0.2, seed = 123)
    m <- set_coef(s, alpha = alpha, beta = beta)
    expect_equal(dim(m), c(s, s))
  }
})

test_that("set_coef errors for unsupported dist in alpha/beta", {
  s <- sample(2:5, 1)
  alpha <- list(dist = "unsupported", mu = 1, seed = 123)
  beta <- list(dist = "constant", mu = 2, seed = 123)
  expect_error(set_coef(s, alpha = alpha, beta = beta), "Unsupported dist in 'alpha'")

  alpha <- list(dist = "constant", mu = 1, seed = 123)
  beta <- list(dist = "unsupported", mu = 2, seed = 123)
  expect_error(set_coef(s, alpha = alpha, beta = beta), "Unsupported dist in 'beta'")
})
