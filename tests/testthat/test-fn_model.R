# tests/testthat/test-fn_model.R

test_that("fn_model returns a function", {
  dyn_bh <- fn_model("bh")
  dyn_ricker <- fn_model("ricker")

  expect_true(is.function(dyn_bh))
  expect_true(is.function(dyn_ricker))
})

test_that("fn_model errors for unsupported model", {
  expect_error(fn_model("unsupported"), "Unsupported model type")
})

test_that("fn_model returns correct length output", {
  s <- 3
  n <- rep(10, s)
  r <- rep(0.1, s)
  m_coef <- diag(-0.05, s, s)
  eps <- rep(0, s)

  dyn_bh <- fn_model("bh")
  dyn_ricker <- fn_model("ricker")

  y_bh <- dyn_bh(r, n, m_coef, eps)
  y_ricker <- dyn_ricker(r, n, m_coef, eps)

  expect_equal(length(y_bh), s)
  expect_equal(length(y_ricker), s)
})

test_that("fn_model stochastic argument produces integer counts", {
  s <- 3
  n <- rep(10, s)
  r <- rep(0.1, s)
  m_coef <- diag(-0.05, s, s)
  eps <- rep(0, s)

  dyn_bh <- fn_model("bh")
  dyn_ricker <- fn_model("ricker")

  y_bh_stoch <- dyn_bh(r, n, m_coef, eps, stochastic = TRUE)
  y_ricker_stoch <- dyn_ricker(r, n, m_coef, eps, stochastic = TRUE)

  expect_true(all(y_bh_stoch >= 0))
  expect_true(all(y_ricker_stoch >= 0))
  expect_true(all(y_bh_stoch == floor(y_bh_stoch)))
  expect_true(all(y_ricker_stoch == floor(y_ricker_stoch)))
})

test_that("fn_model non-stochastic produces numeric output", {
  s <- 3
  n <- rep(10, s)
  r <- rep(0.1, s)
  m_coef <- diag(-0.05, s, s)
  eps <- rep(0, s)

  dyn_bh <- fn_model("bh")
  y <- dyn_bh(r, n, m_coef, eps, stochastic = FALSE)

  expect_true(is.numeric(y))
  expect_equal(length(y), s)
})
