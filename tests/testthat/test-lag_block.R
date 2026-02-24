# tests/testthat/test-lag_block.R

library(tidyverse)

test_that("lag_block works without index", {

  n <- sample(100, 1)

  df <- data.frame(
    t = 1:n,
    n = rpois(n, lambda = 10)
  )

  out <- lag_block(df, y = "n", t = "t")

  expect_equal(nrow(out), n)
  expect_true("n_lag" %in% names(out))
  expect_equal(out$n_lag[-1], df$n[-n])
})
