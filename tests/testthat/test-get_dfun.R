# tests/testthat/test-get_dfun.R

library(tidyverse)

test_that("exp type works with max scaling", {

  x <- c(1, 2, 4)
  theta <- 1

  f <- get_dfun("exp", method = "max")

  out <- f(x, theta)

  # manual calculation
  d <- x / max(x)
  expected <- exp(-theta * d)

  expect_equal(out, expected)
})

test_that("gaussian type works with max scaling", {

  x <- c(1, 2, 4)
  theta <- 2

  f <- get_dfun("gaussian", method = "max")

  out <- f(x, theta)

  d <- x / max(x)
  expected <- exp(- theta^2 * (d^2 / 2))

  expect_equal(out, expected)
})

test_that("unsupported type throws error", {

  expect_error(
    get_dfun("logistic"),
    "Unsupported type"
  )
})
