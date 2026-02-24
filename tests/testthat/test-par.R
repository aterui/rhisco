# tests/testthat/test-par.R

library(tidyverse)

test_that("par.control returns a list with expected names", {
  out <- par.control()
  expect_type(out, "list")
  expect_named(out, c("dist", "pg", "mu", "sd", "min", "max", "intv", "seed"))
})

test_that("par.control default values are correct", {
  out <- par.control()
  expect_equal(out$dist, "constant")
  expect_equal(out$pg, "sync")
  expect_equal(out$mu, 1)
  expect_equal(out$sd, 0.1)
  expect_equal(out$min, 0)
  expect_equal(out$max, 1)
  expect_equal(out$intv, 1)
  expect_true(is.na(out$seed))
})

test_that("par.control matches valid dist and pg arguments", {
  # Test all valid dist options
  dist_options <- c("constant", "exp", "unif", "normal")
  for (d in dist_options) {
    out <- par.control(dist = d)
    expect_equal(out$dist, d)
  }

  # Test all valid pg options
  pg_options <- c("sync", "ord")
  for (p in pg_options) {
    out <- par.control(pg = p)
    expect_equal(out$pg, p)
  }
})

test_that("par.control errors for invalid dist or pg arguments", {
  expect_error(par.control(dist = "invalid"), "should be one of")
  expect_error(par.control(pg = "invalid"), "should be one of")
})

test_that("par.control accepts custom numeric parameters", {
  out <- par.control(mu = 5, sd = 0.5, min = 2, max = 10, intv = 3, seed = 123)
  expect_equal(out$mu, 5)
  expect_equal(out$sd, 0.5)
  expect_equal(out$min, 2)
  expect_equal(out$max, 10)
  expect_equal(out$intv, 3)
  expect_equal(out$seed, 123)
})
