# tests/testthat/test-csim.R

library(tidyverse)

test_that("csim returns tibble with correct columns and rows", {
  s <- sample(2:5, 1)
  ts <- 10
  warmup <- 5
  burnin <- 5

  df <- csim(ts = ts,
             warmup = warmup,
             burnin = burnin,
             s = s,
             r = par.control(dist = "constant", mu = 0.1, seed = 123),
             alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
             beta = par.control(dist = "constant", mu = 0.1, seed = 123),
             model = "ricker",
             sd_eps = 0,
             stochastic = FALSE,
             pg = par.control(pg = "sync", min = 0, max = 1, seed = 123),
             extinct = 0,
             trim = TRUE)

  # Check type and column names
  expect_s3_class(df, "tbl_df")
  expect_equal(colnames(df), c("ts", "species", "density"))

  # Check number of rows (after trimming)
  n_sim <- warmup + burnin + ts
  n_cut <- warmup + burnin
  expect_equal(nrow(df), ts * s)
})

test_that("csim preserves attributes r and coef", {
  s <- 3
  df <- csim(ts = 5, warmup = 2, burnin = 2, s = s,
             r = par.control(dist = "constant", mu = 0.1, seed = 123),
             alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
             beta = par.control(dist = "constant", mu = 0.1, seed = 123),
             model = "bh",
             trim = TRUE)

  expect_true(!is.null(attr(df, "r")))
  expect_true(!is.null(attr(df, "coef")))
  expect_equal(length(attr(df, "r")), s)
  expect_equal(dim(attr(df, "coef")), c(s, s))
})

test_that("csim enforces extinction threshold", {
  s <- 3

  df_det <- csim(ts = 5, warmup = 2, burnin = 2, s = s,
                 r = par.control(dist = "constant", mu = -5, seed = 123),
                 alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
                 beta = par.control(dist = "constant", mu = 0.1, seed = 123),
                 model = "bh",
                 sd_eps = 0,
                 stochastic = FALSE,
                 extinct = 1,
                 trim = TRUE)

  df_stoch <- csim(ts = 5, warmup = 2, burnin = 2, s = s,
                   r = par.control(dist = "constant", mu = -5, seed = 123),
                   alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
                   beta = par.control(dist = "constant", mu = 0.1, seed = 123),
                   model = "bh",
                   sd_eps = 0,
                   stochastic = TRUE,
                   extinct = 1,
                   trim = TRUE)

  expect_true(all(c(df_det$density, df_stoch$density) == 0))
})

test_that("csim trim = FALSE returns full matrix", {
  s <- 3
  ts <- 5
  warmup <- 2
  burnin <- 2

  df <- csim(ts = ts, warmup = warmup, burnin = burnin, s = s,
             r = par.control(dist = "constant", mu = 0.1, seed = 123),
             alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
             beta = par.control(dist = "constant", mu = 0.1, seed = 123),
             model = "ricker",
             trim = FALSE)

  n_sim <- ts + warmup + burnin
  expect_equal(nrow(df), n_sim * s)
})

test_that("csim stochastic produces numeric output", {
  s <- 3
  df <- csim(ts = 5, warmup = 2, burnin = 2, s = s,
             r = par.control(dist = "constant", mu = 0.1, seed = 123),
             alpha = par.control(dist = "constant", mu = 0.1, seed = 123),
             beta = par.control(dist = "constant", mu = 0.1, seed = 123),
             model = "bh",
             stochastic = TRUE)

  expect_true(all(df$density >= 0))
})
