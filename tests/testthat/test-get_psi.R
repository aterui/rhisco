# tests/testthat/test-get_psi.R

library(tidyverse)

test_that("get_psi returns numeric psi with correct attributes", {

  df <- csim(ts = 10) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag)) %>%
    drop_na(y)

  psi <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    data = df,
    x_star = 15,
    theta = seq(0.5, 2, by = 0.5),
    n_sim = 50,
    model = "lmer"
  )

  expect_true(is.numeric(psi))
  expect_true(psi >= 0 && psi <= 1)
  expect_true(!is.null(attr(psi, "theta")))
  expect_true(!is.null(attr(psi, "rmse")))
  expect_true(!is.null(attr(psi, "message")))
})

test_that("get_psi works with glm", {

  df <- csim(ts = 10) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag)) %>%
    drop_na(y)

  psi <- get_psi(
    y ~ n_lag + nt_lag,
    group = "species",
    data = df,
    x_star = 15,
    model = "glm",
    type = "gaussian",
    theta = seq(0.5, 2, by = 0.5),
    n_sim = 50
  )

  expect_true(is.numeric(psi))
})

test_that("get_psi works with lmer/glmer if lme4 installed", {

  df <- csim(ts = 10) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag)) %>%
    drop_na(y)

  psi <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    data = df,
    x_star = 0,
    model = "lmer",
    n_sim = 10
  )

  expect_true(is.numeric(psi))

  psi_glmer <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    data = df,
    x_star = 0,
    model = "glmer",
    group = "species",
    n_sim = 10
  )

  expect_true(is.numeric(psi_glmer))
})

test_that("get_psi works with glmmTMB if installed", {
  skip_if_not_installed("glmmTMB")

  df <- data.frame(
    y = rpois(10, lambda = 2),
    n_lag = 1:10,
    nt_lag = 11:20,
    species = 1
  )

  psi <- get_psi(
    y ~ n_lag + nt_lag,
    data = df,
    x_star = 15,
    model = "glmmTMB",
    n_sim = 10
  )

  expect_true(is.numeric(psi))
})
