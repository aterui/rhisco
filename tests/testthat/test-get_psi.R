# tests/testthat/test-get_psi.R

library(tidyverse)

test_that("get_psi returns numeric psi with correct attributes", {

  df <- csim(ts = 10,
             s = 30,
             r = par.control(dist = "unif",
                             min = 0.5,
                             max = 0.5),
             alpha = par.control(dist = "constant",
                                 mu = 1),
             beta = par.control(dist = "unif",
                                min = 1,
                                max = 1)) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag))

  df_i <- df %>%
    drop_na(y)

  df_t <- df %>%
    distinct(ts, nt, nt_lag) %>%
    mutate(y = log(nt / nt_lag)) %>%
    drop_na(y)

  x_star <- xeq(y ~ nt_lag, data = df_t, theta = 0)

  (psi <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    group = "species",
    data = df_i,
    x_star = x_star,
    model = "lmer",
    rescale = TRUE,
    theta = 0
  ))

  (psi0 <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    group = "species",
    data = df_i,
    x_star = x_star,
    model = "lmer",
    rescale = FALSE,
    theta = 0
  ))

  expect_true(is.numeric(psi))
  expect_true(psi >= 0 && psi <= 1)
  expect_equal(c(psi), c(psi0))
  expect_true(!is.null(attr(psi, "theta")))
  expect_true(!is.null(attr(psi, "rmse")))
  expect_true(!is.null(attr(psi, "message")))
})

test_that("get_psi works with glm", {

  df <- csim(ts = 10,
             s = 6,
             beta = par.control(dist = "unif",
                                min = 1,
                                max = 1)) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag))

  df_i <- df %>%
    drop_na(y)

  df_t <- df %>%
    distinct(ts, nt, nt_lag) %>%
    mutate(y = log(nt / nt_lag)) %>%
    drop_na(y)

  x_star <- xeq(y ~ nt_lag, data = df_t, theta = 0)

  (psi <- get_psi(
    y ~ n_lag + nt_lag,
    group = "species",
    data = df_i,
    x_star = x_star,
    model = "glm",
    type = "gaussian"
  ))

  expect_true(is.numeric(psi))
})

test_that("get_psi works with lmer/glmer if lme4 installed", {

  df <- csim(ts = 10,
             s = 5) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag)) %>%
    drop_na(y)

  psi <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    data = df,
    x_star = 0,
    model = "lmer"
  )

  expect_true(is.numeric(psi))

})

test_that("get_psi works with glmmTMB if installed", {

  df <- csim(ts = 10,
             s = 6) %>%
    rename(n = density) %>%
    lag_block(t = "ts",
              index = "species") %>%
    mutate(y = log(n / n_lag))

  df_i <- df %>%
    drop_na(y)

  df_t <- df %>%
    distinct(ts, nt, nt_lag) %>%
    mutate(y = log(nt / nt_lag)) %>%
    drop_na(y)

  x_star <- xeq(y ~ nt_lag, data = df_t, theta = 0)

  (psi <- get_psi(
    y ~ n_lag + nt_lag + (1 | species),
    data = df_i,
    x_star = x_star,
    model = "glmmTMB"
  ))

  expect_true(is.numeric(psi))
})
