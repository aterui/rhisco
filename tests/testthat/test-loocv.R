# tests/testthat/test-loocv.R

library(tidyverse)

test_that("loocv returns numeric RMSE", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5,
    species = 1
  )

  rmse <- loocv(
    y ~ x,
    data = df,
    theta = 1,
    model = "lm",
    type = "gaussian"
  )

  expect_true(is.numeric(rmse))
  expect_length(rmse, 1)
})

test_that("loocv errors for missing group column", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5
  )

  expect_error(
    loocv(
      y ~ x,
      data = df,
      theta = 1,
      model = "lm",
      group = "species"
    ),
    "species is not found in the dataframe"
  )
})

test_that("loocv errors for no time column", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    species = 1
  )

  expect_error(
    loocv(
      y ~ x,
      data = df,
      theta = 1,
      model = "lm"
    ),
    "No time column found"
  )
})

test_that("loocv errors if operational units do not have full timesteps", {
  df <- data.frame(
    y = rnorm(6),
    x = 1:6,
    t = c(1,1,2,2,3,3),
    species = c(1,2,1,2,1,3) # species 3 missing timestep
  )

  expect_error(
    loocv(
      y ~ x,
      data = df,
      theta = 1,
      model = "lm",
      group = "species"
    ),
    "must contain one observation in each timestep"
  )
})

test_that("loocv supports different weight types and function", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5,
    species = 1
  )

  df_mtype <- expand.grid(
    method = c("max", "mean"),
    type   = c("gaussian", "exp")
  ) %>%
    mutate(across(.cols = everything(),
                  .fns = as.character))

  sapply(
    1:nrow(df_mtype),
    function(i) {
      rmse <- loocv(
        y ~ x,
        data = df,
        theta = 1,
        model = "lm",
        type   = df_mtype[i, "type"],
        method = df_mtype[i, "method"]
      )

      expect_true(is.numeric(rmse))
      expect_length(rmse, 1)
    }
  )
})

test_that("loocv seed produces reproducible results", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5,
    species = 1
  )

  rmse1 <- loocv(
    y ~ x,
    data = df,
    theta = 1,
    model = "lm",
    seed = 123
  )

  rmse2 <- loocv(
    y ~ x,
    data = df,
    theta = 1,
    model = "lm",
    seed = 123
  )

  expect_equal(rmse1, rmse2)
})

test_that("loocv errors for unsupported model type", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5,
    species = 1
  )

  expect_error(
    loocv(
      y ~ x,
      data = df,
      theta = 1,
      model = "unsupported"
    ),
    "Unsupported model type"
  )
})

# --- Non-lm model tests and random-effects formula handling ---

test_that("loocv works with glm", {
  skip_if_not_installed("stats")

  df <- data.frame(
    y = rpois(5, lambda = 3),
    x = 1:5,
    t = 1:5,
    species = 1
  )

  rmse <- loocv(
    y ~ x,
    data = df,
    theta = 1,
    model = "glm",
    type = "gaussian"
  )

  expect_true(is.numeric(rmse))
  expect_length(rmse, 1)
})

test_that("loocv errors if random effects used in lm/glm", {
  df <- data.frame(
    y = rnorm(5),
    x = 1:5,
    t = 1:5,
    g = 1
  )

  expect_error(
    loocv(
      y ~ x | g,
      data = df,
      theta = 1,
      model = "lm"
    ),
    "does not allow random effects"
  )

  expect_error(
    loocv(
      y ~ x | g,
      data = df,
      theta = 1,
      model = "glm"
    ),
    "does not allow random effects"
  )
})

test_that("loocv works with lmer if lme4 installed", {
  skip_if_not_installed("lme4")

  df <- data.frame(
    t = 1:5
  ) %>%
    tidyr::crossing(
      species = 1:5
    ) %>%
    mutate(
      x = rnorm(n = nrow(.)),
      y = rnorm(n = nrow(.))
    )

  rmse <- loocv(
    y ~ x + (1|species),
    data = df,
    theta = 1,
    model = "lmer",
    group = "species"
  )

  expect_true(is.numeric(rmse))
  expect_length(rmse, 1)
})

test_that("loocv works with glmer if lme4 installed", {
  skip_if_not_installed("lme4")

  df <- data.frame(
    t = 1:5
  ) %>%
    tidyr::crossing(species = 1:5) %>%
    mutate(
      x = rnorm(n = nrow(.)),
      y = rpois(n = nrow(.), lambda = 5)
    )

  rmse <- loocv(
    y ~ x + (1|species),
    data = df,
    theta = 1,
    model = "glmer",
    family = "poisson",
    group = "species"
  )

  expect_true(is.numeric(rmse))
  expect_length(rmse, 1)
})

test_that("loocv works with glmmTMB if glmmTMB installed", {
  skip_if_not_installed("glmmTMB")

  df <- data.frame(
    t = 1:5
  ) %>%
    tidyr::crossing(species = 1:5) %>%
    mutate(
      x = rnorm(n = nrow(.)),
      y = rpois(n = nrow(.), lambda = 5)
    )

  rmse <- loocv(
    y ~ x + (1|species),
    data = df,
    theta = 1,
    model = "glmmTMB",
    family = "poisson",
    group = "species"
  )

  expect_true(is.numeric(rmse))
  expect_length(rmse, 1)
})
