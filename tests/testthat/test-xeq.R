# tests/testthat/test-xeq.R

library(tidyverse)

test_that("xeq supports different weight types and functions", {

  df <- data.frame(
    x = runif(10),
    t = 1:10
  ) %>%
    mutate(
      y = rnorm(10, mean = 1 - 2 * x, sd = 0.1)
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
      x_star <- xeq(
        y ~ x,
        data   = df,
        theta  = seq(0.5, 2, by = 0.5),
        type   = df_mtype[i, "type"],
        method = df_mtype[i, "method"]
      )

      expect_true(is.numeric(x_star))
      expect_true(attr(x_star, "gap") < 1)
      expect_true(attr(x_star, "iteration") <= 50)
      expect_true(attr(x_star, "theta") %in% seq(0.5, 2, by = 0.5))
      expect_length(attr(x_star, "rmse"), length(seq(0.5, 2, by = 0.5)))
    }
  )
})

test_that("xeq supports different weight types and functions", {

  df <- data.frame(
    x = runif(10),
    t = 1:10
  ) %>%
    mutate(
      y = rnorm(10, mean = 1 + 2 * x, sd = 0.1)
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
      expect_message(
        xeq(
          y ~ x,
          data   = df,
          theta  = seq(0.5, 2, by = 0.5),
          type   = df_mtype[i, "type"],
          method = df_mtype[i, "method"]
        ),
        "No feasible equilibrium exists."
      )
    })

})
