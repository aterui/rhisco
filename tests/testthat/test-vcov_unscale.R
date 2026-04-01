
test_that("vcov_unscale returns correct variance-covariance", {

  ## lm
  m0 <- lm(Sepal.Length ~ Sepal.Width, iris)
  m1 <- lm(Sepal.Length ~ scale(Sepal.Width), iris)

  expect_equal(
    unname(vcov(m0)),
    vcov_unscale(vcov(m1),
                 means = mean(iris$Sepal.Width),
                 stds = sd(iris$Sepal.Width))
  )

  ## lmer
  m0 <- lme4::lmer(Sepal.Length ~ Sepal.Width + (1 | Species),
                   REML = FALSE,
                   iris)
  m1 <- lme4::lmer(Sepal.Length ~ scale(Sepal.Width) + (1 | Species),
                   REML = FALSE,
                   iris)

  expect_equal(
    unname(as.matrix(vcov(m0))),
    vcov_unscale(as.matrix(vcov(m1)),
                 means = mean(iris$Sepal.Width),
                 stds = sd(iris$Sepal.Width))
  )

  ## glmmTMB
  m0 <- glmmTMB::glmmTMB(Sepal.Length ~ Sepal.Width + (1 | Species),
                         REML = FALSE,
                         iris)
  m1 <- glmmTMB::glmmTMB(Sepal.Length ~ scale(Sepal.Width) + (1 | Species),
                         REML = FALSE,
                         iris)

  expect_equal(
    unname(as.matrix(vcov(m0)$cond)),
    vcov_unscale(as.matrix(vcov(m1)$cond),
                 means = mean(iris$Sepal.Width),
                 stds = sd(iris$Sepal.Width)),
    tolerance = 1e-5
  )
})

