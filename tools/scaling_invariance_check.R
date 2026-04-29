pacman::p_load(lme4, glmmTMB, tidyverse)

scl_iris <- iris %>%
  mutate(across(.cols = c(Sepal.Width, Petal.Width),
                .fns = function(x) (x - mean(x)) / sd(x)))

f <- Sepal.Length ~ Sepal.Width + Petal.Width + (1 + Sepal.Width + Petal.Width | Species)

ctrl <- lmerControl(optCtrl = list(ftol_abs = 1e-12, xtol_abs = 1e-12))
m0 <- lmer(f, iris, REML = F, control = ctrl)
m1 <- lmer(f, scl_iris, REML = F, control = ctrl)

sw0 <- 4
pw0 <- mean(iris$Petal.Width)

sw1 <- (sw0 - mean(iris$Sepal.Width)) / sd(iris$Sepal.Width)
pw1 <- (pw0 - mean(iris$Petal.Width)) / sd(iris$Petal.Width)

const0 <- c("(Intercept)" = 1, "Sepal.Width" = sw0, "Petal.Width" = pw0)
const1 <- c("(Intercept)" = 1, "Sepal.Width" = sw1, "Petal.Width" = pw1)
vf0 <- drop(t(const0) %*% vcov(m0) %*% const0)
vf1 <- drop(t(const1) %*% vcov(m1) %*% const1)

reterm <- rownames(VarCorr(m0)[["Species"]])
rv0 <- drop(t(const0[reterm]) %*% VarCorr(m0)[["Species"]] %*% const0[reterm])
rv1 <- drop(t(const1[reterm]) %*% VarCorr(m1)[["Species"]] %*% const1[reterm])

mu0 <- drop(c(1, sw0, pw0) %*% lme4::fixef(m0))
mu1 <- drop(c(1, sw1, pw1) %*% lme4::fixef(m1))

pnorm(6, mean = mu0, sd = sqrt(vf0 + rv0)) -
pnorm(6, mean = mu1, sd = sqrt(vf1 + rv1))
