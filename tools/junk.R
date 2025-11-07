
library(tidyverse)
library(rhisco)
library(cdyns)

A <- rbind(c(1, 1, 1.5, 0, 0),
           c(1.5, 1, 0, 0, 0),
           c(1.1, 1.5, 1, 0, 0),
           c(0, 0, 0, 1, 0),
           c(0, 0, 0, 0, 1))

# A <- matrix(rexp(25, 2),
#             5, 5)
# diag(A) <- 1

df_lag <- cdynsim(n_species = 5,
                  n_timestep = 20,
                  n_burnin = 0,
                  int_type = "constant",
                  alpha = A,
                  model = "bh",
                  inv_sign = FALSE) %>%
  with(df_dyn) %>%
  rename(y = density) %>%
  rhisco::lag_block(y = "y",
                    t = "timestep",
                    index = "species") %>%
  select(-immigrant) %>%
  mutate(log_r = log(y) - log(y_lag),
         log_rt = log(yt) - log(yt_lag)) %>%
  drop_na(log_r)

x_star <- rhisco::xeq(log_rt ~ yt_lag,
                      data = df_lag)


df_plot <- df_lag %>%
  group_by(species) %>%
  mutate(d = sqrt(y_lag^2 + (yt_lag - x_star)^2),
         max_di = max(d)) %>%
  ungroup() %>%
  mutate(wi = exp(-(d / max_di)),
         w = exp(-(d / max(d))))

df_plot %>%
  ggplot(aes(x = y_lag,
             y = yt_lag,
             color = wi)) +
  geom_point(size = 3) +
  geom_hline(yintercept = x_star) +
  facet_wrap(facets =~ species,
             scales = "free") +
  MetBrewer::scale_color_met_c("Hiroshige",
                               direction = -1)

df_plot %>%
  ggplot(aes(x = y_lag,
             y = log_r,
             color = wi)) +
  geom_point(size = 3) +
  facet_wrap(facets =~ species) +
  MetBrewer::scale_color_met_c("Hiroshige",
                               direction = -1)

(m0 <- lme4::lmer(log_r ~ y_lag + yt_lag + (1 | species),
                  data = df_plot,
                  weights = wi,
                  REML = FALSE))

(m0 <- lm(log_r ~ y_lag + yt_lag,
          data = df_plot,
          weights = wi))

(m1 <- lme4::lmer(log_r ~ y_lag + yt_lag + (1 | species),
                  data = df_plot,
                  weights = w,
                  REML = FALSE))

(m1 <- lm(log_r ~ y_lag + yt_lag,
          data = df_plot,
          weights = w))

# #' @param m Integer. Magnitude of a given link.
# #'
# #' @importFrom stats dpois qpois
# #'
# #' @author Akira Terui, \email{hanabi0111@gmail.com}
# #'
# #' @export

# pdev <- function(fit,
#                  newdata = NULL,
#                  family = NULL,
#                  type = "response") {
#
#   # fit: a fitted model (lm, glm, lmer, glmer)
#   # newdata: optional data.frame for predictions
#   # family: specify "gaussian", "poisson", "binomial"; inferred if NULL
#   # type: type of prediction ("response" by default)
#
#   # Get predicted mean / probability
#   mu_hat <- if (is.null(newdata)) {
#     predict(fit, type = type, re.form = NULL)  # for merMod, includes random effects
#   } else {
#     if ("merMod" %in% class(fit)) {
#       predict(fit, newdata = newdata, type = type, re.form = NULL)
#     } else {
#       predict(fit, newdata = newdata, type = type)
#     }
#   }
#
#   # Observed response
#   if (is.null(newdata)) {
#     y <- model.response(model.frame(fit))
#   } else {
#     y <- newdata[[all.vars(formula(fit))[1]]]
#   }
#
#   # Infer family if not provided
#   if (is.null(family)) {
#     family <- if ("glm" %in% class(fit)) {
#       fit$family$family
#     } else if ("glmerMod" %in% class(fit)) {
#       fit@resp$family$family
#     } else if ("lm" %in% class(fit) || "lmerMod" %in% class(fit)) {
#       "gaussian"
#     } else {
#       stop("Cannot infer family; please specify.")
#     }
#   }
#
#   # Compute pointwise deviance
#   dev <- switch(family,
#                 gaussian = (y - mu_hat)^2,
#                 poisson  = {
#                   y_pos <- y > 0
#                   d <- numeric(length(y))
#                   d[y_pos] <- 2 * (y[y_pos] * log(y[y_pos] / mu_hat[y_pos]) - (y[y_pos] - mu_hat[y_pos]))
#                   d[!y_pos] <- 2 * (- mu_hat[!y_pos])
#                   d
#                 },
#                 binomial = {
#                   # Support binary (n=1) and general (n_i > 1)
#                   if (is.matrix(y) && ncol(y) == 2) {
#                     # y = cbind(successes, failures)
#                     n <- rowSums(y)
#                     s <- y[,1]
#                   } else {
#                     n <- 1
#                     s <- y
#                   }
#                   d <- numeric(length(y))
#                   pos <- s > 0
#                   fail <- (n - s) > 0
#                   d[pos] <- 2 * s[pos] * log(s[pos] / (n[pos] * mu_hat[pos]))
#                   d[fail] <- d[fail] + 2 * (n[fail] - s[fail]) * log((n[fail] - s[fail]) / (n[fail] * (1 - mu_hat[fail])))
#                   d
#                 },
#                 stop("Unsupported family")
#   )
#
#   return(dev)
# }

