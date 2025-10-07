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

#' @param x Vector for embedding.
#' @param k Integer scalar specifying a lag.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

lag_base <- function(x, k = 1L) {
  if (k >= 0) c(rep(NA, k), head(x, -k))
  else c(tail(x, k), rep(NA, -k))
}

