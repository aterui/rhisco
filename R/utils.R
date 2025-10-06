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


#' @param formula Model formula.
#' @param data Data frame.
#' @param theta Numeric vector. Values examined with leave-one-out cross validation.
#' @param model Character. Model type, either \code{lm}, \code{glm}, \code{lmer}, or \code{glmer}.
#' @param ... Additional arguments passed to model functions.
#'
#' @importFrom stats dpois qpois
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

loocv <- function(formula,
                  data,
                  theta,
                  model = "lm",
                  ...) {

  ## get matrix X
  y <- all.vars(formula)[1]
  v_cnm_x <- attr(terms(formula), "term.labels")
  v_cnm_x <- v_cnm_x[!grepl("\\|", v_cnm_x)]

  m_x <- data[, v_cnm_x] |>
    data.matrix()

  ## distance matrix for weighting
  m_dist <- dist(m_x,
                 diag = TRUE,
                 upper = TRUE) |>
    data.matrix()

  v <- seq_len(nrow(data))

  ss <- sapply(v,
               function(i) {

                 ## training data
                 df_train <- data[-i, ]

                 ## calculate weights
                 mu_d <- mean(m_dist[i, -i])
                 df_train$w <- exp(- theta * m_dist[i, -i] / mu_d)

                 lw <- switch(model,
                              lm = lm(formula,
                                      data = df_train,
                                      weights = w,
                                      ...),
                              glm = glm(formula,
                                        data = df_train,
                                        weights = w,
                                        ...),
                              lmer = lme4::lmer(formula,
                                                data = df_train,
                                                weights = w,
                                                ...),
                              glmer = lme4::glmer(formula,
                                                  data = df_train,
                                                  weights = w,
                                                  ...),
                              stop("Unsupported model type")
                 )

                 y0 <- predict(lw, newdata = data[i, , drop = FALSE])
                 y1 <- data[i, y]
                 eps <- (y1 - y0)^2

                 return(eps)
               }) |>
    sum()

  return(sqrt(ss / nrow(m_x)))
}

