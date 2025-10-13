#' Leave-One-Out Cross-Validation for Distance-Weighted Local Regression
#'
#' Performs leave-one-out cross-validation (LOOCV) to evaluate localized regression models
#' that apply distance-based weighting to observations. For each candidate value of the
#' weighting parameter \code{theta}, the function computes the root mean squared error (RMSE)
#' across all held-out samples.
#'
#' @param formula A model formula specifying the response and predictor variables.
#' @param data A data frame containing the variables referenced in \code{formula}.
#' @param theta Numeric scalar. A candidate value of the distance-weighting parameter to be evaluated via LOOCV.
#' @param model Character string specifying the model type to fit. Must be one of
#'   \code{"lm"}, \code{"glm"}, \code{"lmer"}, \code{"glmer"}, \code{"glmmTMB"}, or \code{"inla"}.
#' @param ... Additional arguments passed to the underlying model-fitting function.
#'
#' @details
#' The function fits a localized regression model at each observation using all other
#' data points as training data. Weights are assigned based on pairwise distances among
#' predictor variables using an exponential decay function:
#' \deqn{w_i = \exp(-\theta d_i / \bar{d})}
#' where \eqn{d_i} is the Euclidean distance between the focal and training points, and \eqn{\bar{d}}
#' is the mean distance from the focal point to all others. The mean squared prediction
#' error is computed for each left-out observation, and the overall RMSE is returned.
#'
#' @return
#' A numeric value giving the root mean squared error (RMSE) from leave-one-out cross-validation.
#'
#' @examples
#' \dontrun{
#' ## Example using a simple linear model
#' df <- data.frame(
#'   x1 = rnorm(20),
#'   x2 = rnorm(20),
#'   y  = rnorm(20)
#' )
#' loocv(y ~ x1 + x2, data = df, theta = 0.5, model = "lm")
#' }
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#' @export

loocv <- function(formula,
                  data,
                  theta,
                  model = "lm",
                  ...) {

  ## get matrix X
  y <- all.vars(formula)[1]
  v_cnm_x <- attr(stats::terms(formula), "term.labels")
  v_cnm_x <- v_cnm_x[!grepl("\\|", v_cnm_x)]

  m_x <- data[, v_cnm_x] |>
    data.matrix()

  ## distance matrix for weighting
  m_dist <- stats::dist(m_x,
                        diag = TRUE,
                        upper = TRUE) |>
    data.matrix()

  v <- seq_len(nrow(data))

  ## define model fitting function once
  fit_fun <- switch(model,
                    lm = function(formula, data, w, ...) {
                      stats::lm(formula,
                                data = data,
                                weights = w,
                                ...)
                    },
                    glm = function(formula, data, w, ...) {
                      stats::glm(formula,
                                 data = data,
                                 weights = w,
                                 ...)
                    },
                    lmer = function(formula, data, w, ...) {
                      lme4::lmer(formula,
                                 data = data,
                                 weights = w,
                                 ...)
                    },
                    glmer = function(formula, data, w, ...) {
                      lme4::glmer(formula,
                                  data = data,
                                  weights = w,
                                  ...)
                    },
                    glmmTMB = function(formula, data, w, ...) {
                      glmmTMB::glmmTMB(formula,
                                       data = data,
                                       weights = w,
                                       ...)
                    },
                    inla = function(formula, data, w, ...) {
                      inla_lmer(formula,
                                data = data,
                                weights = w,
                                ...)
                    },
                    stop("Unsupported model type")
  )

  rmse <- sapply(v,
                 function(i) {

                   ## training data
                   df_train <- data[-i, ]
                   data[i, "w"] <- 1

                   ## calculate weights
                   ## note: "inla" is sensitive to scaling
                   mu_d <- mean(m_dist[i, -i])
                   w0 <- exp(-theta * m_dist[i, -i] / mu_d)
                   df_train$w <- (w0 / sum(w0)) * nrow(df_train)


                   lw <- fit_fun(formula,
                                 data = df_train,
                                 w = w,
                                 ...)

                   if (inherits(lw, "inla")) {
                     y0 <- point_predict(lw,
                                         newdata = data[i, , drop = FALSE])
                   } else {
                     y0 <- stats::predict(lw,
                                          newdata = data[i, , drop = FALSE])
                   }

                   y1 <- data[[y]][i]
                   eps <- (y1 - y0)^2

                   return(eps)
                 }) |>
    mean() |>
    sqrt()

  return(rmse)
}


#' Estimate the Equilibrium Density of a Localized Community Model
#'
#' Estimates the equilibrium density of the total community using a localized regression
#' approach with distance-based weighting. The equilibrium point is obtained iteratively,
#' updating the local regression until convergence is reached.
#'
#' @param formula A model formula specifying the response and predictor variables.
#' @param data A data frame containing the variables referenced in \code{formula}.
#' @param theta Numeric vector. Candidate values of the distance-weighting parameter, evaluated using leave-one-out cross-validation to identify the optimal value.
#' @param maxit Integer. Maximum number of iterations used to approximate the equilibrium density.
#' @param tol Numeric. Convergence threshold for the iterative update of the equilibrium estimate.
#' @param ... Additional arguments passed to the underlying model-fitting functions.
#'
#' @details
#' This function estimates the equilibrium density (\eqn{x^*}) of the total community
#' by iteratively fitting localized regression models. The equilibrium point is defined
#' as the value of the predictor (\eqn{x}) for which the predicted response equals zero.
#'
#' For each iteration, the model is refitted with weights assigned according to
#' the exponential distance-decay function:
#' \deqn{w_i = \exp(-\theta_0 d_i / \bar{d})}
#' where \eqn{d_i} is the absolute distance between the focal point and observation \eqn{i},
#' and \eqn{\bar{d}} is the mean distance. The optimal weighting parameter \eqn{\theta_0}
#' is selected using leave-one-out cross-validation via \code{\link{loocv}}.
#'
#' Iterations continue until the change in the equilibrium estimate between successive
#' steps falls below \code{tol}, or until \code{maxit} iterations are reached.
#'
#' @return
#' A numeric value representing the estimated equilibrium density (\eqn{x^*}) of the total community.
#' The returned object also includes the following attributes:
#' \itemize{
#'   \item \code{"gap"} – the final difference between successive estimates (convergence criterion)
#'   \item \code{"iteration"} – the total number of iterations performed
#' }
#'
#' @examples
#' \dontrun{
#' ## Example using a simple regression model
#' df <- data.frame(
#'   x = seq(1, 10, length.out = 20),
#'   y = 5 - 0.5 * seq(1, 10, length.out = 20) + rnorm(20, 0, 0.2)
#' )
#' xeq(y ~ x, data = df, theta = seq(0, 5, by = 0.5))
#' }
#'
#' @seealso \code{\link{loocv}} for model calibration using leave-one-out cross-validation.
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#' @export

xeq <- function(formula,
                data,
                theta = seq(0, 10, by = 0.5),
                maxit = 50,
                tol = 1E-4,
                ...) {

  ## initialize x_hat
  x_hat <- rep(NA, times = maxit)

  ## get the first estimate of x_star
  m0 <- stats::lm(formula,
                  data = data,
                  ...)

  x_hat[1] <- -stats::coef(m0)[1] / stats::coef(m0)[2]

  if (x_hat[1] < 0) {
    message("No feasible equilibrium exists.")
    return(NA)
  }

  ## x label
  x_cnm <- attr(stats::terms(formula), "term.labels")

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = "lm",
                           ...)
                   })

  theta0 <- theta[which.min(v_rmse)]

  for (i in seq_len(maxit - 1)) {
    d <- abs(data[, x_cnm, drop = TRUE] - x_hat[i])
    w0 <- exp(-theta0 * (d / mean(d)))
    data$w <- w0 / sum(w0)

    m <- stats::lm(formula,
                   data = data,
                   weights = w,
                   ...)

    x_hat[i + 1] <- -stats::coef(m)[1] / stats::coef(m)[2]
    gap <- abs(x_hat[i + 1] - x_hat[i])
    if (gap < tol) break
  }

  x_star <- x_hat[which.max(!is.na(x_hat))]
  attr(x_star, "gap") <- gap
  attr(x_star, "iteration") <- i + 1
  attr(x_star, "theta") <- theta0
  attr(x_star, "rmse") <- v_rmse

  return(x_star)
}


#' Estimate Historical Contingency (ψ)
#'
#' Estimates the degree of historical contingency (\eqn{\psi}) in community dynamics
#' based on the equilibrium total community density (\eqn{x^*}). The statistic
#' quantifies the probability that the intrinsic population growth rate at equilibrium
#' is negative, given parameter uncertainty in a localized regression model.
#'
#' @param formula A model formula specifying the response and predictor variables.
#' @param data A data frame containing the variables referenced in \code{formula}.
#' @param x_star Numeric value representing the equilibrium total community density,
#'   typically obtained from \code{\link{xeq}}.
#' @param theta Numeric vector. Candidate values of the distance-weighting parameter,
#'   evaluated via leave-one-out cross-validation to select the optimal value.
#' @param n_lag Character string specifying the column name for lagged population density.
#' @param nt_lag Character string specifying the column name for lagged total community density.
#' @param n_sim Integer. Number of random samples drawn from a multivariate normal
#'   distribution to propagate parameter uncertainty.
#' @param model Character string specifying the model type to fit. Must be one of
#'   \code{"lm"}, \code{"glm"}, \code{"lmer"}, \code{"glmer"}, or \code{"glmmTMB"}.
#' @param rescale Logical. If \code{TRUE}, predictor variables are standardized to mean 0 and SD 1.
#' @param ... Additional arguments passed to the underlying model-fitting functions.
#'
#' @details
#' This function computes historical contingency (\eqn{\psi}) as the proportion of
#' simulated intrinsic growth rates (\eqn{r}) that are negative at the equilibrium
#' total community density (\eqn{x^*}). The procedure is as follows:
#'
#' \enumerate{
#'   \item The optimal distance-weighting parameter (\eqn{\theta_0}) is estimated
#'         using leave-one-out cross-validation (\code{\link{loocv}}).
#'   \item A localized regression model is fit using exponential distance-based
#'         weights:
#'         \deqn{w_i = \exp(-\theta_0 d_i / \bar{d})}
#'         where \eqn{d_i} is the distance from the focal point and \eqn{\bar{d}} is
#'         the mean distance among observations.
#'   \item Parameter uncertainty is propagated by simulating coefficients from a
#'         multivariate normal distribution based on the model’s estimated covariance matrix.
#'   \item For each simulated parameter set, the log-scale intrinsic rate of increase
#'         (\eqn{r}) is computed as:
#'         \deqn{r = \beta_0 + \beta_2 x^*}
#'   \item The value of \eqn{\psi} is estimated as the proportion of simulations
#'         where \eqn{r < 0}.
#' }
#'
#' @return
#' A numeric value representing the estimated historical contingency (\eqn{\psi}),
#' with the following attributes:
#' \itemize{
#'   \item \code{"theta"} – the optimal distance-weighting parameter (\eqn{\theta_0})
#'   \item \code{"message"} – model convergence information (if available)
#' }
#'
#' @examples
#' \dontrun{
#' ## Example with synthetic data
#' df <- data.frame(
#'   log_r = runif(50, -1, 1),
#'   n_lag = rnorm(50, 10, 2),
#'   nt_lag = rnorm(50, 50, 10),
#'   species = rep(1:10, times = 5)
#' )
#'
#' x_star <- 10  # example equilibrium density
#' get_psi(
#'   formula = log_r ~ n_lag + nt_lag,
#'   data = df,
#'   x_star = x_star,
#'   theta = seq(0, 5, by = 0.5),
#'   model = "lm"
#' )
#' }
#'
#' @seealso
#' \code{\link{xeq}} for estimating equilibrium densities, and
#' \code{\link{loocv}} for selecting the optimal weighting parameter.
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#' @export

get_psi <- function(formula,
                    data,
                    x_star,
                    theta = seq(0, 10, by = 0.5),
                    n_lag = attr(stats::terms(formula), "term.labels")[1],
                    nt_lag = attr(stats::terms(formula), "term.labels")[2],
                    n_sim = 1000,
                    model = "lm",
                    rescale = TRUE,
                    ...) {

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = model,
                           ...)
                   })

  theta0 <- theta[which.min(v_rmse)]

  ## distance
  data[["d"]] <- sqrt(data[[n_lag]]^2 + (data[[nt_lag]] - x_star)^2)

  ## mean and sd
  v_mu <- sapply(data[, c(n_lag, nt_lag)], mean)
  v_sigma <- sapply(data[, c(n_lag, nt_lag)], stats::sd)

  ## scaled predictors
  if (rescale) {
    data[[n_lag]] <- scale(data[[n_lag]],
                           center = TRUE,
                           scale = TRUE)

    data[[nt_lag]] <- scale(data[[nt_lag]],
                            center = TRUE,
                            scale = TRUE)
  }

  ## get scaled weight
  w0 <- with(data, exp(-theta0 * (d / mean(d))))

  data$w <- w0 / sum(w0)

  m <- switch(model,
              lm = stats::lm(formula,
                             data = data,
                             weights = w,
                             ...),
              glm = stats::glm(formula,
                               data = data,
                               weights = w,
                               ...),
              lmer = lme4::lmer(formula,
                                data = data,
                                weights = w,
                                ...),
              glmer = lme4::glmer(formula,
                                  data = data,
                                  weights = w,
                                  ...),
              glmmTMB = glmmTMB::glmmTMB(formula,
                                         data = data,
                                         weights = w,
                                         ...),
              stop("Unsupported model type")
  )

  ## mean and var-covar matrix for simulated parameters
  v_id <- c("(Intercept)", n_lag, nt_lag)

  if (any(class(m) %in% c("lm", "glm"))) {

    v_mu <- stats::coef(m)[v_id]
    m_sigma <- stats::vcov(m)[v_id, v_id]

  } else if (any(class(m) %in% c("lmerMod", "glmerMod"))) {

    v_mu <- lme4::fixef(m)[v_id]
    m_sigma <- stats::vcov(m)[v_id, v_id]

  } else if (any(class(m) %in% c("glmmTMB"))) {

    v_mu <- glmmTMB::fixef(m)$cond[v_id]
    m_sigma <- stats::vcov(m)$cond[v_id, v_id]

  }

  ## get simulated parameters
  m_sim <- MASS::mvrnorm(n = n_sim,
                         mu = v_mu,
                         Sigma = m_sigma)

  if (rescale) {
    m_sim <- apply(m_sim,
                   MARGIN = 1,
                   FUN = function(z) {

                     ## intercept original
                     b0 <- z[1] -
                       (z[2] / v_sigma[n_lag]) * v_mu[n_lag] -
                       (z[3] / v_sigma[nt_lag]) * v_mu[nt_lag]

                     ## n0's effect original
                     b1 <- z[2] / v_sigma[n_lag]

                     ## nt0's effect original
                     b2 <- z[3] / v_sigma[nt_lag]

                     return(c(b0, b1, b2))
                   }) |>
      t()
  }

  ## get a vector of simulated log-scale r
  cnm <- colnames(m_sim)
  key <- paste("[Ii]ntercept", nt_lag, sep = "|")
  v_r <- apply(m_sim[, grepl(key, cnm)],
               MARGIN = 1,
               FUN = function(b) {
                 b[1] + b[2] * x_star
               })

  ## estimate psi
  psi <- mean(v_r < 0)

  if (class(m) %in% c("lm", "glm")) {
    note <- m$converged
  } else if (any(class(m) %in% c("lmerMod", "glmerMod"))) {
    note <- summary(m)$optinfo$conv$lme4$messages
  } else if (any(class(m) %in% c("glmmTMB"))) {
    note <- m$fit$message
  }

  attr(psi, "message") <- note
  attr(psi, "theta") <- theta0

  return(psi)
}
