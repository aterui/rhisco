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
#'   \code{"lm"}, \code{"glm"}, \code{"lmer"}, \code{"glmer"}, or \code{"glmmTMB"}.
#' @param size Integer scalar. The number of subset data points for leave-one-out cross validation.
#' @param seed Integer scalar specifying a random seed.
#' @param type Type of distance weighting function. Either `"exp"` or `"gaussian"`.
#' @param method Character string specifying the scaling method. Options are `"mean"` or `"max"`.
#' @param ... Additional arguments passed to the underlying model-fitting function.
#'
#' @details
#' This function evaluates localized regression models using leave-one-out cross-validation (LOOCV)
#' with distance-based weighting. For each observation selected in the subset, the model is trained
#' on all remaining data points, and the held-out observation is used for prediction. The squared
#' prediction errors across all test points are averaged to compute the overall root mean squared error (RMSE).
#'
#' Distance-based weights are computed using one of two kernel functions:
#' \describe{
#'   \item{Exponential:}{\deqn{w_i = \exp(-\theta d_i / \bar{d})}}
#'   \item{Gaussian:}{\deqn{w_i = \exp(- (d_i / \bar{d})^2 / (2\theta^2))}}
#' }
#' where \eqn{d_i} is the Euclidean distance between the focal (test) point and each training point,
#' and \eqn{\bar{d}} is the mean or maximum distance from the focal point to all others. Each weight is normalized
#' such that the total weight equals the number of training samples, maintaining comparability across fits.
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
                  group = "species",
                  tc = c("t", "time", "timestep", "ts"),
                  size = NULL,
                  seed = NULL,
                  type = "gaussian",
                  method = "max",
                  ...) {

  ## define resample
  resample <- function(x, ...) x[sample.int(length(x), ...)]

  ## check if `group` column exists in the data, then sort
  if (!(group %in% colnames(data)))
    stop(paste(group, "is not found in the dataframe."))

  data <- data[order(data[[group]]), ] |>
    transform(index = seq_len(nrow(data)))

  ## check if "t" column exists
  tcol <- intersect(tc, colnames(data))[1]
  if (is.na(tcol)) stop("No time column found.")

  n_unq_t <- unique(table(data[[tcol]]))
  n_gr <- unique(data[[group]]) |>
    length()

  if (n_unq_t != n_gr)
    stop(paste("All constituent elements in",
               sQuote(group),
               "must contain one observation in each timestep"))

  ## define dfun
  dfun <- get_dfun(type = type,
                   method = method)

  ## get matrix X and distance
  y <- all.vars(formula)[1]
  v_cnm_x <- attr(stats::terms(formula), "term.labels")
  v_cnm_x <- v_cnm_x[!grepl("\\|", v_cnm_x)]

  list_dist <- split(data, data[[group]]) |>
    lapply(FUN = function(x) {

      df_x <- x[, v_cnm_x] |>
        stats::dist(diag = TRUE,
                    upper = TRUE) |>
        data.matrix() |>
        as.data.frame()

      colnames(df_x) <- x[[tcol]]
      df_x$t <- x[[tcol]]
      df_x$g <- x[[group]]
      return(df_x)
    })

  df_dist <- do.call(rbind,
                     list_dist)

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
                    stop("Unsupported model type")
  )

  ## select random subset
  if (!is.null(seed))
    set.seed(seed)

  if (is.null(size))
    size <- length(unique(data[[tcol]]))

  v_t <- unique(data[[tcol]]) |>
    resample(size = size,
             replace = FALSE) |>
    sort()

  data$w <- 1

  rmse <- sapply(seq_len(length(v_t)),
                 function(i) {

                   ## get vector indices for "leave-out"
                   v_idx <- which(data[[tcol]] == v_t[i])

                   ## training data
                   df_train <- data[-v_idx, ]

                   ## calculate weights
                   ## - remove "leave-out" data points
                   cnm <- c(as.character(v_t[i]), "t", "g")
                   df_dist_i <- df_dist[-v_idx, cnm] |>
                     setNames(c("d", "t", "g"))

                   ## - calculate raw weights by group
                   w0 <- with(df_dist_i,
                              ave(d, g, FUN = function(x) dfun(x = x,
                                                               theta = theta)))

                   ## - normalize raw weights by group
                   df_train$w0 <- w0
                   df_train$w <- with(df_train,
                                      ave(w0, df_train[[group]],
                                          FUN = function(x) x / sum(x)))

                   ## fit model
                   lw <- fit_fun(formula,
                                 data = df_train,
                                 w = w,
                                 ...)

                   y0 <- stats::predict(lw,
                                        newdata = data[v_idx, , drop = FALSE],
                                        type = "response")

                   y1 <- data[[y]][v_idx]
                   eps <- sum((y1 - y0)^2)
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
#' @param type Type of distance weighting function. Either `"exp"` or `"gaussian"`.
#' @param method Character string specifying the scaling method. Options are `"mean"` or `"max"`.
#' @param ... Additional arguments passed to the underlying model-fitting functions.
#'
#' @details
#' This function estimates the equilibrium density (\eqn{x^*}) of the total community
#' by iteratively fitting localized regression models that apply distance-based weighting
#' to observations. The equilibrium point is defined as the predictor value (\eqn{x})
#' at which the predicted response equals zero, representing a steady-state condition
#' in the modeled system.
#'
#' For each iteration, the model is refitted with weights determined by a distance-decay
#' function, reflecting the influence of nearby data points on local model estimation.
#' Two types of weighting kernels are supported:
#' \describe{
#'   \item{Exponential:}{\deqn{w_i = \exp(-\theta_0 d_i / \bar{d})}}
#'   \item{Gaussian:}{\deqn{w_i = \exp(-(d_i / \bar{d})^2 / (2\theta_0^2))}}
#' }
#' where \eqn{d_i} is the absolute distance between the focal point and observation \eqn{i},
#' and \eqn{\bar{d}} is the mean pairwise distance. Weights are normalized so that their
#' total equals the number of training samples, ensuring comparability across fits.
#'
#' The optimal weighting parameter \eqn{\theta_0} is determined via
#' leave-one-out cross-validation (LOOCV) using \code{\link{loocv}}. This procedure identifies
#' the value of \eqn{\theta_0} that minimizes the root mean squared prediction error (RMSE),
#' balancing model smoothness and locality.
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
                theta = seq(0.5, 10, by = 0.5),
                maxit = 50,
                tol = 1E-4,
                type = "gaussian",
                method = "max",
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
                           type = type,
                           method = method,
                           ...)
                   })

  theta0 <- theta[which.min(v_rmse)]

  ## define dfun
  dfun <- get_dfun(type = type,
                   method = method)

  for (i in seq_len(maxit - 1)) {
    d <- abs(data[, x_cnm, drop = TRUE] - x_hat[i])
    w0 <- dfun(d, theta = theta0)
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
#' @param size Integer. Sub-sample size for leave-one-out cross validation.
#' @param seed Integer. Random seed for leave-one-out cross validation.
#' @param type Type of distance weighting function. Either `"exp"` or `"gaussian"`.
#' @param method Character string specifying the scaling method. Options are `"mean"` or `"max"`.
#' @param ... Additional arguments passed to the underlying model-fitting functions.
#'
#' @details
#' This function quantifies historical contingency (\eqn{\psi})—the degree to which
#' current community states depend on past conditions—based on the equilibrium total
#' community density (\eqn{x^*}). Specifically, \eqn{\psi} represents the probability
#' that the intrinsic population growth rate at equilibrium is negative, accounting
#' for parameter uncertainty in a localized regression model.
#'
#' The procedure involves the following steps:
#' \enumerate{
#'   \item The optimal distance-weighting parameter (\eqn{\theta_0}) is identified
#'         via leave-one-out cross-validation (\code{\link{loocv}}), minimizing
#'         prediction error across candidate \eqn{\theta} values.
#'   \item A localized regression model is then fitted using distance-based weights
#'         derived from one of two kernel functions:
#'         \describe{
#'           \item{Exponential:}{\deqn{w_i = \exp(-\theta_0 d_i / \bar{d})}}
#'           \item{Gaussian:}{\deqn{w_i = \exp(-(d_i / \bar{d})^2 / (2\theta_0^2))}}
#'         }
#'         where \eqn{d_i} is the Euclidean distance between the focal (test) point and
#'         observation \eqn{i}, and \eqn{\bar{d}} is the mean pairwise distance.
#'         Weights are normalized to sum to one to maintain scale consistency.
#'   \item Parameter uncertainty is propagated by simulating regression coefficients
#'         from a multivariate normal distribution using the model’s estimated mean
#'         vector and covariance matrix.
#'   \item For each simulated parameter set, the intrinsic growth rate (\eqn{r})
#'         at the equilibrium density (\eqn{x^*}) is computed as:
#'         \deqn{r = \beta_0 + \beta_2 x^*}
#'   \item The historical contingency metric (\eqn{\psi}) is then estimated as the
#'         proportion of simulations where \eqn{r < 0}, reflecting the probability that
#'         equilibrium conditions yield a declining population under parameter uncertainty.
#' }
#'
#' The LOOCV step can optionally use a random subset of the data (\code{size}) for
#' computational efficiency. The choice of \code{type} ("exp" or "gaussian") determines
#' how spatial proximity affects weighting. Supported model types include
#' \code{"lm"}, \code{"glm"} (from \pkg{stats}),
#' \code{"lmer"}, \code{"glmer"} (from \pkg{lme4}), and
#' \code{"glmmTMB"} (from \pkg{glmmTMB}).
#'
#' This method provides a flexible framework to evaluate the strength of historical
#' contingency in community dynamics, integrating spatially localized regression,
#' cross-validation–based parameter optimization, and uncertainty propagation.
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
                    theta = seq(0.5, 10, by = 0.5),
                    n_lag = attr(stats::terms(formula), "term.labels")[1],
                    nt_lag = attr(stats::terms(formula), "term.labels")[2],
                    n_sim = 1000,
                    model = "lm",
                    rescale = TRUE,
                    size = min(100, nrow(data)),
                    seed = NULL,
                    type = "gaussian",
                    method = "max",
                    ...) {

  # reformat data -----------------------------------------------------------

  ## Euclidean distance to the point of approximation
  data[["d"]] <- sqrt(data[[n_lag]]^2 + (data[[nt_lag]] - x_star)^2)

  ## mean and sd
  v_mu <- sapply(data[, c(n_lag, nt_lag)], mean)
  v_sigma <- sapply(data[, c(n_lag, nt_lag)], stats::sd)

  ## scaled predictors
  if (rescale) {
    data[[n_lag]] <- scale(data[[n_lag]],
                           center = TRUE,
                           scale = TRUE) |>
      c()

    data[[nt_lag]] <- scale(data[[nt_lag]],
                            center = TRUE,
                            scale = TRUE) |>
      c()
  }

  # fit ---------------------------------------------------------------------

  if (is.null(seed)) seed <- stats::rpois(1, 100)

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = model,
                           size = size,
                           seed = seed,
                           type = type,
                           method = method,
                           ...)
                   })

  names(v_rmse) <- theta
  theta0 <- theta[which.min(v_rmse)]

  ## define dfun
  dfun <- get_dfun(type = type,
                   method = method)

  ## get scaled weight
  w0 <- with(data, dfun(d, theta = theta0))
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

    v_beta <- stats::coef(m)[v_id]
    m_sigma <- stats::vcov(m)[v_id, v_id]

  } else if (any(class(m) %in% c("lmerMod", "glmerMod"))) {

    v_beta <- lme4::fixef(m)[v_id]
    m_sigma <- stats::vcov(m)[v_id, v_id]

  } else if (any(class(m) %in% c("glmmTMB"))) {

    v_beta <- glmmTMB::fixef(m)$cond[v_id]
    m_sigma <- stats::vcov(m)$cond[v_id, v_id]

  }

  ## get simulated parameters
  m_sim <- MASS::mvrnorm(n = n_sim,
                         mu = v_beta,
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

  # Convert everything to character
  if (is.list(note)) {
    note <- paste(unlist(note), collapse = " | ")
  } else if (length(note) != 0) {
    note <- as.character(note)  # assign back to note
  } else {
    note <- NA_character_
  }

  attr(psi, "message") <- note
  attr(psi, "theta") <- theta0
  attr(psi, "rmse") <- v_rmse

  return(psi)
}
