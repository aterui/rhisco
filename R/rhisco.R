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
#' @param group Character string specifying the grouping variable in `data`.
#'   Each group is assumed to represent a species unit that is repeatedly
#'   observed across time. The function requires that every group contains exactly one
#'   observation per timestep (see `tc`). Defaults to \code{NULL}.
#' @param tc Character vector of possible column names representing the time index.
#'   The function searches for the first matching name in `data`. This column must define
#'   the temporal ordering of observations and must be consistent across all groups,
#'   such that every group has one observation in each timestep. Defaults to
#'   `c("t", "time", "timestep", "ts")`.
#' @param size Integer specifying the number of timesteps to include in the
#'   cross-validation subset. Instead of performing LOOCV across *all* timesteps,
#'   the function randomly selects `size` distinct time points and leaves out one
#'   timestep at a time within this subset. This reduces computation while maintaining
#'   representative temporal coverage. If `NULL` (default), all timesteps are used.
#' @param seed Integer scalar specifying a random seed.
#' @param type Type of distance weighting function. Either `"exp"` or `"gaussian"`.
#' @param method Character string specifying the scaling method. Options are `"mean"` or `"max"`.
#' @param REML Logical. If \code{TRUE}, Restricted Maximum Likelihood is used for estimation. Defaults to \code{FALSE}.
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
#'   \item{Exponential:}{\deqn{w_i = \exp(-\theta d_i)}}
#'   \item{Gaussian:}{\deqn{w_i = \exp(-\theta^2 (d_i)^2 / 2)}}
#' }
#' where \eqn{d_i} is the Euclidean distance between the focal (test) point and each training point,
#' scaled either by the mean or maximum distance from the focal point to all others (see argument \code{method}). Each weight is normalized
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
                  group = NULL,
                  tc = c("t", "time", "timestep", "ts"),
                  size = NULL,
                  seed = NULL,
                  type = "gaussian",
                  method = "max",
                  REML = FALSE,
                  ...) {

  ## define resample
  resample <- function(x, ...) x[sample.int(length(x), ...)]

  ## check if `group` column exists in the data, then sort
  if (is.null(group)) {

    group <- "species"
    data$species <- 1

  } else {

    if (!(group %in% colnames(data)))
      stop(paste(group, "is not found in the dataframe."))

    data <- data[order(data[[group]]), ] |>
      transform(index = seq_len(nrow(data)))

  }

  ## check if "t" column exists
  tcol <- intersect(tc, colnames(data))[1]
  if (is.na(tcol)) stop("No time column found.")

  n_unq_t <- unique(table(data[[tcol]]))
  n_gr <- unique(data[[group]]) |>
    length()

  if (n_unq_t != n_gr)
    stop(paste("All operational units in 'group' must contain one observation in each timestep"))

  ## define dfun
  dfun <- get_dfun(type = type,
                   method = method)

  ## get matrix X and distance
  y <- all.vars(formula)[1]
  v_cnm_x0 <- attr(stats::terms(formula), "term.labels")
  v_cnm_x <- v_cnm_x0[!grepl("\\|", v_cnm_x0)]

  if (model %in% c("lm", "glm") && any(grepl("\\|", v_cnm_x0)))
    stop(paste(sQuote(model), "does not allow random effects"))

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
                                 REML = REML,
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
                                       REML = REML,
                                       ...)
                    },
                    stop("Unsupported model type")
  )

  ## select random subset
  if (is.null(size))
    size <- length(unique(data[[tcol]]))

  v_t <- gen_values(seed, function() {
    unique(data[[tcol]]) |>
      resample(size = size,
               replace = FALSE) |>
      sort()
  })

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
                     stats::setNames(c("d", "t", "g"))

                   ## - calculate raw weights by group
                   w0 <- with(
                     df_dist_i,
                     ave(d,
                         g,
                         FUN = function(x) dfun(x = x, theta = theta))
                   )

                   ## - normalize raw weights by group
                   df_train$w0 <- w0
                   df_train$w <- with(
                     df_train,
                     ave(w0,
                         df_train[[group]],
                         FUN = function(x) x / sum(x))
                   )

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
#' @inheritParams loocv
#' @param theta Numeric vector. Candidate values of the distance-weighting parameter,
#'   evaluated using leave-one-out cross-validation to identify the optimal value.
#' @param maxit Integer. Maximum number of iterations used to approximate the equilibrium density.
#' @param tol Numeric. Convergence threshold for the iterative update of the equilibrium estimate.
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
#' The optimal weighting parameter \eqn{\theta_0} is selected via leave-one-out cross-validation
#' using \code{\link{loocv}}.
#'
#' Iterations continue until the change in the equilibrium estimate between successive
#' steps falls below \code{tol}, or until \code{maxit} iterations are reached.
#'
#' @return
#' A numeric value representing the estimated equilibrium density (\eqn{x^*}) of the total community.
#' The returned object also includes:
#' \itemize{
#'   \item \code{"gap"} – final difference between successive estimates (convergence criterion)
#'   \item \code{"iteration"} – number of iterations performed
#'   \item \code{"theta"} – the selected optimal weighting parameter
#'   \item \code{"rmse"} – the LOOCV error profile used for parameter selection
#' }
#'
#' @seealso \code{\link{loocv}} for model calibration.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   x = seq(1, 10, length.out = 20),
#'   y = 5 - 0.5 * seq(1, 10, length.out = 20) + rnorm(20, 0, 0.2)
#' )
#' xeq(y ~ x, data = df, theta = seq(0, 5, by = 0.5))
#' }
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#' @export

xeq <- function(formula,
                data,
                tc = c("t", "time", "timestep", "ts"),
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

  if (length(x_cnm) > 1)
    stop(paste(sQuote(formula)[3], "can contain only one predictor"))

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = "lm",
                           tc = tc,
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
#' based on the equilibrium total community density (\eqn{x^*}). The statistic quantifies
#' the probability that the intrinsic population growth rate at equilibrium is negative,
#' given parameter uncertainty in a localized regression model.
#'
#' @inheritParams loocv
#' @param x_star Numeric value representing the equilibrium total community density,
#'   typically obtained from \code{\link{xeq}}.
#' @param theta Numeric vector. Candidate values of the distance-weighting parameter,
#'   evaluated via leave-one-out cross-validation to identify the optimal value.
#' @param n_lag Character string specifying the column name for the lagged population density.
#' @param nt_lag Character string specifying the column name for the lagged total community density.
#' @param rescale Logical. If \code{TRUE}, predictor variables are standardized to mean 0 and SD 1.
#' @param dist.method Character string specifying the scaling method. Options are `"mean"` or `"max"`.
#' @param psi.method Character string specifying how species-level variation is
#'   incorporated into \eqn{\psi}. Options are \code{"marginal"} or \code{"mean"}.
#'   \code{"marginal"} integrates over species-level variation in the intercept
#'   (i.e., \code{var_r} from the random effect) when computing the probability
#'   that the invasion growth rate is negative, so that \eqn{\psi} reflects both
#'   estimation uncertainty and genuine between-species variation.
#'   \code{"mean"} ignores species-level variation and computes \eqn{\psi} based
#'   solely on the mean invasion growth rate and its estimation
#'   uncertainty, equivalent to asking whether the average species
#'   cannot invade. Defaults to \code{"marginal"}.
#'
#' @details
#' This function quantifies historical contingency (\eqn{\psi}) by evaluating uncertainty
#' in the intrinsic growth rate at the estimated equilibrium density (\eqn{x^*}). The
#' optimal distance-weighting parameter (\eqn{\theta_0}) is selected via leave-one-out
#' cross-validation using \code{\link{loocv}}.
#'
#' @return
#' A numeric value representing historical contingency (\eqn{\psi}). The returned object
#' also includes:
#' \itemize{
#'   \item \code{"theta"} – selected optimal distance-weighting parameter
#'   \item \code{"message"} – model or convergence diagnostics (if present)
#' }
#'
#' @examples
#' \dontrun{
#' ## this is fake data with no ecological meaning
#' ## users must prepare the data appropriately
#' df <- data.frame(
#'   log_r = runif(50, -1, 1),
#'   n_lag = rnorm(50, 10, 2),
#'   nt_lag = rnorm(50, 50, 10),
#'   species = rep(1:10, each = 5)
#' )
#'
#' x_star <- 10
#' get_psi(
#'   formula = log_r ~ n_lag + nt_lag + (1 | species),
#'   data = df,
#'   x_star = x_star,
#'   theta = seq(0, 5, by = 0.5),
#'   model = "lmer"
#' )
#' }
#'
#' @seealso \code{\link{xeq}} for equilibrium estimation; \code{\link{loocv}} for
#' cross-validation–based parameter selection.
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#' @export

get_psi <- function(formula,
                    data,
                    x_star,
                    group = extract_group(formula),
                    tc = c("t", "time", "timestep", "ts"),
                    theta = seq(0, 10, by = 0.5),
                    n_lag = attr(stats::terms(formula), "term.labels")[1],
                    nt_lag = attr(stats::terms(formula), "term.labels")[2],
                    model = "lm",
                    rescale = TRUE,
                    size = NULL,
                    seed = NULL,
                    type = "gaussian",
                    dist.method = "max",
                    psi.method = c("marginal", "mean"),
                    REML = FALSE,
                    ...) {

  # check input -------------------------------------------------------------

  ## psi.method matching
  psi.method <- match.arg(psi.method)

  ## random effect arguments
  if (is.null(group))
    stop("Missing 'group' argument.")

  ## check if the random effect is single and intercept only
  re_check(formula)

  ## do random effects exist?
  re_idx <- !is.null(extract_group(formula))

  ## number of species
  nsp <- length(unique(data[[group]]))

  ## number of terms, term-name mismatch?
  v_terms <- reformulas::nobars(formula) |>
    stats::terms() |>
    attr("term.labels")

  if (length(v_terms) != 2)
    stop(paste(length(v_terms), "fixed effects detected. The number of predictors must be two."))

  if (length(intersect(v_terms, c(n_lag, nt_lag))) != 2)
    stop(paste(n_lag, "and/or", nt_lag, "are not found in the formula."))

  # reformat data -----------------------------------------------------------

  ## mean and sd
  v_mu <- sapply(data[, c(n_lag, nt_lag)], mean)
  v_sigma <- sapply(data[, c(n_lag, nt_lag)], stats::sd)
  x0 <- 0

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

    x0 <- -v_mu[n_lag] / v_sigma[n_lag]
    x_star <- (x_star - v_mu[nt_lag]) / v_sigma[nt_lag]
  }

  ## variable names
  v_id <- c("(Intercept)", n_lag, nt_lag)

  # define theta ------------------------------------------------------------

  seen <- new.env(parent = emptyenv())
  seen$warnings <- character(0)
  seen$messages <- character(0)

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {

                     withCallingHandlers(
                       loocv(formula,
                             data = data,
                             theta = x,
                             model = model,
                             group = group,
                             tc = tc,
                             size = size,
                             seed = seed,
                             type = type,
                             method = dist.method,
                             REML = REML,
                             ...),

                       warning = function(w) {
                         msg <- conditionMessage(w)
                         if (!(msg %in% seen$warnings))
                           seen$warnings <- c(seen$warnings, msg)
                         invokeRestart("muffleWarning")
                       },

                       message = function(m) {
                         msg <- conditionMessage(m)
                         if (!(msg %in% seen$messages))
                           seen$messages <- c(seen$messages, msg)
                         invokeRestart("muffleMessage")
                       }

                     ) # withCallingHandlers

                   }) # sapply

  names(v_rmse) <- theta
  theta0 <- theta[which.min(v_rmse)]

  ## define dfun
  dfun <- get_dfun(type = type,
                   method = dist.method)

  # get X* - <K_0> ----------------------------------------------------------

  ## euclidean distance to the point of approximation; then, get scaled weight
  v_d0 <- abs(data[[nt_lag]] - x_star)
  data[["w0"]] <- dfun(v_d0, theta = theta0)

  ## model fit
  m0 <- switch(model,
               lm = stats::lm(formula,
                              data = data,
                              weights = w0,
                              ...),
               lmer = lme4::lmer(formula,
                                 data = data,
                                 weights = w0,
                                 REML = REML,
                                 ...),
               glmmTMB = glmmTMB::glmmTMB(formula,
                                          data = data,
                                          weights = w0,
                                          REML = REML,
                                          ...),
               stop("Unsupported model type")
  )

  ## get mean beta, vcov, random var
  l_m0_par <- get_params(m0, idx = v_id, re = re_idx)

  ## estimate X*_{-i}
  if (re_idx) {
    m_b0 <- l_m0_par[["rbeta"]]
  } else {
    m_b0 <- l_m0_par[["beta"]]
  }

  if (rescale) {
    # obtain unscaled parameters
    # - v_intra; m_b[2, ] represents intra - inter
    # - thus, m_b[2, ] + m_b[3, ] (inter) transforms it back to pure intra
    np <- length(v_id)
    m_std <- matrix(0, np, np)
    diag(m_std) <- c(1, 1 / v_sigma)
    m_std[1, seq(2, np, by = 1)] <- -v_mu / v_sigma
    m_b <- m_std %*% m_b0
    v_intra <- m_b[2, ] + m_b[3, ]

    # get k0 with unscaled x_star (= x_star_usc)
    x_star_usc <- x_star * v_sigma[nt_lag] + v_mu[nt_lag]
    v_k0 <- -(m_b[1, ] + m_b[3, ] * x_star_usc) / v_intra
    v_k0 <- ifelse(v_k0 > 0, v_k0, 0)
    x_star_i_usc <- mean(x_star_usc - v_k0)
    x_star_i <- (x_star_i_usc - v_mu[nt_lag]) / v_sigma[nt_lag]
  } else {
    m_b <- m_b0
    v_intra <- m_b[2, ] + m_b[3, ]

    v_k0 <- -(m_b[1, ] + m_b[3, ] * x_star) / v_intra
    v_k0 <- ifelse(v_k0 > 0, v_k0, 0)
    x_star_i <- mean(x_star - v_k0)
  }

  # psi ---------------------------------------------------------------------

  ## a vector of constants; assigned as matrix
  m_const <- matrix(c(1, x0, x_star_i),
                    nrow = 3,
                    ncol = 1)
  rownames(m_const) <- v_id

  ## euclidean distance to the point of approximation; then, get scaled weight
  v_d <- sqrt((data[[n_lag]] - x0)^2 + (data[[nt_lag]] - x_star_i)^2)
  data[["w"]] <- dfun(v_d, theta = theta0)

  m <- switch(model,
              lm = stats::lm(formula,
                             data = data,
                             weights = w,
                             ...),
              lmer = lme4::lmer(formula,
                                data = data,
                                weights = w,
                                REML = REML,
                                ...),
              glmmTMB = glmmTMB::glmmTMB(formula,
                                         data = data,
                                         weights = w,
                                         REML = REML,
                                         ...),
              stop("Unsupported model type")
  )

  ## get mean beta, vcov, random var
  l_m_par <- get_params(m, idx = v_id, re = re_idx)

  ## mean r
  mu_r <- drop(t(m_const) %*% l_m_par[["beta"]])

  ## sd r, fixed effect variance + random effect variance [species]
  ## - fixed terms
  m_sigma <- l_m_par[["sigma"]]
  var_x <- drop(t(m_const) %*% m_sigma %*% m_const)

  ## - random terms
  if (!re_idx) {
    var_r0 <- 0
    psi.method <- "mean"
    message("No random effect in the model; psi.method has been set to 'mean' automatically.")
  } else {
    # select random-effect terms present in the model
    m_rvar <- l_m_par[["rvar"]]
    re_terms <- intersect(v_id, rownames(m_rvar))
    m_const_re <- m_const[match(re_terms, v_id), , drop = FALSE]

    # compute variance contribution from random effects
    var_r0 <- drop(t(m_const_re) %*% m_rvar[re_terms, re_terms] %*% m_const_re)
  }

  var_r <- if (psi.method == "mean") 0 else var_r0

  ## - total SD
  sd_tot <- sqrt(var_x + var_r)

  ## psi, integrating species effects out
  psi <- stats::pnorm(0, mean = mu_r, sd = sd_tot)

  ## export
  if (any(class(m) %in% c("lm", "glm"))) {
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
  attr(psi, "sd") <- c(sd_mu = sqrt(var_x),
                       sd_r = if (re_idx) sqrt(var_r0) else NA)

  return(psi)
}

