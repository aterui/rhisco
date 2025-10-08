#' Leave-One-Out Cross-Validation for Localized Regression
#'
#' Computes the root mean squared error (RMSE) using leave-one-out cross-validation
#' for localized regression models with distance-based weighting.
#'
#' @param formula A model formula specifying the response and predictor variables.
#' @param data A data frame containing the variables used in \code{formula}.
#' @param theta Numeric vector. Candidate values for the weighting parameter, evaluated using leave-one-out cross-validation.
#' @param model Character. Type of model to fit: one of \code{"lm"}, \code{"glm"}, \code{"lmer"}, or \code{"glmer"}.
#' @param ... Additional arguments passed to the underlying model-fitting functions.
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
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

  ## define model fitting function once
  fit_fun <- switch(model,
                    lm = function(formula, data, w, ...) {
                      lm(formula,
                         data = data,
                         weights = w,
                         ...)
                    },
                    glm = function(formula, data, w, ...) {
                      glm(formula,
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
                    stop("Unsupported model type")
  )

  rmse <- sapply(v,
                 function(i) {

                   ## training data
                   df_train <- data[-i, ]

                   ## calculate weights
                   mu_d <- mean(m_dist[i, -i])
                   df_train$w <- exp(-theta * m_dist[i, -i] / mu_d)

                   lw <- fit_fun(formula,
                                 data = df_train,
                                 w = w,
                                 ...)

                   y0 <- predict(lw, newdata = data[i, , drop = FALSE])
                   y1 <- data[i, y]
                   eps <- (y1 - y0)^2

                   return(eps)
                 }) |>
    mean() |>
    sqrt()

  return(rmse)
}


#' Estimate the Equilibrium Density of the Total Community
#'
#' Computes the equilibrium density of the total community based on a localized regression model.
#'
#' @param formula A model formula specifying the response and predictors.
#' @param data A data frame containing the variables in \code{formula}.
#' @param maxit Integer. Maximum number of iterations used to approximate the localized regression near the equilibrium.
#' @param tol Numeric. Tolerance threshold for determining convergence of the equilibrium estimate.
#' @param ... Additional arguments passed to other methods.
#' @inheritParams loocv
#'
#' @return Numeric value representing the estimated equilibrium density of the total community.
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#'
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
  m0 <- lm(formula,
           data = data,
           ...)

  x_hat[1] <- -coef(m0)[1] / coef(m0)[2]

  if (x_hat[1] < 0) stop("No feasible equilibrium exists.")

  ## x label
  x_cnm <- attr(terms(formula), "term.labels")

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = "lm")
                   })

  theta0 <- theta[which.min(v_rmse)]

  for (i in seq_len(maxit - 1)) {
    d <- abs(data[, x_cnm, drop = TRUE] - x_hat[i])
    w0 <- exp(-theta0 * (d / mean(d)))
    data$w <- w0 / sum(w0)

    m <- lm(formula,
            data = data,
            weights = w,
            ...)

    x_hat[i + 1] <- -coef(m)[1] / coef(m)[2]
    gap <- abs(x_hat[i + 1] - x_hat[i])
    if (gap < tol) break
  }

  x_star <- x_hat[which.max(!is.na(x_hat))]
  attr(x_star, "gap") <- gap
  attr(x_star, "iteration") <- i + 1

  return(x_star)
}


#' Estimate Historical Contingency (psi)
#'
#' Computes historical contingency (psi) based on equilibrium community density estimates.
#'
#' @param x_star Numeric. Equilibrium total community density. The estimate should be obtained from \code{xeq}.
#' @param n_sim Integer. Number of random samples drawn from a multivariate normal distribution.
#' @param ... Additional arguments passed to other methods.
#' @inheritParams xeq
#'
#' @author Akira Terui (\email{hanabi0111@gmail.com})
#'
#' @return A numeric value representing the estimated historical contingency (psi).
#'
#' @export

get_psi <- function(formula,
                    data,
                    x_star,
                    theta = seq(0, 10, by = 0.5),
                    n_sim = 1000,
                    model = "lm",
                    ...) {

  ## estimate best theta with leave-one-out cross validation
  v_rmse <- sapply(theta,
                   function(x) {
                     loocv(formula,
                           data = data,
                           theta = x,
                           model = model)
                   })

  theta0 <- theta[which.min(v_rmse)]

  ## scaled predictors
  data <- transform(data,
                    scl_n0 = scale(n0,
                                   center = TRUE,
                                   scale = TRUE),
                    scl_nt0 = scale(nt0,
                                    center = TRUE,
                                    scale = TRUE),
                    d = sqrt(n0^2 + (nt0 - x_star)^2)
  )

  v_mu <- sapply(data[, c("n0", "nt0")], mean)
  v_sigma <- sapply(data[, c("n0", "nt0")], sd)

  ## get scaled weight
  w0 <- with(data,
             exp(-theta0 * (d / mean(d))))

  data$w <- w0 / sum(w0)

  m <- switch(model,
              lm = lm(formula,
                      data = data,
                      weights = w,
                      ...),
              glm = glm(formula,
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
              stop("Unsupported model type")
  )

  if (any(class(m) %in% c("lm", "glm"))) {
    m_sim <- MASS::mvrnorm(n = n_sim,
                           mu = coef(m),
                           Sigma = vcov(m)) |>
      apply(MARGIN = 1,
            function(z) {

              ## intercept original
              b0 <- z[1] - (z[2] / v_sigma["n0"]) * v_mu["n0"] - (z[3] / v_sigma["nt0"]) * v_mu["n0"]

              ## n0's effect original
              b1 <- z[2] / v_sigma["n0"]

              ## nt0's effect original
              b2 <- z[3] / v_sigma["nt0"]

              return(c(b0, b1, b2))
            }) |>
      t()
  } else if (any(class(m) %in% c("lmerMod", "glmerMod"))) {
    m_sim <- MASS::mvrnorm(n = n_sim,
                           mu = m@beta,
                           Sigma = vcov(m)) |>
      apply(MARGIN = 1,
            function(z) {

              ## intercept original
              b0 <- z[1] -
                (z[2] / v_sigma["n0"]) * v_mu["n0"] -
                (z[3] / v_sigma["nt0"]) * v_mu["nt0"]

              ## n0's effect original
              b1 <- z[2] / v_sigma["n0"]

              ## nt0's effect original
              b2 <- z[3] / v_sigma["nt0"]

              return(c(b0, b1, b2))
            }) |>
      t()
  }

  ## get a vector of simulated log-scale r
  cnm <- colnames(m_sim)
  v_r <- apply(m_sim[, grepl("[Ii]ntercept|nt0", cnm)],
               MARGIN = 1,
               FUN = function(b) {
                 b[1] + b[2] * x_star
               })

  ## estimate psi
  psi <- mean(v_r < 0)

  if (class(m) %in% c("lm", "glm")) {
    note <- m$converged
  } else {
    note <- summary(m)$optinfo$conv$lme4$messages
  }
  attr(psi, "message") <- note
  attr(psi, "theta") <- theta0

  return(psi)
}
