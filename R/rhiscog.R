#' XXX
#'
#' @param formula
#' @param data
#' @param x_star
#' @param theta
#' @param n_sim
#' @param model
#' @param ...
#'
#' @importFrom stats dpois qpois
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
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
                   function(x) loocv(formula,
                                     data = data,
                                     theta = x,
                                     model = model)
  )

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
              b0 <- z[1] - (z[2] / v_sigma["n0"]) * v_mu["n0"] - (z[3] / v_sigma["nt0"]) * v_mu["n0"]

              ## n0's effect original
              b1 <- z[2] / v_sigma["n0"]

              ## nt0's effect original
              b2 <- z[3] / v_sigma["nt0"]

              return(c(b0, b1, b2))
            }) |>
      t()
  }

  ## get a vector of simulated log-scale r
  cnm <- colnames(msim)
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
