#' XXX
#'
#' @param data
#' @param n_rep
#' @param n_sim
#' @param theta
#' @param maxit
#' @param method
#'
#' @importFrom stats dpois qpois
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

get_psi <- function(data,
                    n_rep = 100,
                    n_sim = 1000,
                    theta = 1,
                    maxit = 100,
                    method = "lm") {

  ## get best theta through S-map variant
  v_rmse1 <- with(data, {
    sapply(theta,
           function(a) loocv(data_c,
                             y = "log_r",
                             x = "nt0",
                             model = "rlm",
                             maxit = maxit,
                             theta = a))
  })

  theta1 <- theta[which.min(v_rmse1)]

  ## get estimate for community regulation
  ## repeated weighted regression to approximate a partial derivative
  list_c <- with(data, {
    m0 <- MASS::rlm(log_r ~ nt0,
                    data = data_c,
                    maxit = maxit)
    x_hat <- -coef(m0)[1] / coef(m0)[2]

    for (i in 1:n_rep) {
      d <- sqrt((data_c$log_r)^2 + (data_c$nt0 - x_hat[i])^2)
      w0 <- exp(-theta1 * (d / mean(d)))
      w <- w0 / sum(w0)

      m0 <- MASS::rlm(log_r ~ nt0,
                      data = data_c,
                      weights = w,
                      maxit = maxit)
      x_hat[i + 1] <- -coef(m0)[1] / coef(m0)[2]
    }

    mu_x_hat <- mean(x_hat[ceiling(n_rep * 0.5):n_rep])

    return(list(m0 = m0,
                x_hat = mu_x_hat))
  })

  x_star <- list_c$x_hat
  m0 <- list_c$m0
  if (x_star <= 0) return(NA)

  ## get estimate for population regulation
  ## partial derivative around X* (community equilibrium)
  data$data_i <- with(data, {
    sigma_i <<- sd(data_i$n0)
    sigma_c <<- sd(data_i$nt0)
    mu_c <<- sd(data_i$nt0)
    mu_i <<- sd(data_i$nt0)

    data_i <- data_i %>%
      mutate(scl_n0 = (n0 - mu_i) / sigma_i,
             scl_nt0 = (nt0 - mu_c) / sigma_c)

    return(data_i)
  })

  v_rmse2 <- with(data, {
    sapply(theta,
           function(a) loocv(data_i,
                             theta = a,
                             y = "log_r",
                             x = c("n0", "nt0"),
                             random = "species",
                             model = "lmer"))
  })

  theta2 <- theta[which.min(v_rmse2)]

  m <- with(data, {
    d1 <- abs(data_i$nt0 - x_star)
    d2 <- abs(data_i$n0 - 0)
    d <- sqrt(d1^2 + d2^2)

    if (method == "jags") {
      m <- jagslmer(log_r ~ scl_n0 + scl_nt0,
                    random = "species",
                    data = data_i,
                    distance = d / mean(d))
    } else {

      w0 <- exp(-theta2 * (d / mean(d)))
      w <- w0 / sum(w0)
      m <- lme4::lmer(log_r ~ scl_n0 + scl_nt0 + (1 | species),
                      data = data_i,
                      weights = w)
    }

    return(m)
  })

  ## - delta, species level
  if (method == "jags") {
    cnm <- colnames(m$mcmc)
    v_delta <- m$mcmc[, str_detect(cnm, "beta\\[1\\]|beta\\[2\\]|beta\\[3\\]")] %>%
      apply(MARGIN = 1, FUN = function(b) {
        (b[1] - (b[2] / sigma_i) * mu_i - (b[3] / sigma_c) * mu_c) +
          (b[3] / sigma_c) * x_star
      })

  } else {
    msim <- MASS::mvrnorm(n = n_sim,
                          mu = m@beta,
                          Sigma = vcov(m))
    cnm <- colnames(msim)
    v_delta <- msim[, str_detect(cnm, "Intercept|n0|nt0")] %>%
      apply(MARGIN = 1, FUN = function(b) {
        (b[1] - (b[2] / sigma_i) * mu_i - (b[3] / sigma_c) * mu_c) +
          (b[3] / sigma_c) * x_star
      })
  }

  ## competitive exceedance psi
  psi <- mean(v_delta < 0)

  if (method == "jags") {
    attr(psi, "theta") <- c(theta1,
                            with(m, post$`50%`[rownames(post) == "theta"]))
  } else {
    attr(psi, "message") <- summary(m)$optinfo$conv$lme4$messages
    attr(psi, "theta") <- c(theta1, theta2)
  }

  return(psi)
}
