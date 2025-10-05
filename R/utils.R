#' @param m Integer. Magnitude of a given link.
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
                              "lm" = lm(formula,
                                        data = df_train,
                                        weights = w,
                                        ...),
                              "glm" = glm(formula,
                                          data = df_train,
                                          weights = w,
                                          ...),
                              "lmer" = lme4::lmer(formula,
                                                  data = df_train,
                                                  weights = w,
                                                  ...),
                              "glmer" = lme4::glmer(formula,
                                                    data = df_train,
                                                    weights = w,
                                                    ...),
                              stop("Unsupported model type")
                 )

                 y0 <- predict(lw, newdata = data[i, ])
                 y1 <- data[i, y]

                 eps <- (y1 - y0)^2
                 return(eps)
               }) |>
    sum()

  return(sqrt(ss / nrow(m_x)))
}
