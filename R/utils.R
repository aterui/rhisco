#'
#'
#' @param m Integer. Magnitude of a given link.
#'
#' @importFrom stats dpois qpois
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

loocv <- function(data,
                  theta,
                  x = "x",
                  y = "y",
                  model = "lm",
                  random = NULL,
                  maxit = 100) {

  ## select relevant columns
  df0 <- data %>%
    mutate(id = row_number()) %>%
    mutate(across(.cols = which(colnames(.) %in% x),
                  .fns = function(x) c(scale(x))))

  ## get matrix x & y
  mxy <- df0[, c(x, y)] %>%
    data.matrix()

  ## distance matrix for weighting
  dmat <- mxy %>%
    dist(diag = TRUE,
         upper = TRUE) %>%
    data.matrix()

  ## summed squared deviance
  if (model == "lm") {
    ss <- sapply(1:(nrow(mxy)),
                 function(i) {
                   ## training data
                   df_train <- df0[-i, ] %>%
                     dplyr::select(all_of(c(y, x)))

                   ## calculate weights
                   mu_d <- mean(dmat[i, -i])
                   w <- exp(- theta * dmat[i, -i] / mu_d)

                   ## prediction
                   f <- parse(text = paste0("lm(", y, " ~ ",
                                            paste(x, collapse = " + "), ", ",
                                            "data = df_train, ",
                                            "weights = w)"))
                   lw <- eval(f)

                   y1 <- mxy[i, y]
                   y0 <- drop(c(1, mxy[i, x]) %*% coef(lw))

                   eps <- (y1 - y0)^2
                   return(eps)
                 }) %>%
      sum()
  }

  if (model == "rlm") {
    ss <- sapply(1:(nrow(mxy)),
                 function(i) {
                   ## training data
                   df_train <- df0[-i, ] %>%
                     dplyr::select(all_of(c(y, x)))

                   ## calculate weights
                   mu_d <- mean(dmat[i, -i])
                   w <- exp(- theta * dmat[i, -i] / mu_d)

                   ## prediction
                   f <- parse(text = paste0("MASS::rlm(", y, " ~ ",
                                            paste(x, collapse = " + "), ", ",
                                            "data = df_train, ",
                                            "maxit = maxit, ",
                                            "weights = w)"))
                   lw <- eval(f)

                   y1 <- mxy[i, y]
                   y0 <- drop(c(1, mxy[i, x]) %*% coef(lw))

                   eps <- (y1 - y0)^2
                   return(eps)
                 }) %>%
      sum()
  }

  if (model == "lmer") {
    if (is.null(random))
      stop("'random' is required for 'lmer'")

    ss <- sapply(1:(nrow(mxy)),
                 function(i) {
                   ## training data
                   df_train <- df0[-i, ] %>%
                     dplyr::select(all_of(c(y, x, random)))

                   ## calculate weights
                   mu_d <- mean(dmat[i, -i])
                   w <- exp(- theta * dmat[i, -i] / mu_d)

                   ## prediction
                   f <- parse(text = paste0("lme4::lmer(", y, " ~ ",
                                            paste(x, collapse = " + "),
                                            paste0(" + ", "(1 |", random, ")"),
                                            ", ",
                                            "data = df_train, ",
                                            "control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = 'ignore',  tol = 1e-4)), ",
                                            "weights = w)"))
                   lw <- eval(f)

                   y1 <- mxy[i, y]
                   y0 <- drop(c(1, mxy[i, x]) %*% lw@beta)

                   eps <- (y1 - y0)^2
                   return(eps)
                 }) %>%
      sum()
  }

  return(sqrt(ss / nrow(mxy)))
}
