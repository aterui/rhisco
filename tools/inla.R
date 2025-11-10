#' Fit Linear Mixed Models Using INLA with lmer-Style Formulas
#'
#' This function allows you to fit linear mixed models using the
#' \code{INLA} package while specifying formulas in the familiar
#' \code{lme4::lmer} style (e.g., \code{y ~ x + (1|group)}). Random
#' intercepts are automatically converted to \code{INLA}'s
#' \code{f(group, model='iid')} format.
#'
#' @param formula A formula object in \code{lmer} style, e.g.,
#'   \code{y ~ x + (1|group1) + (1|group2)}.
#' @param data A data frame containing all variables used in the
#'   formula.
#' @param family The likelihood family for \code{INLA}. Defaults to
#'   \code{"gaussian"}.
#' @param control.compute \code{control.compute} argument in \code{INLA::inla()}.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
#'
#' @return An \code{INLA} model object returned by
#'   \code{INLA::inla()}, including posterior summaries and
#'   fitted values.
#'
#' @details
#' This wrapper function handles:
#' \itemize{
#'   \item Conversion of fixed effects and random intercepts
#'     from \code{lmer}-style formulas to \code{INLA} format.
#'   \item Multiple random intercepts.
#' }
#'
#' Limitations:
#' \itemize{
#'   \item Does not currently support random slopes (e.g., \code{(x|group)}).
#'   \item Requires \code{lme4} and \code{INLA} packages.
#' }
#'
#' @examples
#' \dontrun{
#' library(INLA)
#' library(lme4)
#'
#' # Simulate example data
#' set.seed(123)
#' dat <- data.frame(
#'   y = rnorm(20),
#'   x = rnorm(20),
#'   group = rep(1:4, each = 5)
#' )
#'
#' # Fit model using wrapper
#' fit <- inla_lmer(y ~ x + (1|group), data = dat)
#' summary(fit)
#' }
#'
#' @export

inla_lmer <- function(formula,
                      data,
                      family = "gaussian",
                      control.compute = list(config = TRUE),
                      ...) {

  # Convert lmer-style random effects to INLA format
  # Example: y ~ x + (1|group) -> y ~ x + f(group, model="iid")

  # Extract terms
  tt <- stats::terms(formula, keep.order = TRUE)
  v_bars <- lme4::findbars(formula)

  res <- attr(tt, "variables")[[2]]
  allt <- attr(tt, "term.labels")
  fixt <- allt[!grepl("\\|", allt)]
  fixt_inla <- paste(fixt, collapse = " + ")

  # Convert random effects
  formula_char <- paste(res, "~",
                        fixt_inla)

  if (!is.null(v_bars)) {
    rant <- v_bars |>
      sapply(FUN = function(x) as.character(x[[3]]))

    rant_inla <- paste0("f(", rant, ", model = 'iid')",
                        collapse = " + ")

    formula_char <- paste(formula_char,
                          " + ",
                          rant_inla)
  }

  formula_inla <- stats::as.formula(formula_char)

  # Run INLA
  res <- INLA::inla(formula_inla,
                    data = data,
                    family = family,
                    control.compute = control.compute,
                    ...)

  return(res)
}

#' Extract Posterior Samples of Fixed Effects from an INLA Model
#'
#' This function extracts posterior samples for all fixed-effect coefficients
#' from a fitted \code{INLA} model object. It provides a convenient interface
#' for accessing parameter uncertainty in a format similar to
#' \code{fixef()} in mixed-effects modeling frameworks.
#'
#' @param m A fitted \code{INLA} model object, typically returned by
#'   \code{INLA::inla()} or a wrapper such as \code{inla_lmer()}.
#' @param n Integer. The number of posterior samples to draw from the fitted
#'   model. Default is \code{1000}.
#'
#' @return A numeric matrix of posterior samples with \code{n} rows (samples)
#'   and one column per fixed-effect parameter. Column names correspond to
#'   the fixed-effect names from the model.
#'
#' @details
#' The function uses \code{INLA::inla.posterior.sample()} to draw samples
#' from the joint posterior distribution and then evaluates each sample’s
#' fixed-effect component using
#' \code{INLA::inla.posterior.sample.eval()}.
#'
#' This allows users to summarize uncertainty in regression coefficients,
#' visualize posterior densities, or propagate parameter uncertainty in
#' downstream predictions.
#'
#' @examples
#' \dontrun{
#' library(INLA)
#'
#' # Example model
#' set.seed(123)
#' dat <- data.frame(
#'   y = rnorm(20),
#'   x = rnorm(20),
#'   group = rep(1:4, each = 5)
#' )
#' fit <- inla(y ~ x + f(group, model = "iid"),
#'             data = dat,
#'             family = "gaussian")
#'
#' # Extract 500 posterior samples of fixed effects
#' post <- fixef_posterior(fit, n = 500)
#' }
#'
#' @export

fixef_posterior <- function(m,
                            n = 1000) {

  if (!inherits(m, what = "inla"))
    stop("m must be class 'inla'")

  fixt <- m$names.fixed
  post <- INLA::inla.posterior.sample(n = n,
                                      result = m)

  post_sample <- INLA::inla.posterior.sample.eval(fun = fixt,
                                                  samples = post) |>
    t()

  colnames(post_sample) <- fixt
  rownames(post_sample) <- NULL

  return(post_sample)
}

#' Point predictions from INLA model fits
#'
#' Generate linear predictor estimates (posterior means) for fitted
#' \code{INLA} model objects, optionally using new data.
#'
#' @param m A fitted \code{INLA} model object.
#' @param newdata Optional data frame containing new observations
#'   for which predictions are desired. Must contain all covariates
#'   and random-effect grouping factors used in the model formula.
#'   If \code{NULL}, predictions are computed for the data used to
#'   fit the model.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric vector of posterior mean predictions for the linear
#' predictor (\eqn{\eta = X\beta + Zu}). Predictions are made at the
#' population level for unseen random-effect levels (i.e., random effects
#' not present in the training data are assumed to have mean zero).
#'
#' @details
#' This function reconstructs the linear predictor by combining
#' posterior means of the fixed and random effects from a fitted INLA
#' model. Random effects corresponding to levels not present in the
#' training data are set to zero, representing the population-level mean.
#'
#' Note that this function provides point estimates only; it does not
#' propagate posterior uncertainty. For predictive uncertainty or
#' posterior predictive distributions, use
#' \code{\link[INLA]{inla.posterior.sample}} and compute predictions
#' from sampled parameters.
#'
#' @seealso
#' \code{\link[INLA]{inla}}, \code{\link[INLA]{inla.posterior.sample}}
#'
#' @examples
#' \dontrun{
#' library(INLA)
#'
#' set.seed(123)
#' dat <- data.frame(
#'   y = rnorm(20),
#'   x = rnorm(20),
#'   group = rep(1:4, each = 5)
#' )
#'
#' fit <- inla(y ~ x + f(group, model = "iid"), data = dat)
#'
#' # Predictions for existing data
#' head(point_predict(fit))
#'
#' # Predictions for new data (group 5 unseen → population mean)
#' newdat <- data.frame(x = c(-1, 0, 1), group = c(1, 2, 5))
#' point_predict(fit, newdata = newdat)
#' }
#'
#' @export

point_predict <- function(m,
                          newdata = NULL,
                          ...) {

  # if no newdata, use model data
  if (is.null(newdata)) {
    newdata <- m$.args$data
  }

  # fixed effect
  v_b <- m$summary.fixed[["mean"]]
  v_cnm <- m$names.fixed
  fix_form <- paste("~",
                    paste(v_cnm[!grepl("[Ii]ntercept", v_cnm)],
                          collapse = " + ")) |>
    stats::as.formula()

  X <- stats::model.matrix(fix_form,
                           data = newdata)

  # random effect
  ranef_names <- names(m$summary.random)
  re_list <- list()

  if (length(ranef_names) > 0) {

    # loop across different random effects
    for (rn in ranef_names) {
      re <- m$summary.random[[rn]]
      idx <- match(newdata[[rn]], re$ID)
      idx[is.na(idx)] <- 0  # unseen levels → NA/0
      re_mean <- ifelse(idx > 0, re$mean[idx], 0)
      re_list[[rn]] <- re_mean
    }

    Z <- Reduce(`+`, re_list)
  } else {
    Z <- 0
  }

  # linear predictor for each sample
  v_eta <- drop(X %*% v_b + Z)
  names(v_eta) <- NULL

  return(v_eta)
}

