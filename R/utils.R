#' Generate Lagged or Lead Versions of a Vector
#'
#' This function creates a lagged (or lead) version of a numeric or character vector by shifting its elements
#' forward or backward by a specified number of positions.
#' Missing values (\code{NA}) are used to fill in positions that fall outside the original range after shifting.
#'
#' @param x A numeric, character, or logical vector. The vector to be lagged or led.
#' @param k An integer scalar specifying the number of positions to shift.
#'   - A positive value of \code{k} produces a **lag**, moving elements forward and filling the first \code{k} elements with \code{NA}.
#'   - A negative value of \code{k} produces a **lead**, moving elements backward and filling the last \code{abs(k)} elements with \code{NA}.
#'
#' @return A vector of the same length as \code{x}, with elements shifted according to \code{k}.
#'
#' @examples
#' x <- 1:5
#' lag_base(x, k = 1)   # returns c(NA, 1, 2, 3, 4)
#' lag_base(x, k = -1)  # returns c(2, 3, 4, 5, NA)
#' lag_base(x, k = 0)   # returns c(1, 2, 3, 4, 5)
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

lag_base <- function(x, k = 1L) {
  if (k >= 0) c(rep(NA, k), utils::head(x, -k))
  else c(utils::tail(x, k), rep(NA, -k))
}

#' Generate Lagged Values Within Groups
#'
#' This function computes lagged values of a specified variable within each group
#' defined by a grouping (index) column. It first orders data within each group by
#' a time or sequence variable, then applies a one-step lag to the target variable.
#' The resulting data frame includes a new column containing lagged values.
#'
#' @param data A data frame containing at least the grouping, time, and response variables.
#' @param y Character string specifying the name of the response variable to be lagged.
#'   Default is \code{"n"}.
#' @param t Character string specifying the name of the time or sequence variable
#'   used to order observations within each group. Default is \code{"t"}.
#' @param index Character string specifying the grouping variable (e.g., species, site, ID).
#'   Default is \code{"species"}.
#' @param total Logical. If \code{TRUE}, total community density will be calculated.
#'
#' @return A data frame identical to \code{data} but with an additional column named
#'   \code{<y>_lag}, containing the one-step lagged values of the response variable
#'   computed within each group.
#'
#' @examples
#' data <- data.frame(
#'   species = rep(c("A", "B"), each = 4),
#'   t = rep(1:4, 2),
#'   n = c(5, 6, 7, 8, 2, 3, 4, 5)
#' )
#'
#' lag_block(data, index = "species")
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @seealso \code{\link{lag_base}}, \code{\link{split}}
#'
#' @export

lag_block <- function(data,
                      y = "n",
                      t = "t",
                      index = NULL,
                      total = TRUE) {

  if (is.null(index)) {

    data <- data[order(data[[t]]), ]
    data[[paste0(y, "_lag")]] <- lag_base(data[[y]], k = 1)
    df_lag <- data

    if (any(table(data[[t]]) > 1))
      stop("Detected multiple records for the same time point. If you intended to create blocks separately for different groups, please specify the `index` argument.")

  } else {

    if (total) {
      nt <- paste0(y, "t")
      nt_lag <- paste0(y, "t_lag")

      data[[nt]] <- stats::ave(data[[y]],
                               data[[t]],
                               FUN = sum)
    }

    list_lag <- split(x = data,
                      f = data[[index]]) |>
      lapply(function(d) {
        d <- d[order(d[[t]]), ]  # sort by time column
        d[[paste0(y, "_lag")]] <- lag_base(d[[y]], k = 1)

        if (total) {
          d[[nt_lag]] <- lag_base(d[[nt]], k = 1)
        }

        return(d)
      })

    df_lag <- do.call(rbind,
                      c(list_lag, make.row.names = FALSE))

  }

  return(df_lag)
}

#' Generate a distance-based weighting function
#'
#' Returns a function that computes weights based on distance, using either
#' an exponential or Gaussian kernel.
#'
#' @param type Character string specifying the kernel type. Options are
#'   `"exp"` (exponential) or `"gaussian"`.
#' @param method Character string specifying the scaling method. Options are
#'   `"mean"` or `"max"`.
#'
#' @return A function of the form \code{f(x, theta)} that computes weights
#'   for a numeric vector \code{x} given parameter \code{theta}.
#'
#' @examples
#' dfun <- get_dfun("exp")
#' dfun(1:5, theta = 0.5)
#'
#' @export

get_dfun <- function(type,
                     method = "max") {

  f <- switch(method,
              max = max,
              mean = mean,
              stop("Unsupported method"))

  switch(type,
         exp = function(x, theta) {
           d <- x / f(x)
           w <- exp(-theta * d)
           return(w)
         },
         gaussian = function(x, theta) {
           d <- x / f(x)
           w <- exp(- theta^2 * (d^2 / 2))
           return(w)
         },
         stop("Unsupported type"))

}

#' Run an expression with an optional local random seed
#'
#' @param seed Numeric or NA. If numeric, sets a local RNG seed for reproducibility.
#'             If NA, runs the expression using the current RNG state.
#' @param fn A function with no arguments that returns the value to generate.
#'
#' @return The result of evaluating \code{fn()} with or without the local seed.
#'
#' @export

gen_values <- function(seed, fn) {
  if (!(is.na(seed) || is.null(seed))) {
    withr::with_seed(seed, fn())
  } else {
    fn()
  }
}

#' Extract grouping factor from a mixed-effects formula
#'
#' Returns the first grouping factor (text after `|`) in a mixed-effects
#' model formula. If no random-effect term is present, `NULL` is returned.
#'
#' @param formula A model formula potentially containing random effects,
#'   e.g. `y ~ x + (1 | species)`.
#'
#' @return A character string giving the grouping factor name,
#'   or `NULL` if no random effect is found.
#'
#' @export

extract_group <- function(formula) {
  re_terms <- reformulas::findbars(formula)

  if (length(re_terms) == 0) return(NULL)

  as.character(re_terms[[1]][[3]])
}

#' Recover the full variance-covariance matrix for unscaled predictors
#'
#' Transforms a variance covariance matrix from a model with centered and scaled predictors
#' back to its original scale, including the intercept and all covariance terms.
#'
#' @param v The variance-covariance matrix from the scaled model, where the first row/column is the intercept.
#' @param means A numeric vector of the means used to center the predictors.
#' @param stds A numeric vector of the standard deviations used to scale the predictors.
#'
#' @return A matrix representing the VCV on the original scale.
#'
#' @export

vcov_unscale <- function(v, means, stds) {
  # Input validation
  if (!is.numeric(means) || !is.numeric(stds))
    stop("means and stds must be numeric vectors")

  if (length(means) != length(stds))
    stop("means and stds must have the same length")

  p <- length(means)

  if (!is.matrix(v) || !isTRUE(all.equal(dim(v), c(p + 1L, p + 1L))))
    stop("v must be a (length(means)+1) x (length(means)+1) matrix")

  if (!isSymmetric(v, tol = 1e-8))
    stop("v must be a symmetric matrix")

  if (any(abs(stds) < .Machine$double.eps * 100))
    stop("stds contains zero or near-zero values; division undefined")

  # Build transformation matrix
  A <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  A[1, 1] <- 1
  A[1, 2:(p + 1)] <- -means / stds
  A[cbind(2:(p+1), 2:(p+1))] <- 1 / stds

  # Transform and force symmetry to absorb floating-point error
  result <- A %*% v %*% t(A)
  (result + t(result)) / 2
}

#' Check Random Effects Structure
#'
#' Validates that a formula contains at most one random intercept and no random
#' slopes. Intended for use in models where random slopes cause scale-dependent
#' results.
#'
#' @param formula A two-sided formula, potentially containing random effects.
#'
#' @export

re_check <- function(formula) {
  re_terms <- reformulas::findbars(formula)

  if (is.null(re_terms)) return(invisible(NULL))

  if (length(re_terms) > 1)
    stop("Only one random effect is allowed, e.g., (1 | group).")

  has_slope <- sapply(re_terms, function(x) {
    slope_terms <- paste("~", deparse(x[[2]])) |>
      stats::as.formula() |>
      stats::terms() |>
      attr("term.labels")

    length(slope_terms) > 0
  })

  if (any(has_slope))
    stop("Random slopes are not allowed. Use random intercept only, e.g., (1 | group).")

  invisible(NULL)
}
