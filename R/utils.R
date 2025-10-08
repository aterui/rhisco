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
  if (k >= 0) c(rep(NA, k), head(x, -k))
  else c(tail(x, k), rep(NA, -k))
}

