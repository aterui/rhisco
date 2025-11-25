#' Generate control parameters for value generation
#'
#' Creates a list of control parameters used when generating random or constant
#' values in functions such as `set_coef` or `set_r`. This helper function
#' provides sensible defaults while allowing users to customize the distribution
#' type and its associated parameters.
#'
#' @param dist Character. The distribution from which values will be generated.
#'   Options include:
#'   \describe{
#'     \item{constant}{Use a fixed value equal to `mu`.}
#'     \item{exp}{Draw values from an exponential distribution with rate = 1 / `mu`.}
#'     \item{unif}{Draw values from a uniform distribution on `min` and `max`.}
#'     \item{normal}{Draw values from a normal distribution with mean `mu` and
#'                   standard deviation `sd`.}
#'   }
#'
#' @param pg Character. Pattern of propagule introduction. Options are:
#'   \describe{
#'     \item{sync}{Introduce propagules synchronously across all species.}
#'     \item{ord}{Introduce propagules in order.}
#'   }
#'
#' @param mu Numeric. Mean or fixed value. Used for `constant`, `exp`, and
#'   `normal` types. Default is 1.
#'
#' @param sd Numeric. Standard deviation for the `normal` distribution.
#'   Default is 0.1.
#'
#' @param min Numeric. Minimum value for the `unif` distribution. Default is 0.
#' @param max Numeric. Maximum value for the `unif` distribution. Default is 1.
#'
#' @param intv Numeric. Interval or step size used in certain functions
#'   (if applicable). Default is 1.
#'
#' @param seed Integer. Set seed for random number generation.
#'
#' @return A list containing the control settings:
#'   `dist`, `pg`, `mu`, `sd`, `min`, `max`, and `intv`.
#'
#' @export

par.control <- function(dist = c("constant", "exp", "unif", "normal"),
                        pg = c("sync", "ord"),
                        mu = 1,
                        sd = 0.1,
                        min = 0,
                        max = 1,
                        intv = 1,
                        seed = NA) {

  dist <- match.arg(dist)
  pg   <- match.arg(pg)

  list(
    dist = as.character(dist),
    pg   = as.character(pg),
    mu   = mu,
    sd   = sd,
    min  = min,
    max  = max,
    intv = intv,
    seed = seed
  )
}

#' Generate Propagule Introduction Schedule for Species
#'
#' Constructs a matrix \code{m_pg} that specifies when each species receives
#' non-zero propagule input during a warmup period. This function supports either
#' synchronous introduction (all species introduced together at regular intervals)
#' or sequential introduction (species introduced one at a time in randomized
#' order, separated by a fixed interval).
#'
#' The output is a matrix with \code{n_sim} rows (time steps) and \code{s}
#' columns (species). Values are zero before a species is introduced. When a
#' species is introduced, its propagule value at that time step is drawn from a
#' uniform distribution defined by \code{pg$min} and \code{pg$max}.
#'
#' @param warmup Integer. Number of warmup time steps during which propagule
#' introduction occurs.
#' @param n_sim Integer. Total number of simulation time steps.
#' @param s Integer. Number of species.
#' @param pg A list of control settings created with \code{\link{par.control}}.
#'   Defines introduction type, bounds of the uniform distribution, and the
#'   introduction interval.
#'
#' @details
#' Two introduction modes are available:
#'
#' \describe{
#'   \item{\code{pg = "sync"}}{
#'   All species receive propagules simultaneously at regular intervals
#'   determined by \code{intv}.}
#'
#'   \item{\code{pg = "ord"}}{
#'   Species receive propagules one by one in a random order. The interval
#'   between introductions is controlled by \code{intv}. If the final
#'   introduction time exceeds \code{warmup}, remaining species will not be
#'   introduced and a warning is issued.}
#' }
#'
#' @return
#' A numeric matrix of dimension \code{n_sim x s}. Entries are zero until a
#' species is introduced, at which point propagule values are assigned at the
#' introduction time step.
#'
#' @examples
#' \dontrun{
#'
#' ## Synchronous propagule introduction
#' set_pg(warmup = 100, n_sim = 150, s = 4,
#'        pg = par.control(pg = "sync"))
#'
#' ## Sequential propagule introduction every 5 steps
#' set_pg(warmup = 100, n_sim = 150, s = 4,
#'        pg = par.control(pg = "ord", min = 0.2, max = 0.8, intv = 5))
#'
#' }
#'
#' @export

set_pg <- function(warmup,
                   n_sim,
                   s,
                   pg = par.control()) {

  ## define resample
  resample <- function(x, ...) x[sample.int(length(x), ...)]

  if (!is.list(pg))
    stop("'pg' must be a list.")

  par_pg <- utils::modifyList(par.control(), pg)
  m_pg <- matrix(0,
                 nrow = n_sim,
                 ncol = s)

  with(par_pg, {

    if (!is.na(seed)) set.seed(seed)

    switch(pg,
           sync = {
             ## time index for introduction
             idx <- seq(from = 1,
                        to = warmup,
                        by = intv)

             m_pg[idx, ] <- stats::runif(n = s * length(idx),
                                         min = min,
                                         max = max)

             return(m_pg)
           },
           ord = {
             ## time index for introduction
             idx <- seq(from = 1,
                        by = intv,
                        length.out = s)

             if (max(idx) > warmup) {
               idx <- idx[idx < warmup]

               warning("Introduction times exceed 'warmup'. Species after'warmup' time limit will not be introduced.")
             }

             ## random species introduction order
             intro_order <- resample(seq_len(s), size = s)
             intro_order <- intro_order[seq_len(length(idx))]

             m_pg[cbind(idx, intro_order)] <- stats::runif(n = length(idx),
                                                           min = min,
                                                           max = max)

             return(m_pg)
           },
           stop("Unsupported pg type."))
  })
}

#' Generate Diagonal and Interaction Coefficients
#'
#' Creates an \eqn{s \times s} coefficient matrix whose diagonal elements
#' represent self-regulation terms and whose off-diagonal elements represent
#' interaction strengths. The distribution and parameters used to generate these
#' coefficients are controlled via `alpha` and `beta`, which are
#' lists produced by [par.control()].
#'
#' The diagonal elements are generated from `alpha`, while the full
#' matrix (before overwriting the diagonal) is generated from `beta`.
#' Each control list specifies the type of distribution (`"constant"`, `"exp"`,
#' `"unif"`, or `"normal"`) and its associated parameters.
#'
#' @param s Integer. The dimension of the coefficient matrix.
#' @param alpha A list generated by [par.control()] that specifies how
#'   to generate the \eqn{s} diagonal coefficients. User-supplied values override
#'   defaults through partial matching via [utils::modifyList()].
#' @param beta A list generated by [par.control()] that specifies how
#'   to generate the \eqn{s \times s} matrix of off-diagonal coefficients.
#'   User-supplied values override defaults through [utils::modifyList()].
#' @param negative Logical. If `TRUE` (default), all coefficients are multiplied
#'   by `-1`, making them negative. If `FALSE`, coefficients remain positive.
#'
#' @details
#' For `dist = "constant"`, `mu` may be either a scalar or a vector of length
#' \eqn{s} for `alpha`, and either a scalar or an \eqn{s \times s}
#' matrix for `beta`. For all other types, random draws are generated
#' from the appropriate distribution.
#'
#' The diagonal of the resulting matrix is always replaced by the values
#' generated from `alpha`, ensuring that diagonal and off-diagonal
#' coefficients can be specified independently. The `negative` argument allows
#' the user to control the sign of the coefficients.
#'
#' @return
#' A numeric \eqn{s \times s} matrix of coefficients.
#'
#' @examples
#' \dontrun{
#'
#' # Default: negative constant self-effects and interaction effects
#' set_coef(s = 4)
#'
#' # Uniform variation in diagonal coefficients
#' set_coef(s = 4, alpha = list(dist = "unif", min = 0.5, max = 2))
#'
#' # Normal variation in off-diagonal interaction strengths
#' set_coef(s = 4, beta = list(dist = "normal", mu = 0, sd = 0.2))
#'
#' # Keep coefficients positive
#' set_coef(s = 4, negative = FALSE)
#'
#' }
#'
#' @export

set_coef <- function(s,
                     alpha = par.control(),
                     beta = par.control(),
                     negative = TRUE) {

  # Merge user values with defaults
  if (!is.list(alpha))
    stop("'alpha' must be a list.")
  if (!is.list(beta))
    stop("'beta' must be a list.")

  par_alpha <- utils::modifyList(par.control(), alpha)
  par_beta <- utils::modifyList(par.control(), beta)

  v_alpha <- with(par_alpha, {

    if (!is.na(seed)) set.seed(seed)

    switch(dist,

           constant = {
             if (length(mu) == 1) {
               rep(mu, times = s)
             } else {
               if (length(mu) != s)
                 stop("'mu' in 'alpha' must be scalar or vector of length ", s)
               mu
             }
           },

           exp = {
             if (mu <= 0) stop("'mu' must be positive for 'exp' type")

             stats::rexp(n = s,
                         rate = 1 / mu)
           },

           unif = {
             stats::runif(n = s,
                          min = min,
                          max = max)
           },

           normal = {
             stats::rnorm(n = s,
                          mean = mu,
                          sd = sd)
           },

           stop("Unsupported dist in 'alpha'.")
    )
  })

  m_coef <- with(par_beta, {

    if (!is.na(seed)) set.seed(seed)

    switch(dist,

           constant = {
             if (length(mu) == 1) {
               matrix(mu,
                      nrow = s,
                      ncol = s)
             } else {
               if (!is.matrix(mu) || any(dim(mu) != c(s, s)))
                 stop("'mu' in 'beta' must be scalar or s x s matrix when dist = 'constant'.")

               mu
             }
           },

           exp = {
             if (mu <= 0) stop("'mu' must be positive for 'exp' type")
             matrix(stats::rexp(n = s * s,
                                rate = 1 / mu),
                    nrow = s,
                    ncol = s)
           },

           unif = {
             matrix(stats::runif(n = s * s,
                                 min = min,
                                 max = max),
                    nrow = s,
                    ncol = s)
           },

           normal = {
             matrix(stats::rnorm(n = s * s,
                                 mean = mu,
                                 sd = sd),
                    nrow = s,
                    ncol = s)
           },

           stop("Unsupported dist in 'beta'.")
    )
  })

  diag(m_coef) <- v_alpha

  m_coef <- negative * (-m_coef) + (1 - negative) * m_coef

  return(m_coef)
}

#' Generate Species Growth Rates
#'
#' Generates a numeric vector of intrinsic growth rates (\eqn{r}) for a set of
#' species. The length of the vector is specified by `s`, and the distribution
#' and parameters used to generate these rates are controlled via `r`,
#' which is a list produced by [par.control()].
#'
#' @param s Integer. The number of species, i.e., the length of the resulting
#'   growth rate vector.
#' @param r A list generated by [par.control()] that specifies how to
#'   generate the vector of growth rates. User-supplied values override defaults
#'   through partial matching via [utils::modifyList()].
#'
#' @details
#' The `r` list specifies the type of distribution (`"constant"`,
#' `"exp"`, `"unif"`, or `"normal"`) and its associated parameters (`mu`, `sd`,
#' `min`, `max`). For `dist = "constant"`, `mu` may be a scalar or a vector of
#' length \eqn{s}. For other types, random draws are generated from the specified
#' distribution. For `exp` type, `mu` must be positive.
#'
#' @return
#' A numeric vector of length \eqn{s} containing the growth rates.
#'
#' @examples
#' \dontrun{
#' # Default: constant growth rates
#' set_r(s = 5)
#'
#' # Exponentially distributed growth rates
#' set_r(s = 5, r = list(dist = "exp", mu = 0.5))
#'
#' # Uniformly distributed growth rates
#' set_r(s = 5, r = list(dist = "unif", min = 0, max = 1))
#'
#' # Normally distributed growth rates
#' set_r(s = 5, r = list(dist = "normal", mu = 0.5, sd = 0.1))
#' }
#'
#' @export

set_r <- function(s,
                  r = par.control()) {

  # Ensure r is a list
  if (!is.list(r))
    stop("'r' must be a list.")

  # Merge user-supplied control with defaults
  par_r <- utils::modifyList(par.control(), r)

  # Generate vector of r values
  v_r <- with(par_r, {

    if (!is.na(seed)) set.seed(seed)

    switch(dist,

           constant = {
             if (length(mu) == 1) {
               rep(mu, times = s)
             } else {
               if (length(mu) != s)
                 stop("'mu' in 'r' must be scalar or vector of length ", s)
               mu
             }
           },

           exp = {
             if (mu <= 0) stop("'mu' must be positive for 'exp' type")

             stats::rexp(n = s,
                         rate = 1 / mu)
           },

           unif = {
             stats::runif(n = s,
                          min = min,
                          max = max)
           },

           normal = {
             stats::rnorm(n = s,
                          mean = mu,
                          sd = sd)
           },

           stop("Unsupported dist in 'r'.")
    )
  })

  return(v_r)
}

#' Population Dynamics Function
#'
#' This function returns an update function corresponding to a chosen
#' population dynamics model. The returned function takes current abundances and
#' model parameters to compute the next time-step abundances, with an option for
#' demographic stochasticity. The function currently supports two models:
#' \itemize{
#'   \item \code{"bh"}: Beverton–Holt model with density dependence.
#'   \item \code{"ricker"}: Ricker model with density dependence.
#' }
#'
#' Both models incorporate a vector of environmental stochastic effects
#' (\code{eps}) and a density-dependence matrix (\code{m_coef}).
#' When \code{stochastic = TRUE}, abundances are drawn
#' from a Poisson distribution; otherwise, expected values are returned.
#'
#' @param model A character string specifying the model to use. Supported values
#'   are \code{"bh"} for Beverton–Holt and \code{"ricker"} for the Ricker model.
#'
#' @return A function that computes one time-step update of population
#' abundances. The returned function has arguments:
#' \code{r}, \code{n}, \code{m_coef}, \code{eps}, and \code{stochastic}.
#'
#' @export

fn_model <- function(model) {

  dyn <- switch(model,
                bh = function(r,
                              n,
                              m_coef,
                              eps,
                              stochastic = FALSE) {

                  n_bar <- (n * exp(r)) / (1 + m_coef %*% n)
                  n_hat <- n_bar * exp(eps)

                  if (stochastic) {
                    y <- stats::rpois(n = length(n_hat),
                                      lambda = n_hat)
                  } else {
                    y <- drop(n_hat)
                  }

                  return(y)
                },

                ricker = function(r,
                                  n,
                                  m_coef,
                                  eps,
                                  stochastic = FALSE) {

                  n_bar <- n * exp(r + m_coef %*% n)
                  n_hat <- n_bar * exp(eps)

                  if (stochastic) {
                    y <- stats::rpois(n = length(n_hat),
                                      lambda = n_hat)
                  } else {
                    y <- drop(n_hat)
                  }

                  return(y)
                },

                stop("Unsupported model type")
  )

  return(dyn)
}

#' Simulate Community Dynamics
#'
#' Simulates multispecies community dynamics using either the Ricker or
#' Beverton-Holt population model, with species interactions and
#' environmental stochasticity. Species can be initialized and introduced
#' according to a user-specified propagule schedule during warmup.
#'
#' The model runs for `warmup + burnin + ts` timesteps, and returns species
#' densities in a tidy data format. Warmup and burn-in periods can be excluded
#' from the output using `trim = TRUE`.
#'
#' @param ts Integer. Number of timesteps to retain after warmup and burn-in.
#' @param warmup Integer. Number of warmup steps with propagule introduction (excluded from output if `trim = TRUE`).
#' @param burnin Integer. Number of burn-in steps without propagule introduction (excluded from output if `trim = TRUE`).
#' @param s Integer. Number of species in the community.
#' @param r List or scalar passed to `par.control()`. Controls intrinsic growth rates.
#' @param alpha List or scalar passed to `par.control()`. Controls competition coefficients.
#' @param beta List or scalar passed to `par.control()`. Controls facilitation coefficients.
#' @param model Character. Population model: `"ricker"` or `"bh"`.
#' @param negative Logical. Whether interaction coefficients are set negative by default.
#'   If `NULL`, defaults to `TRUE` for Ricker and `FALSE` for Beverton-Holt.
#' @param sd_eps Numeric. Standard deviation of environmental stochasticity.
#' @param stochastic Logical. If `TRUE`, applies environmental stochasticity to growth.
#' @param pg Object passed to `par.control()`. Defines the propagule introduction schedule.
#' @param extinct Numeric. Density threshold below which species density is set to zero.
#' @param trim Logical. If `TRUE` (default), returns only post warmup+burn-in timesteps.
#'
#' @return A tibble with columns:
#'   * `ts` timestep
#'   * `species` species ID
#'   * `density` species density
#'
#' Returned object includes attributes:
#'   * `"r"` intrinsic growth parameters
#'   * `"coef"` species interaction matrix
#'
#' @export

csim <- function(ts = 1000,
                 warmup = 100,
                 burnin = 100,
                 s = 10,
                 r = par.control(),
                 alpha = par.control(),
                 beta = par.control(),
                 model = c("ricker", "bh"),
                 negative = NULL,
                 sd_eps = 0.1,
                 stochastic = FALSE,
                 pg = par.control(),
                 extinct = 0,
                 trim = TRUE) {

  # model type --------------------------------------------------------------

  model <- match.arg(model)
  if (is.null(negative))
    negative <- ifelse(model == "ricker", TRUE, FALSE)

  dyn <- fn_model(model = model)

  # parms -------------------------------------------------------------------

  # number of total simulation steps and final trimming index
  n_sim <- warmup + burnin + ts
  n_cut <- warmup + burnin

  ## parameter: intrinsic growth
  v_r <- set_r(s = s,
               r = r)

  ## parameter: species interaction
  m_coef <- set_coef(s = s,
                     alpha = alpha,
                     beta = beta,
                     negative = negative)

  # variables ---------------------------------------------------------------

  cnm <- c("ts",
           "species",
           "density")

  m_dyn <- matrix(NA,
                  nrow = n_sim * s,
                  ncol = length(cnm))

  colnames(m_dyn) <- cnm

  st_row <- seq(from = 1,
                to = nrow(m_dyn),
                by = s)

  # propagule introduction schedule (time x species)
  m_pg <- set_pg(warmup = warmup,
                 n_sim = n_sim,
                 s = s,
                 pg = pg)

  ## environmental stochasticity
  m_eps <- matrix(stats::rnorm(n = n_sim * s,
                               mean = 0,
                               sd = sd_eps),
                  nrow = n_sim,
                  ncol = s,
                  byrow = TRUE) # time x species matrix

  # dynamics ----------------------------------------------------------------

  ## initialize m_dyn[,]
  m_dyn[1:s, ] <- cbind(rep(1, s), # time step
                        seq_len(s), # species ID
                        m_pg[1, ] # initial density
  )

  ## initialize density
  v_n <- rep(0, times = s)

  for (i in 1:(n_sim - 1)) {

    ## add propagule
    v_n <- v_n + m_pg[i, ]

    ## update density
    v_n <- dyn(r = v_r,
               n = v_n,
               m_coef = m_coef,
               eps = m_eps[i, ],
               stochastic = stochastic)

    v_n[v_n < extinct] <- 0

    row_id <- seq(from = st_row[i + 1],
                  to = st_row[i + 1] + (s - 1),
                  by = 1)

    m_dyn[row_id, ] <- cbind(rep(i + 1, s), # timestep
                             seq_len(s), # species ID
                             v_n # density
    )

  }

  # return ------------------------------------------------------------------

  if (trim) m_dyn <- m_dyn[(n_cut * s + 1):(n_sim * s), ]

  df_dyn <- dplyr::as_tibble(m_dyn)

  attr(df_dyn, "r") <- v_r
  attr(df_dyn, "coef") <- m_coef

  return(df_dyn)

}
