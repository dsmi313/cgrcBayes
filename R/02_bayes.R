# Normal-Inverse-Gamma posterior for each stratum mean.
#
#   y_ij | mu_j, sigma2_j ~ N(mu_j, sigma2_j)
#   mu_j | sigma2_j       ~ N(m0, sigma2_j / k0)      <- NOTE: depends on sigma2
#   sigma2_j              ~ InverseGamma(a0, b0)
#
# Posterior update:
#   k_n = k0 + n
#   m_n = (k0*m0 + n*ybar) / k_n
#   a_n = a0 + n/2
#   b_n = b0 + SS/2 + k0*n*(ybar - m0)^2 / (2*k_n)
#
# k0 is a prior sample size in pseudo-observations: with k0 = 1e-6 against
# n = 48 real observations, m_n differs from ybar by ~4e-07.
#
# We draw sigma2 FIRST and then mu | sigma2 because that is the factorisation
# of the joint posterior: p(mu, sigma2 | y) = p(sigma2 | y) p(mu | sigma2, y).
# The marginal p(sigma2 | y) is Inverse-Gamma in closed form; the conditional
# p(mu | sigma2, y) is Normal. Drawing in that order gives exact joint draws.
#
# Iterations are INDEPENDENT of one another. Within a single iteration mu and
# sigma2 are DEPENDENT: a draw with large sigma2 comes paired with a more
# dispersed mu. Averaging over iterations therefore propagates variance
# uncertainty into the mean, which is the whole point.

#' Normal-Inverse-Gamma posterior draws for a stratum mean
#'
#' Exact joint draws from the conjugate posterior of a normal mean and
#' variance: first \eqn{\sigma^2} from its Inverse-Gamma marginal, then
#' \eqn{\mu \mid \sigma^2} from a normal. Draws are independent across
#' iterations; within a draw, the mean and variance are dependent. The default
#' prior is vague (`k0 = 1e-6` pseudo-observations for the mean).
#'
#' @param y Numeric outcomes for one stratum.
#' @param n_draws Number of posterior draws. This controls Monte Carlo error in
#'   posterior summaries only; it is not a sample size.
#' @param m0,k0 Prior mean, and prior sample size for the mean.
#' @param a0,b0 Inverse-Gamma shape and scale for the variance.
#' @param return_sigma2 If `TRUE`, also return the variance draws and the
#'   posterior hyperparameters.
#' @return A numeric vector of draws of the mean or, if `return_sigma2 = TRUE`,
#'   a list with `mu`, `sigma2` and `hyper`.
#' @examples
#' set.seed(1)
#' summary(nig_draws(rnorm(40, mean = 10, sd = 3), n_draws = 2000))
nig_draws <- function(y, n_draws = 20000, m0 = 0, k0 = 1e-6,
                      a0 = 1e-3, b0 = 1e-3, return_sigma2 = FALSE) {
  n <- length(y); ybar <- mean(y); ss <- sum((y - ybar)^2)

  kn <- k0 + n
  mn <- (k0 * m0 + n * ybar) / kn
  an <- a0 + n / 2
  bn <- b0 + 0.5 * ss + (k0 * n * (ybar - m0)^2) / (2 * kn)

  sigma2 <- bn / stats::rgamma(n_draws, shape = an, rate = 1)  # InvGamma(an, bn)
  mu     <- mn + stats::rnorm(n_draws) * sqrt(sigma2 / kn)

  if (return_sigma2) {
    return(list(mu = mu, sigma2 = sigma2,
                hyper = c(n = n, ybar = ybar, ss = ss,
                          kn = kn, mn = mn, an = an, bn = bn)))
  }
  mu
}

#' Conjugate posterior CGR curve
#'
#' Gives each of the four strata an independent Normal-Inverse-Gamma posterior
#' ([nig_draws()]) and pushes the draws through the CGRC reweighting at every
#' grid point. No MCMC is involved.
#'
#' @inheritParams estimand
#' @param n_draws Posterior draws per stratum. This controls Monte Carlo error
#'   in the summaries only; it is not a sample size.
#' @param prior Named list of prior hyperparameters for [nig_draws()] (`m0`,
#'   `k0`, `a0`, `b0`). An empty list uses the vague defaults.
#' @param direction `1` if higher outcome values are favourable, `-1` if lower
#'   values are.
#' @return A data frame with one row per grid point: `cgr`, `method`, the
#'   posterior mean `est`, `sd`, 95% credible limits `lo` and `hi`, `p_fav`
#'   (the posterior probability that the effect is favourable) and the Monte
#'   Carlo standard error `mcse`.
cgr_conjugate <- function(df, grid = seq(0, 1, length.out = 101),
                          n_draws = 20000, prior = list(),
                          direction = 1) {
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))

  d <- lapply(grid, cgr_delta, mu = mu, r = rat$r, s = rat$s)
  sdv <- vapply(d, stats::sd, numeric(1))

  out <- data.frame(
    cgr = grid, method = "conjugate",
    est = vapply(d, mean, numeric(1)), sd = sdv,
    lo = vapply(d, function(x) unname(stats::quantile(x, 0.025)), numeric(1)),
    hi = vapply(d, function(x) unname(stats::quantile(x, 0.975)), numeric(1)),
    # P(favourable): direction = +1 if higher is better, -1 if lower is better
    p_fav = vapply(d, function(x) mean(direction * x > 0), numeric(1)),
    stringsAsFactors = FALSE
  )
  out$mcse <- sdv / sqrt(n_draws)   # iid draws: no autocorrelation correction
  out
}
