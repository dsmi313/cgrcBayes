#' Region of practical equivalence (ROPE) decomposition
#'
#' Splits the posterior of the adjusted contrast at each CGR into three
#' exhaustive regions relative to the band `[-delta, delta]`: meaningful harm,
#' practically negligible, and meaningful benefit. The three probabilities sum
#' to one. Unlike a single tail probability, this separates a posterior
#' concentrated near zero (a negligible effect) from one too diffuse to say
#' anything. `delta` must be declared, which puts the question of practical
#' significance in the open; the default is `0.1 * sd(value)`.
#'
#' `cgr_rope_sensitivity()` shows how the conclusion at one CGR moves as
#' `delta` widens. The `cgr_unknown_` versions apply the same decomposition to
#' the six-stratum UNKNOWN extension.
#'
#' A Bayes factor against a point null is deliberately not offered: under the
#' vague default prior it depends on the prior width without limit (the
#' Jeffreys-Lindley paradox).
#'
#' @param df A data frame with columns `condition`, `guess` and `value`, coded
#'   as for [cgr_strata()] (binary) or [cgr_unknown_strata()] (UNKNOWN).
#' @inheritParams cgr_conjugate
#' @param u_target Target UNKNOWN rate. `NULL` uses the observed rate.
#' @param delta Half-width of the band, in outcome units. `NULL` uses
#'   `delta_sd_frac * sd(value)`.
#' @param delta_sd_frac Half-width as a fraction of the outcome SD.
#' @param at_cgr CGR at which to evaluate the sensitivity.
#' @param fracs Half-widths to try, as fractions of the outcome SD.
#' @return `cgr_rope()`: a data frame with one row per grid point, giving
#'   `delta`, the posterior mean, 95% and 50% intervals, and `p_harm`,
#'   `p_negligible` and `p_benefit`. `cgr_rope_sensitivity()`: a data frame
#'   with one row per fraction.
cgr_rope <- function(df, grid = seq(0, 1, length.out = 101),
                     n_draws = 20000, delta = NULL, delta_sd_frac = 0.1,
                     direction = 1, prior = list()) {
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(df$value)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))
  d   <- lapply(grid, cgr_delta, mu = mu, r = rat$r, s = rat$s)
  q   <- function(x, p) unname(stats::quantile(x, p))

  data.frame(
    cgr = grid, delta = delta,
    est   = vapply(d, mean, numeric(1)),
    lo95  = vapply(d, function(x) q(x, 0.025), numeric(1)),
    hi95  = vapply(d, function(x) q(x, 0.975), numeric(1)),
    lo50  = vapply(d, function(x) q(x, 0.250), numeric(1)),
    hi50  = vapply(d, function(x) q(x, 0.750), numeric(1)),
    # region probabilities are oriented by `direction` so + always = benefit
    p_harm      = vapply(d, function(x) mean(direction * x < -delta), numeric(1)),
    p_negligible= vapply(d, function(x) mean(abs(x) <= delta),        numeric(1)),
    p_benefit   = vapply(d, function(x) mean(direction * x >  delta), numeric(1)),
    stringsAsFactors = FALSE)
}

#' @rdname cgr_rope
cgr_rope_sensitivity <- function(df, at_cgr = 0.5,
                                 fracs = c(0.05, 0.10, 0.20, 0.30),
                                 n_draws = 12000, direction = 1) {
  sdy <- stats::sd(df$value)
  do.call(rbind, lapply(fracs, function(f) {
    z <- cgr_rope(df, grid = at_cgr, n_draws = n_draws,
                  delta_sd_frac = f, direction = direction)
    data.frame(delta_in_SD = f, delta = f * sdy,
               p_negligible = z$p_negligible, p_benefit = z$p_benefit,
               stringsAsFactors = FALSE)
  }))
}
