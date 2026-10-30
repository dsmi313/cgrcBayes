# One-call entry points: cgrc() adjusts a trial and cgrc_headline() summarises
# it. cgr_operating() in R/05_sim.R checks a design before adjustment.

#' CGR-adjusted analysis of a trial
#'
#' The main entry point. Fits the conjugate posterior curve
#' ([cgr_conjugate()]) on a grid that includes the exact observed CGR, then
#' summarises it at the observed CGR (unadjusted) and at a target CGR of 0.50
#' (guessing at chance). The 0.50 target is a sensitivity analysis, not proof
#' of perfect blinding.
#'
#' Check that adjustment is reliable for your design first with
#' [cgr_operating()].
#'
#' @inheritParams cgr_conjugate
#' @param seed Optional RNG seed. When given, it is set before drawing and
#'   stored in the result.
#' @param outcome Optional outcome label for the summary table.
#' @param x A `cgrc` object.
#' @param ... For `plot()`, passed to [cgr_plot()]. Ignored by `print()`.
#' @return An object of class `cgrc`: a list with the posterior `curve`, the
#'   `summary` table from [cgr_summary_table()], `observed_cgr` and `seed`.
#' @examples
#' set.seed(1)
#' trial <- sim_aeb(n = 230, p_cg = 0.7, dte_on = TRUE, aeb_on = TRUE)
#' fit <- cgrc(trial, n_draws = 4000)
#' fit
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
cgrc <- function(df, n_draws = 20000, direction = 1, prior = list(), seed = NULL,
                 outcome = "") {
  if (!is.null(seed)) set.seed(seed)
  o    <- cgr_observed(cgr_strata(df))
  grid <- sort(unique(c(seq(0, 1, length.out = 101), o)))
  cur  <- cgr_conjugate(df, grid, n_draws = n_draws,
                        direction = direction, prior = prior)
  structure(list(curve = cur,
                 summary = cgr_summary_table(cur, o, label = outcome),
                 observed_cgr = o, seed = seed),
            class = "cgrc")
}

#' @rdname cgrc
print.cgrc <- function(x, ...) {
  cat(sprintf("CGRC-adjusted analysis  (observed CGR = %.4f)\n", x$observed_cgr))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @rdname cgrc
plot.cgrc <- function(x, ...) {
  cgr_plot(x$curve, obs_cgr = x$observed_cgr, ...)
}

#' Two-probability headline for a CGR-adjusted analysis
#'
#' Summarises a trial with two plain probabilities, each before and after
#' reweighting to a CGR of 0.50: the probability that the effect is favourable,
#' and the probability that it exceeds a meaningful difference `delta`. No
#' significance cut-off is applied. The default `delta` is half the outcome SD,
#' the minimal important difference suggested by Norman et al. (2003).
#'
#' @inheritParams cgrc
#' @param delta Meaningful-difference threshold, in outcome units. `NULL` uses
#'   `delta_sd_frac * sd(value)`.
#' @param delta_sd_frac Threshold as a fraction of the outcome SD.
#' @param x A `cgrc_headline` object.
#' @param ... Ignored.
#' @return An object of class `cgrc_headline`: a list with the observed CGR,
#'   `delta`, the adjusted estimate and its 95% credible interval
#'   (`adj_est`, `adj_lo`, `adj_hi`), the four probabilities (`p_dir_obs`,
#'   `p_dir_blind`, `p_meaningful_obs`, `p_meaningful_blind`) and a
#'   one-paragraph `text` summary.
#' @examples
#' set.seed(1)
#' cgrc_headline(sim_aeb(n = 230, p_cg = 0.7, dte_on = TRUE), n_draws = 4000)
cgrc_headline <- function(df, direction = 1, delta = NULL, delta_sd_frac = 0.5,
                          n_draws = 20000, prior = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  st  <- cgr_strata(df); rat <- cgr_ratios(st)
  if (is.null(delta)) delta <- delta_sd_frac * stats::sd(df$value)
  mu  <- lapply(st, function(y)
    do.call(nig_draws, c(list(y = y, n_draws = n_draws), prior)))
  o   <- cgr_observed(st)
  d_obs <- cgr_delta(o,   mu, rat$r, rat$s)
  d_bl  <- cgr_delta(0.5, mu, rat$r, rat$s)
  eff_obs <- direction * d_obs; eff_bl <- direction * d_bl   # + = favourable
  q <- function(x, p) unname(stats::quantile(x, p))
  res <- list(
    observed_cgr = o, delta = delta, direction = direction, seed = seed,
    adj_est = mean(d_bl), adj_lo = q(d_bl, 0.025), adj_hi = q(d_bl, 0.975),
    p_dir_obs        = mean(eff_obs > 0),     p_dir_blind        = mean(eff_bl > 0),
    p_meaningful_obs = mean(eff_obs > delta), p_meaningful_blind = mean(eff_bl > delta))
  dtxt <- trimws(sprintf("%.2g", delta))
  unit <- if (identical(dtxt, "1")) "point" else "points"
  res$text <- sprintf(paste0(
    "Raw, this trial shows a %.0f%% probability of a favourable effect and %.0f%% ",
    "that it is meaningful (beyond %s %s). Reweighted to a correct-guess rate of ",
    "0.50 (guessing at chance), those become %.0f%% and %.0f%% (adjusted effect ",
    "%.2f, 95%% CrI %.2f to %.2f)."),
    100 * res$p_dir_obs, 100 * res$p_meaningful_obs, dtxt, unit,
    100 * res$p_dir_blind, 100 * res$p_meaningful_blind,
    res$adj_est, res$adj_lo, res$adj_hi)
  structure(res, class = "cgrc_headline")
}

#' @rdname cgrc_headline
print.cgrc_headline <- function(x, ...) { cat(x$text, "\n"); invisible(x) }
