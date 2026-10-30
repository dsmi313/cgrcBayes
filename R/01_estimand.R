# The four-stratum CGRC estimand. Base R only, so a masked generic such as
# dplyr::filter() over stats::filter() cannot change behaviour.

#' Four-stratum CGRC estimand
#'
#' Building blocks of the Correct Guess Rate Curve estimand of Szigeti et al.
#' (2023). Crossing assigned arm with guessed arm gives four strata, written
#' arm first: `ACAC`, `ACPL`, `PLAC`, `PLPL` (the `STRATA` constant).
#' Reweighting to a target correct-guess rate `c` preserves the observed arm
#' split within the correct-guess class (`r`, placebo's share) and within the
#' incorrect-guess class (`s`, active's share), then renormalises each arm by
#' its own total weight.
#'
#' * `cgr_strata()` splits outcomes into the four strata. It errors on an empty
#'   stratum, where the estimand is undefined.
#' * `cgr_observed()` returns the observed correct-guess rate.
#' * `cgr_ratios()` returns the within-class arm shares `r` and `s`.
#' * `cgr_weights()` returns the target weights `w_ACAC = c(1 - r)`,
#'   `w_ACPL = (1 - c)s`, `w_PLAC = (1 - c)(1 - s)` and `w_PLPL = cr`.
#'   Correct-guess weights sum to `c`; incorrect-guess weights sum to `1 - c`.
#' * `cgr_delta()` returns the reweighted active-minus-placebo contrast
#'   \eqn{\Delta(c)}.
#' * `cgr_sizes()` returns integer stratum sizes at a target CGR, as used by
#'   the original resampling procedure. It reproduces the supplementary worked
#'   example of Szigeti et al. (65/35/25/75 becomes 54/58/42/46 at CGR 0.50).
#' * `cgr_analytic()` evaluates \eqn{\Delta(c)} from the stratum sample means
#'   over a grid.
#'
#' At the observed CGR the reweighting is a no-op, so \eqn{\Delta(c_{obs})}
#' equals the raw difference in arm means (the observed-CGR identity).
#'
#' @param df A data frame with one row per participant and columns
#'   `condition` (`"AC"` or `"PL"`), `guess` (`"AC"` or `"PL"`) and numeric
#'   `value`.
#' @param st A list of stratum outcome vectors from `cgr_strata()`.
#' @param legacy_round If `TRUE`, round stratum proportions to two decimals
#'   first, as the original Python implementation does.
#' @param c,cgr Target correct-guess rate, between 0 and 1.
#' @param r,s Within-class arm shares from `cgr_ratios()`.
#' @param mu Named list of stratum means: four numbers, or four equal-length
#'   vectors of posterior draws.
#' @param grid Correct-guess rates at which to evaluate the curve.
#' @return `cgr_strata()`: a named list of four numeric vectors.
#'   `cgr_observed()`: a number. `cgr_ratios()`: a list with `r` and `s`.
#'   `cgr_weights()` and `cgr_sizes()`: named numeric vectors. `cgr_delta()`: a
#'   number, or a vector of draws when `mu` holds draws. `cgr_analytic()`: a
#'   data frame with columns `cgr`, `method` and `est`.
#' @references Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D (2023). The
#'   difference between 'placebo group' and 'placebo control': a case study in
#'   psychedelic microdosing. *Scientific Reports* 13:12107.
#'   \doi{10.1038/s41598-023-34938-7}
#' @examples
#' set.seed(1)
#' d   <- sim_aeb(n = 200, p_cg = 0.7, dte_on = TRUE)
#' st  <- cgr_strata(d)
#' rat <- cgr_ratios(st)
#' cgr_observed(st)
#' cgr_weights(0.5, rat$r, rat$s)
#' cgr_delta(0.5, lapply(st, mean), rat$r, rat$s)
#' @format NULL
#' @name estimand
STRATA <- c("ACAC", "ACPL", "PLAC", "PLPL")

#' @rdname estimand
cgr_strata <- function(df) {
  key <- factor(paste0(df$condition, df$guess), levels = STRATA)
  if (anyNA(key)) stop("condition/guess must be coded AC/PL", call. = FALSE)
  out <- split(df$value, key)
  empty <- STRATA[lengths(out) == 0L]
  if (length(empty)) {
    stop("empty stratum: ", paste(empty, collapse = ", "),
         "; the estimand is undefined", call. = FALSE)
  }
  out
}

#' @rdname estimand
cgr_observed <- function(st) {
  (length(st$ACAC) + length(st$PLPL)) / sum(lengths(st))
}

#' @rdname estimand
cgr_ratios <- function(st, legacy_round = FALSE) {
  n <- sum(lengths(st)); rho <- lengths(st) / n
  if (legacy_round) rho <- round(rho, 2)   # reproduces Szigeti's shipped code
  list(r = unname(rho[["PLPL"]] / (rho[["PLPL"]] + rho[["ACAC"]])),
       s = unname(rho[["ACPL"]] / (rho[["ACPL"]] + rho[["PLAC"]])))
}

#' @rdname estimand
cgr_weights <- function(c, r, s) {
  c(ACAC = c * (1 - r), ACPL = (1 - c) * s,
    PLAC = (1 - c) * (1 - s), PLPL = c * r)
}

#' @rdname estimand
cgr_delta <- function(c, mu, r, s) {
  w <- cgr_weights(c, r, s)
  den_ac <- w[["ACAC"]] + w[["ACPL"]]
  den_pl <- w[["PLPL"]] + w[["PLAC"]]
  if (den_ac <= 0 || den_pl <= 0) {
    stop("degenerate weights at cgr = ", c, call. = FALSE)
  }
  (w[["ACAC"]] * mu$ACAC + w[["ACPL"]] * mu$ACPL) / den_ac -
    (w[["PLPL"]] * mu$PLPL + w[["PLAC"]] * mu$PLAC) / den_pl
}

#' @rdname estimand
cgr_sizes <- function(st, cgr, legacy_round = FALSE) {
  n <- sum(lengths(st)); rat <- cgr_ratios(st, legacy_round)
  n_corr <- round(cgr * n); n_inc <- n - n_corr
  plpl <- round(n_corr * rat$r); acpl <- round(n_inc * rat$s)
  c(ACAC = n_corr - plpl, ACPL = acpl, PLAC = n_inc - acpl, PLPL = plpl)
}

#' @rdname estimand
cgr_analytic <- function(df, grid = seq(0, 1, length.out = 101)) {
  st <- cgr_strata(df); rat <- cgr_ratios(st); mu <- lapply(st, mean)
  data.frame(cgr = grid, method = "analytic",
             est = vapply(grid, cgr_delta, numeric(1),
                          mu = mu, r = rat$r, s = rat$s),
             stringsAsFactors = FALSE)
}

#' Posterior draws of the within-class arm shares
#'
#' An optional sensitivity analysis. The CGRC conditions on the observed
#' within-class arm shares; these functions instead draw them from Beta
#' posteriors, propagating sampling error in the stratum composition. That is
#' a different, larger uncertainty statement than the original estimand, so it
#' is off by default everywhere.
#'
#' In `cgr_unknown_ratio_draws()`, a share whose class has one empty arm cell
#' is held at its plug-in value (0 or 1), and the share of an empty class is
#' `NA`.
#'
#' @param st Stratum list from [cgr_strata()] or [cgr_unknown_strata()].
#' @param n_draws Number of draws.
#' @param alpha Beta prior pseudo-count added to each cell.
#' @return A list of draws of `r` and `s`, and of `t` for the UNKNOWN version.
cgr_ratio_draws <- function(st, n_draws, alpha = 1) {
  n <- lengths(st)
  list(r = stats::rbeta(n_draws, n[["PLPL"]] + alpha, n[["ACAC"]] + alpha),
       s = stats::rbeta(n_draws, n[["ACPL"]] + alpha, n[["PLAC"]] + alpha))
}

#' Reference-line test
#'
#' Checks a correct-guess-rate reference line against the observed-CGR
#' identity. At the observed CGR the reweighting is a no-op, so
#' \eqn{\Delta(c_{obs})} equals the raw active-minus-placebo mean difference:
#' the unadjusted estimate a paper reports. A reference line drawn at some
#' claimed CGR is only in the right place if \eqn{\Delta} there equals that
#' value. A large `err_at_orig` with a near-zero `err_at_obs` is the signature
#' of a misplaced line.
#'
#' `cgr_unknown_reference_line_test()` is the six-stratum version, evaluated at
#' the observed directional CGR with the UNKNOWN rate held at its observed
#' value.
#'
#' @param df A data frame with columns `condition`, `guess` and `value`, coded
#'   as for [cgr_strata()] (binary) or [cgr_unknown_strata()] (UNKNOWN).
#' @param orig_cgr The claimed reference CGR. In the UNKNOWN version, `NULL`
#'   uses the observed directional CGR.
#' @param published_unadj Optional published unadjusted estimate to compare
#'   against.
#' @return A one-row data frame: the computed observed CGR, \eqn{\Delta} at
#'   the observed and claimed CGRs, the raw mean difference, and the errors
#'   against `published_unadj`.
cgr_reference_line_test <- function(df, orig_cgr, published_unadj = NA_real_) {
  st  <- cgr_strata(df); rat <- cgr_ratios(st); mu <- lapply(st, mean)
  o   <- cgr_observed(st)
  raw <- mean(df$value[df$condition == "AC"]) -
         mean(df$value[df$condition == "PL"])
  d_obs  <- cgr_delta(o,        mu, rat$r, rat$s)
  d_orig <- cgr_delta(orig_cgr, mu, rat$r, rat$s)
  data.frame(
    computed_obs_cgr = o,
    D_at_obs         = d_obs,
    raw_mean_diff    = raw,
    published_unadj  = published_unadj,
    orig_cgr         = orig_cgr,
    D_at_orig_cgr    = d_orig,
    err_at_obs       = d_obs  - published_unadj,
    err_at_orig      = d_orig - published_unadj,
    stringsAsFactors = FALSE)
}
