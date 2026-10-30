# Reimplementation of the original Szigeti et al. (2023) resampling procedure,
# for comparison. The paper fits a KDE per stratum with scikit-learn's
# KernelDensity, "all parameters left at default value": a Gaussian kernel with
# a fixed bandwidth of 1.0, which is narrow on a 0-100 VAS. Because the
# estimand uses only stratum means and KDE smoothing preserves means, the
# bandwidth does not affect the point estimate. It would matter for a quantile
# or tail-based estimand.

rkde <- function(y, m, bw = 1.0, round_scores = FALSE) {
  # sample from a Gaussian KDE: pick a datum, add N(0, bw^2)
  out <- y[sample.int(length(y), m, replace = TRUE)] + stats::rnorm(m, 0, bw)
  # Python's round() and R's round() both use ties-to-even on ordinary numeric
  # inputs. Ties have probability zero after Gaussian KDE sampling.
  if (round_scores) out <- round(out)
  out
}

#' Original KDE resampling procedure
#'
#' Reimplements the resampling procedure of Szigeti et al. (2023) for
#' comparison with the analytic and Bayesian versions. Each stratum gets a
#' Gaussian kernel density estimate with fixed bandwidth `bw` (scikit-learn's
#' default of 1.0). Pseudo-samples are drawn so the combined data have the
#' target CGR, the arms are compared with a pooled-variance t-test, and
#' estimates and p-values are averaged over `n_rep` resamples.
#'
#' With `source_faithful = TRUE`, stratum proportions are rounded to two
#' decimals and pseudo-scores to integers, as in the original Python code.
#' With `FALSE`, exact ratios and unrounded scores are used.
#'
#' * `cgr_kde()` runs the procedure at one CGR.
#' * `cgr_kde_curve()` runs it over a grid and returns the mean and
#'   across-resample SD of the estimate and p-value, for [szigeti_panel()]. The
#'   default 13-point grid matches the original code's configuration.
#' * `cgr_kde_ladder()` compares the KDE average at several resample counts
#'   with the analytic value, separating Monte Carlo noise from any real
#'   difference.
#'
#' Averaging p-values across resamples has no accepted inferential
#' interpretation. It is reproduced here to characterise the original method.
#'
#' @inheritParams estimand
#' @param cgr Target correct-guess rate.
#' @param n_rep Number of resamples.
#' @param reps Resample counts to compare.
#' @param bw Kernel bandwidth.
#' @param source_faithful Reproduce the original code's rounding.
#' @return `cgr_kde()`: a list with the mean estimate `est`, its Monte Carlo
#'   standard error `est_mcse`, the mean p-value `p`, the per-resample `draws`,
#'   the `target_sizes`, and flags. `cgr_kde_curve()`: a data frame with `cgr`,
#'   `est`, `est_sd`, `p` and `p_sd`. `cgr_kde_ladder()`: a data frame with one
#'   row per resample count.
#' @references Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D (2023). The
#'   difference between 'placebo group' and 'placebo control': a case study in
#'   psychedelic microdosing. *Scientific Reports* 13:12107.
#'   \doi{10.1038/s41598-023-34938-7}
cgr_kde <- function(df, cgr, n_rep = 100, bw = 1.0,
                    source_faithful = TRUE) {
  st <- cgr_strata(df)
  # The shipped Python first rounds stratum proportions to two decimals and
  # rounds every KDE pseudo-observation to an integer. These source artefacts
  # are deliberately isolated from the exact analytic/Bayesian estimand.
  nk <- cgr_sizes(st, cgr, legacy_round = source_faithful)

  est <- numeric(n_rep); pval <- numeric(n_rep); all_integer <- TRUE
  for (i in seq_len(n_rep)) {
    ys <- list(); trt <- list()
    for (k in STRATA) {
      m <- nk[[k]]
      if (m <= 0) next
      ys[[k]]  <- rkde(st[[k]], m, bw, round_scores = source_faithful)
      all_integer <- all_integer && all(ys[[k]] == round(ys[[k]]))
      trt[[k]] <- rep(substr(k, 1, 2) == "AC", m)
    }
    y <- unlist(ys, use.names = FALSE)
    t <- unlist(trt, use.names = FALSE)
    tt <- stats::t.test(y[t], y[!t], var.equal = TRUE)
    est[i]  <- unname(diff(rev(tt$estimate)))
    pval[i] <- tt$p.value
  }
  list(est = mean(est), est_mcse = stats::sd(est) / sqrt(n_rep),
       p = mean(pval), n_rep = n_rep, draws = est,
       target_sizes = nk, source_faithful = source_faithful,
       scores_integer = all_integer)
}

#' @rdname cgr_kde
cgr_kde_curve <- function(df, grid = seq(0, 1, length.out = 13),
                          n_rep = 80, bw = 1.0,
                          source_faithful = TRUE) {
  st <- cgr_strata(df)
  do.call(rbind, lapply(grid, function(cc) {
    nk <- cgr_sizes(st, cc, legacy_round = source_faithful)
    e <- p <- numeric(n_rep)
    for (i in seq_len(n_rep)) {
      ys <- list(); tr <- list()
      for (k in STRATA) {
        m <- nk[[k]]; if (m <= 0) next
        ys[[k]]  <- rkde(st[[k]], m, bw, round_scores = source_faithful)
        tr[[k]] <- rep(substr(k, 1, 2) == "AC", m)
      }
      y <- unlist(ys, use.names = FALSE); t <- unlist(tr, use.names = FALSE)
      if (sum(t) < 2 || sum(!t) < 2) { e[i] <- NA; p[i] <- NA; next }
      e[i] <- mean(y[t]) - mean(y[!t])
      p[i] <- stats::t.test(y[t], y[!t], var.equal = TRUE)$p.value
    }
    data.frame(cgr = cc, est = mean(e, na.rm = TRUE),
               est_sd = stats::sd(e, na.rm = TRUE),
               p = mean(p, na.rm = TRUE), p_sd = stats::sd(p, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

#' @rdname cgr_kde
cgr_kde_ladder <- function(df, cgr, reps = c(100, 1000, 10000), bw = 1.0,
                           source_faithful = FALSE) {
  st <- cgr_strata(df); rat <- cgr_ratios(st)
  analytic <- cgr_delta(cgr, lapply(st, mean), rat$r, rat$s)
  out <- do.call(rbind, lapply(reps, function(m) {
    k <- cgr_kde(df, cgr, n_rep = m, bw = bw,
                 source_faithful = source_faithful)
    data.frame(n_rep = m, kde_est = k$est, kde_mcse = k$est_mcse,
               mean_p = k$p, analytic = analytic,
               diff_from_analytic = k$est - analytic)
  }))
  rownames(out) <- NULL
  out
}
