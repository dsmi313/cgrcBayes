# Activated-expectancy-bias (AEB) generative model (Szigeti et al. 2023,
# Fig 1 and Eq 1-4) and repeated-simulation operating characteristics.
#
#   Eq 1  p(B_TRT = 1) = 0.5
#   Eq 2  p(B_PT = B_TRT) = p_CG
#   Eq 3  B_TE = B_PT
#   Eq 4  OUT = N_NH + B_TRT*N_DTE + B_TE*N_AEB
#
# Eq 4 does not say whether the DTE and AEB noise terms apply to the whole
# sample or only to the affected subgroup. noise = "all" (whole sample) is the
# default because it most closely reproduces the published results:
#   noise = "all": Hedges g = 0.4011 (published 0.40); Table 1 unadjusted
#                  significance 0.056/0.876/0.826/0.992 (published
#                  0.05/0.86/0.78/0.99)
#   noise = "arm": Hedges g = 0.5023; rates 0.056/0.970/0.930/0.998
# The variance scope has not been confirmed against the original simulation
# code.

#' Simulate a trial under the activated-expectancy-bias model
#'
#' `sim_aeb()` draws one trial from the activated-expectancy-bias (AEB) model
#' of Szigeti et al. (2023): random 1:1 allocation, a guess that is correct
#' with probability `p_cg`, and an outcome made of a baseline term, an optional
#' direct treatment effect (DTE) in the active arm, and an optional expectancy
#' effect for participants who believe they received the active treatment.
#'
#' `sim_aeb_unknown()` adds an "I don't know" response under two explicit
#' assumptions: each participant answers UNKNOWN with probability `u`,
#' independent of arm, and an UNKNOWN responder carries no expectancy. Among
#' directional responders the guess is correct with probability `p_cg`.
#'
#' @param n Number of participants.
#' @param p_cg Probability of a correct guess (among directional responders, in
#'   `sim_aeb_unknown()`).
#' @param u Probability of an UNKNOWN response.
#' @param dte_on Include the direct treatment effect.
#' @param aeb_on Include the expectancy effect.
#' @param noise `"all"` applies the DTE and AEB noise terms to every
#'   participant; `"arm"` applies them only within the affected subgroup.
#' @param mu_nh,sd_nh Mean and SD of the baseline term.
#' @param mu_dte,sd_dte Mean and SD of the direct treatment effect.
#' @param mu_aeb,sd_aeb Mean and SD of the expectancy effect.
#' @return A data frame with columns `condition`, `guess` and `value`, ready
#'   for [cgrc()] or [cgrc_unknown()].
#' @references Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D (2023). The
#'   difference between 'placebo group' and 'placebo control': a case study in
#'   psychedelic microdosing. *Scientific Reports* 13:12107.
#'   \doi{10.1038/s41598-023-34938-7}
#' @examples
#' set.seed(1)
#' head(sim_aeb(n = 100, p_cg = 0.7, dte_on = TRUE, aeb_on = TRUE))
#' table(sim_aeb_unknown(n = 100, u = 0.2)$guess)
sim_aeb <- function(n = 230, p_cg = 0.7, dte_on = FALSE, aeb_on = FALSE,
                    noise = c("all", "arm"),
                    mu_nh = 10, sd_nh = 4, mu_dte = 3, sd_dte = 6.2,
                    mu_aeb = 7.7, sd_aeb = 6.2) {
  noise <- match.arg(noise)
  trt     <- stats::runif(n) < 0.5
  correct <- stats::runif(n) < p_cg
  pt      <- ifelse(correct, trt, !trt)
  y <- stats::rnorm(n, mu_nh, sd_nh)
  if (dte_on) y <- y + if (noise == "arm") trt * stats::rnorm(n, mu_dte, sd_dte)
                       else trt * mu_dte + stats::rnorm(n, 0, sd_dte)
  if (aeb_on) y <- y + if (noise == "arm") pt * stats::rnorm(n, mu_aeb, sd_aeb)
                       else pt * mu_aeb + stats::rnorm(n, 0, sd_aeb)
  data.frame(condition = ifelse(trt, "AC", "PL"),
             guess = ifelse(pt, "AC", "PL"), value = y,
             stringsAsFactors = FALSE)
}

#' Operating characteristics of the CGR-adjusted estimator
#'
#' Simulates `n_trials` independent trials from [sim_aeb()] in each of the four
#' DTE by AEB scenarios and summarises how the adjusted estimate
#' \eqn{\Delta(0.5)} performs. This is the pre-flight check Szigeti et al.
#' (2023) recommend before applying CGR adjustment to a trial.
#'
#' Under this model the target \eqn{\Delta(0.5)} equals the direct treatment
#' effect: at a 50% correct-guess rate the guess distribution is the same in
#' both arms, so the expectancy term cancels. That identity is specific to the
#' simulation model, not a causal identification result for real trials.
#'
#' Trials with an empty stratum are skipped and counted in
#' `empty_stratum_rate`. A high rate means CGR adjustment is fragile at these
#' design parameters.
#'
#' `cgr_unknown_operating()` does the same for the six-stratum estimator
#' \eqn{\Delta(0.5, u_{obs})}, simulating with [sim_aeb_unknown()].
#'
#' @param n_trials Number of simulated trials per scenario.
#' @inheritParams sim_aeb
#' @param n_draws Posterior draws per stratum in each trial.
#' @param seed RNG seed.
#' @param mu_dte Direct treatment effect when DTE is on.
#' @param mu_aeb Expectancy effect when AEB is on.
#' @param direction `1` if higher outcome values are favourable, `-1` if lower
#'   values are.
#' @return A data frame with one row per scenario: the true effect; unadjusted
#'   and adjusted means and biases; adjusted RMSE; 95% credible-interval
#'   coverage; the rates of `P(favourable) > 0.95` for the adjusted and
#'   unadjusted Bayesian analyses (`p_fav_gt_95`, `unadj_p_fav_gt_95`) and of
#'   the stricter `P(favourable) > 0.975` (`p_fav_gt_975`); the two-sided
#'   t-test rejection rate `freq_sig`; `empty_stratum_rate`; and the number of
#'   valid trials.
#' @examples
#' cgr_operating(n_trials = 20, n = 120, p_cg = 0.75, n_draws = 500)
cgr_operating <- function(n_trials = 500, n = 230, p_cg = 0.7,
                          noise = "all", n_draws = 4000, seed = 1,
                          mu_dte = 3, mu_aeb = 7.7, direction = 1) {
  set.seed(seed)
  cfg <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  rows <- lapply(cfg, function(z) {
    dte <- as.logical(z[1]); aeb <- as.logical(z[2])
    truth <- if (dte) mu_dte else 0
    adj <- unadj <- numeric(n_trials)
    cov <- fav <- fav975 <- ufav <- sig <- valid <- logical(n_trials)
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb(n, p_cg, dte, aeb, noise, mu_dte = mu_dte, mu_aeb = mu_aeb)
      # At a high correct-guess rate with small n, a wrong-guess stratum can come
      # up empty, and the estimand is then undefined. Skip that trial rather than
      # crash - the RATE of such trials is itself a warning that CGR adjustment
      # is fragile at these parameters, and it is reported below.
      if (!all(STRATA %in% paste0(d$condition, d$guess))) next
      valid[i] <- TRUE
      st <- cgr_strata(d); rat <- cgr_ratios(st)
      mu <- lapply(st, nig_draws, n_draws = n_draws)
      dd <- cgr_delta(0.5, mu, rat$r, rat$s)
      q  <- stats::quantile(dd, c(0.025, 0.975))
      pd <- mean(direction * dd > 0)        # posterior prob. of a FAVOURABLE effect
      adj[i] <- mean(dd)
      cov[i] <- q[1] <= truth && truth <= q[2]
      fav[i]    <- pd > 0.95                # ADJUSTED Bayesian flag (standard)
      fav975[i] <- pd > 0.975               # stricter level, kept for research only
      a <- d$value[d$condition == "AC"]; b <- d$value[d$condition == "PL"]
      unadj[i] <- mean(a) - mean(b)
      # UNADJUSTED BAYESIAN posterior of the raw arm contrast (muAC - muPL) under
      # the same NIG model and the SAME P(favourable)>0.95 criterion as the
      # adjusted flag, so the trade-off compares like with like. Drawn on a
      # saved-and-restored RNG state so it leaves the other columns unchanged.
      .seed_save <- get(".Random.seed", envir = .GlobalEnv)
      udd <- nig_draws(a, n_draws = n_draws) - nig_draws(b, n_draws = n_draws)
      ufav[i] <- mean(direction * udd > 0) > 0.95
      assign(".Random.seed", .seed_save, envir = .GlobalEnv)
      # Conventional (two-sided) unadjusted t-test: a familiar frequentist
      # reference reported in the text only, NOT the plotted criterion.
      sig[i] <- stats::t.test(a, b, var.equal = TRUE)$p.value < 0.05
    }
    v <- which(valid); nv <- length(v)
    m <- function(x) if (nv) mean(x[v]) else NA_real_
    data.frame(DTE = z[1], AEB = z[2], true = truth,
               unadj_mean = m(unadj), unadj_bias = m(unadj) - truth,
               adj_mean = m(adj), adj_bias = m(adj) - truth,
               adj_rmse = if (nv) sqrt(mean((adj[v] - truth)^2)) else NA_real_,
               coverage95 = m(cov),
               p_fav_gt_95 = m(fav),            # ADJUSTED Bayesian P(fav)>0.95
               unadj_p_fav_gt_95 = m(ufav),     # UNADJUSTED Bayesian P(fav)>0.95
               p_fav_gt_975 = m(fav975),        # research only; not shown in the app
               freq_sig = m(sig),               # conventional two-sided t-test ref
               empty_stratum_rate = round(mean(!valid), 4), n_valid = nv)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

#' Design feasibility helpers
#'
#' Closed-form quantities for planning a trial. No simulation is needed.
#'
#' * `cgr_expected_strata()`: expected size of each of the four strata under
#'   1:1 allocation, with the same correct-guess rate in both arms.
#' * `cgr_min_stratum()`: expected size of the smallest stratum,
#'   `n * 0.5 * min(p_cg, 1 - p_cg)`. Below about 15, the wrong-guess strata
#'   are thin and simulated trials start to come up empty.
#' * `cgr_unknown_min_stratum()`: the same for the six strata of the UNKNOWN
#'   extension. Structurally empty classes (`u = 0` or `u = 1`) are ignored, so
#'   at `u = 0` it equals `cgr_min_stratum()`.
#' * `cgr_aeb_inflation()`: the expected inflation of an unadjusted estimate by
#'   expectancy under the AEB model, `mu_aeb * (2 * p_cg - 1)`.
#'
#' @param n Planned number of participants.
#' @param p_cg Expected correct-guess rate.
#' @param u Expected UNKNOWN-response rate.
#' @param mu_aeb Expectancy effect, in outcome units.
#' @return A number, or for `cgr_expected_strata()` a named vector of four.
#' @examples
#' cgr_expected_strata(120, 0.85)
#' cgr_min_stratum(120, 0.85)
#' cgr_aeb_inflation(7.7, 0.85)
#' @name design_helpers
cgr_min_stratum <- function(n, p_cg) n * 0.5 * pmin(p_cg, 1 - p_cg)

#' @rdname design_helpers
cgr_aeb_inflation <- function(mu_aeb, p_cg) mu_aeb * (2 * p_cg - 1)

#' @rdname sim_aeb
sim_aeb_unknown <- function(n = 230, p_cg = 0.7, u = 0.2,
                            dte_on = FALSE, aeb_on = FALSE, noise = c("all", "arm"),
                            mu_nh = 10, sd_nh = 4, mu_dte = 3, sd_dte = 6.2,
                            mu_aeb = 7.7, sd_aeb = 6.2) {
  noise <- match.arg(noise)
  trt     <- stats::runif(n) < 0.5
  is_unk  <- stats::runif(n) < u                 # A1: arm-independent UNKNOWN rate
  correct <- stats::runif(n) < p_cg
  pt      <- ifelse(correct, trt, !trt)          # perceived treatment (directional)
  b_te    <- ifelse(is_unk, 0, as.numeric(pt))   # A2: UNKNOWN -> no expectancy
  y <- stats::rnorm(n, mu_nh, sd_nh)
  if (dte_on) y <- y + if (noise == "arm") trt * stats::rnorm(n, mu_dte, sd_dte)
                       else trt * mu_dte + stats::rnorm(n, 0, sd_dte)
  if (aeb_on) y <- y + if (noise == "arm") b_te * stats::rnorm(n, mu_aeb, sd_aeb)
                       else b_te * mu_aeb + stats::rnorm(n, 0, sd_aeb)
  data.frame(condition = ifelse(trt, "AC", "PL"),
             guess = ifelse(is_unk, "UNKNOWN", ifelse(pt, "AC", "PL")),
             value = y, stringsAsFactors = FALSE)
}

#' @rdname cgr_operating
cgr_unknown_operating <- function(n_trials = 500, n = 230, p_cg = 0.7, u = 0.2,
                                  noise = "all", n_draws = 4000, seed = 1,
                                  mu_dte = 3, mu_aeb = 7.7, direction = 1) {
  set.seed(seed)
  cfg <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  rows <- lapply(cfg, function(z) {
    dte <- as.logical(z[1]); aeb <- as.logical(z[2])
    truth <- if (dte) mu_dte else 0
    adj <- unadj <- numeric(n_trials)
    cov <- fav <- fav975 <- ufav <- sig <- valid <- logical(n_trials)
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb_unknown(n, p_cg, u, dte, aeb, noise, mu_dte = mu_dte, mu_aeb = mu_aeb)
      st <- cgr_unknown_strata(d); o <- cgr_unknown_observed(st)
      # estimand defined at c = 0.5, u = u_obs: both directional classes present,
      # and the UNKNOWN class present if u_obs > 0. Otherwise skip and count it.
      if (o$n_correct == 0 || o$n_incorrect == 0 ||
          (o$u_obs > 0 && o$n_unknown == 0)) next
      valid[i] <- TRUE
      rat <- cgr_unknown_ratios(st)
      mu  <- stats::setNames(lapply(UNKNOWN_STRATA, function(nm) {
        y <- st[[nm]]; if (length(y)) nig_draws(y, n_draws = n_draws) else NA_real_
      }), UNKNOWN_STRATA)
      dd <- cgr_unknown_delta(0.5, o$u_obs, mu, rat$r, rat$s, rat$t)
      q  <- stats::quantile(dd, c(0.025, 0.975)); pd <- mean(direction * dd > 0)
      adj[i] <- mean(dd)
      cov[i] <- q[1] <= truth && truth <= q[2]
      fav[i] <- pd > 0.95; fav975[i] <- pd > 0.975   # 975: research only, not shown
      a <- d$value[d$condition == "AC"]; b <- d$value[d$condition == "PL"]
      unadj[i] <- mean(a) - mean(b)
      # Unadjusted BAYESIAN posterior of the raw arm contrast, same P(fav)>0.95
      # criterion as the adjusted flag (see cgr_operating); saved/restored RNG.
      .seed_save <- get(".Random.seed", envir = .GlobalEnv)
      udd <- nig_draws(a, n_draws = n_draws) - nig_draws(b, n_draws = n_draws)
      ufav[i] <- mean(direction * udd > 0) > 0.95
      assign(".Random.seed", .seed_save, envir = .GlobalEnv)
      # Conventional two-sided unadjusted t-test: familiar frequentist reference
      # only, reported in the text, not the plotted criterion.
      sig[i] <- stats::t.test(a, b, var.equal = TRUE)$p.value < 0.05
    }
    v <- which(valid); nv <- length(v)
    m <- function(x) if (nv) mean(x[v]) else NA_real_
    data.frame(DTE = z[1], AEB = z[2], true = truth, u = u,
               unadj_mean = m(unadj), unadj_bias = m(unadj) - truth,
               adj_mean = m(adj), adj_bias = m(adj) - truth,
               adj_rmse = if (nv) sqrt(mean((adj[v] - truth)^2)) else NA_real_,
               coverage95 = m(cov),
               p_fav_gt_95 = m(fav),            # ADJUSTED Bayesian P(fav)>0.95
               unadj_p_fav_gt_95 = m(ufav),     # UNADJUSTED Bayesian P(fav)>0.95
               p_fav_gt_975 = m(fav975),        # research only; not shown in the app
               freq_sig = m(sig),               # conventional two-sided t-test ref
               empty_stratum_rate = round(mean(!valid), 4), n_valid = nv)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

#' @rdname design_helpers
cgr_unknown_min_stratum <- function(n, p_cg, u) {
  # thinnest OCCUPIED cell: a structurally-empty class (u = 0 or u = 1) is not
  # counted, so at u = 0 this reduces to the binary cgr_min_stratum().
  cands <- c(0.5 * (1 - u) * p_cg, 0.5 * (1 - u) * (1 - p_cg), 0.5 * u)
  n * min(cands[cands > 0])
}
