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
#' `sim_aeb_unknown()` adds an "I don't know" response. By default it makes two
#' explicit assumptions: (A1) each participant answers UNKNOWN with probability
#' `u`, independent of arm, and (A2) an UNKNOWN responder carries no
#' expectancy. Among directional responders the guess is correct with
#' probability `p_cg`. `u_diff` and `aeb_unknown` relax A1 and A2, and
#' `u_wrong` and `u_outcome` make the UNKNOWN response depend on the unstated
#' belief or on the outcome, for stress tests (see [cgr_unknown_stress()]).
#'
#' @param n Number of participants.
#' @param p_cg Probability of a correct guess (among directional responders, in
#'   `sim_aeb_unknown()`).
#' @param u Probability of an UNKNOWN response, averaged over the two arms.
#' @param u_diff Difference in UNKNOWN rates, active minus placebo. The arms
#'   answer UNKNOWN with probabilities `u + u_diff / 2` and `u - u_diff / 2`.
#' @param aeb_unknown Fraction of the expectancy effect carried by UNKNOWN
#'   responders, applied through the belief they did not state.
#' @param u_wrong Extra probability of answering UNKNOWN for participants who
#'   would have guessed wrong, relative to those who would have guessed right.
#'   Applied so that each arm's average UNKNOWN rate is unchanged.
#' @param u_outcome Change in the probability of answering UNKNOWN per
#'   standard deviation of the outcome before expectancy. Positive values make
#'   participants who are doing better more likely to answer UNKNOWN.
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
                            mu_aeb = 7.7, sd_aeb = 6.2,
                            u_diff = 0, aeb_unknown = 0, u_wrong = 0,
                            u_outcome = 0) {
  noise <- match.arg(noise)
  u_ac <- u + u_diff / 2; u_pl <- u - u_diff / 2
  if (min(u_ac, u_pl) < 0 || max(u_ac, u_pl) > 1)
    stop("u +/- u_diff / 2 must lie in [0, 1].", call. = FALSE)
  trt     <- stats::runif(n) < 0.5
  r_unk   <- stats::runif(n)
  correct <- stats::runif(n) < p_cg
  pt      <- ifelse(correct, trt, !trt)          # perceived treatment (directional)
  y <- stats::rnorm(n, mu_nh, sd_nh)
  if (dte_on) y <- y + if (noise == "arm") trt * stats::rnorm(n, mu_dte, sd_dte)
                       else trt * mu_dte + stats::rnorm(n, 0, sd_dte)
  # Probability of answering UNKNOWN: the arm rate (A1 when u_diff = 0), shifted
  # for would-be-wrong guessers (u_wrong, mean-preserving within arm) and by the
  # standardised outcome before expectancy (u_outcome), clamped to [0, 1].
  p_unk <- ifelse(trt, u_ac, u_pl) + u_wrong * ((!correct) - (1 - p_cg))
  if (u_outcome != 0) p_unk <- p_unk + u_outcome * (y - mean(y)) / stats::sd(y)
  is_unk <- r_unk < pmin(pmax(p_unk, 0), 1)
  # A2 when aeb_unknown = 0: an UNKNOWN responder's latent belief pt carries a
  # fraction aeb_unknown of the expectancy effect.
  b_te    <- ifelse(is_unk, aeb_unknown * pt, as.numeric(pt))
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

#' Stress test of the UNKNOWN-preserving estimator
#'
#' Checks how the six-stratum estimator \eqn{\Delta(0.5, u_{obs})} behaves when
#' its two simulation assumptions fail, and compares it with the obvious
#' alternative: dropping UNKNOWN responders and running the binary CGRC on the
#' rest (complete case). Trials are simulated with [sim_aeb_unknown()] with
#' expectancy on, over a grid of
#'
#' * `u_diff`, the difference in UNKNOWN rates between arms (A1 holds at 0), and
#' * `aeb_unknown`, the fraction of the expectancy effect that UNKNOWN
#'   responders carry through the belief they did not state (A2 holds at 0),
#' * `u_wrong`, how much more often would-be-wrong guessers answer UNKNOWN, and
#' * `u_outcome`, how strongly answering UNKNOWN depends on the outcome before
#'   expectancy,
#'
#' with and without a direct treatment effect. With no direct effect, the rate
#' of `P(favourable) > 0.95` is a false-positive rate; with one, it is power.
#'
#' @param u_diff,aeb_unknown,u_wrong,u_outcome Values to sweep; see
#'   [sim_aeb_unknown()].
#' @inheritParams cgr_operating
#' @param u Average UNKNOWN-response rate.
#' @return A data frame with one row per combination of the swept values,
#'   scenario and
#'   method (`"six-stratum"` or `"complete-case"`): the true direct effect,
#'   bias, RMSE, 95% coverage, the rate of `P(favourable) > 0.95`, and the
#'   number of trials in which the estimand was defined.
#' @examples
#' cgr_unknown_stress(u_diff = c(0, 0.2), aeb_unknown = c(0, 1),
#'                    n_trials = 20, n_draws = 500)
cgr_unknown_stress <- function(u_diff = c(-0.2, -0.1, 0, 0.1, 0.2),
                               aeb_unknown = c(0, 0.25, 0.5, 1),
                               u_wrong = 0, u_outcome = 0,
                               n_trials = 500, n = 300, p_cg = 0.7, u = 0.2,
                               mu_dte = 3, mu_aeb = 7.7, n_draws = 2000,
                               seed = 1, direction = 1) {
  set.seed(seed)
  grid <- expand.grid(dte = c(0, 1), aeb_unknown = aeb_unknown, u_diff = u_diff,
                      u_wrong = u_wrong, u_outcome = u_outcome)
  rows <- lapply(seq_len(nrow(grid)), function(g) {
    dte <- grid$dte[g] == 1; truth <- if (dte) mu_dte else 0
    est <- matrix(NA_real_, n_trials, 2, dimnames = list(NULL, c("six", "cc")))
    cov <- fav <- est
    for (i in seq_len(n_trials)) {
      d  <- sim_aeb_unknown(n, p_cg, u, dte_on = dte, aeb_on = TRUE,
                            mu_dte = mu_dte, mu_aeb = mu_aeb,
                            u_diff = grid$u_diff[g],
                            aeb_unknown = grid$aeb_unknown[g],
                            u_wrong = grid$u_wrong[g],
                            u_outcome = grid$u_outcome[g])
      st <- cgr_unknown_strata(d); o <- cgr_unknown_observed(st)
      binary <- st[STRATA]
      if (any(lengths(binary) == 0) || o$n_unknown == 0 ||
          length(st$ACU) == 0 || length(st$PLU) == 0) next
      draws <- function(y) nig_draws(y, n_draws = n_draws)
      rat <- cgr_unknown_ratios(st)
      mu6 <- lapply(st[UNKNOWN_STRATA], draws)
      d6  <- cgr_unknown_delta(0.5, o$u_obs, mu6, rat$r, rat$s, rat$t)
      rb  <- cgr_ratios(binary)
      dcc <- cgr_delta(0.5, mu6[STRATA], rb$r, rb$s)   # same draws, UNKNOWN dropped
      for (m in c("six", "cc")) {
        dd <- if (m == "six") d6 else dcc
        q  <- stats::quantile(dd, c(0.025, 0.975))
        est[i, m] <- mean(dd)
        cov[i, m] <- q[1] <= truth && truth <= q[2]
        fav[i, m] <- mean(direction * dd > 0) > 0.95
      }
    }
    v <- !is.na(est[, 1])
    data.frame(u_diff = grid$u_diff[g], aeb_unknown = grid$aeb_unknown[g],
               u_wrong = grid$u_wrong[g], u_outcome = grid$u_outcome[g],
               DTE = grid$dte[g], true = truth,
               method = c("six-stratum", "complete-case"),
               bias = colMeans(est[v, , drop = FALSE]) - truth,
               rmse = sqrt(colMeans((est[v, , drop = FALSE] - truth)^2)),
               coverage95 = colMeans(cov[v, , drop = FALSE]),
               p_fav_gt_95 = colMeans(fav[v, , drop = FALSE]),
               n_valid = sum(v), row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Approximate bias of the UNKNOWN-preserving estimator
#'
#' A sensitivity calculation: how far \eqn{\Delta(0.5, u_{obs})} is expected to
#' move from the direct effect when its assumptions fail, under the
#' activated-expectancy-bias model,
#' \deqn{\mathrm{bias} \approx u \lambda \mu_{AEB} (2 p_{cg} - 1)
#'   - (1 - \lambda) \frac{u_{diff} \mu_{AEB}}{2},}
#' where \eqn{\lambda} is `aeb_unknown`. The first term is expectancy carried
#' by UNKNOWN responders, whose unstated beliefs still track their arm; the
#' second is the expectancy missing from whichever arm answers UNKNOWN more
#' often. Both vanish when A1 and A2 hold.
#'
#' This is a heuristic approximation, not an exact result. It agrees with
#' [cgr_unknown_stress()] to within simulation error at n = 300,
#' `p_cg = 0.7`, `u = 0.2` and at n = 400, `p_cg = 0.8`, `u = 0.35`; see
#' `docs/unknown-stress-test.md`. Use it to ask how large a violation would
#' have to be to change a conclusion, with `u_diff` set to the observed
#' difference in UNKNOWN rates and a range of plausible `aeb_unknown`.
#'
#' @param u UNKNOWN-response rate, averaged over arms.
#' @param u_diff Difference in UNKNOWN rates, active minus placebo.
#' @param aeb_unknown Fraction of the expectancy effect carried by UNKNOWN
#'   responders, between 0 and 1.
#' @param mu_aeb Expectancy effect, in outcome units.
#' @param p_cg Correct-guess rate among directional responders.
#' @return The approximate bias, in outcome units (vectorised).
#' @examples
#' # Bias if UNKNOWN responders keep none, half or all of a 7.7-point
#' # expectancy effect, at a 20% UNKNOWN rate and 70% directional guessing
#' cgr_unknown_bias_approx(u = 0.2, aeb_unknown = c(0, 0.5, 1), mu_aeb = 7.7,
#'                         p_cg = 0.7)
cgr_unknown_bias_approx <- function(u, u_diff = 0, aeb_unknown = 0, mu_aeb,
                                    p_cg) {
  u * aeb_unknown * mu_aeb * (2 * p_cg - 1) -
    (1 - aeb_unknown) * u_diff * mu_aeb / 2
}

#' @rdname design_helpers
cgr_unknown_min_stratum <- function(n, p_cg, u) {
  # thinnest OCCUPIED cell: a structurally-empty class (u = 0 or u = 1) is not
  # counted, so at u = 0 this reduces to the binary cgr_min_stratum().
  cands <- c(0.5 * (1 - u) * p_cg, 0.5 * (1 - u) * (1 - p_cg), 0.5 * u)
  n * min(cands[cands > 0])
}
