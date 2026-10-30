# JAGS backend. The prior is the Normal-Inverse-Gamma model of the conjugate
# backend: under NIG, Var(mu | sigma2) = sigma2 / k0, so the precision of mu is
# k0 * tau and the prior is written mu[j] ~ dnorm(m0, k0 * tau[j]). The
# hyperparameters are passed as data so the two backends cannot drift apart;
# cgr_check_backends() confirms they agree.

#' JAGS backend for the CGR curve
#'
#' Fits the four-stratum model in JAGS and evaluates \eqn{\Delta(c)} inside
#' the model. With `likelihood = "normal"`, the prior is the same
#' Normal-Inverse-Gamma model as [cgr_conjugate()]; [cgr_check_backends()]
#' verifies that the two agree. `likelihood = "t"` is a robust Student-t
#' sensitivity analysis with estimated degrees of freedom
#' (`nu ~ dexp(0.1)` truncated to 2 to 100).
#'
#' `jags_model_string()` returns the model code.
#'
#' Requires the \pkg{rjags} package and a system JAGS installation.
#'
#' @inheritParams cgr_conjugate
#' @param likelihood `"normal"` or `"t"`.
#' @param n_iter Sampling iterations per chain.
#' @param n_burn Burn-in iterations per chain.
#' @param n_chains Number of chains.
#' @param seed Base RNG seed. Chain `i` uses `seed + i`.
#' @param prior Named list with `m0`, `k0`, `a0` and `b0`.
#' @return `cgr_jags()`: a data frame like that of [cgr_conjugate()], with
#'   added `ess` and `rhat` columns and an autocorrelation-corrected `mcse`. For
#'   the t likelihood, the posterior mean of `nu` is attached as attribute
#'   `"nu"`. `jags_model_string()`: a character string.
jags_model_string <- function(likelihood = c("normal", "t")) {
  likelihood <- match.arg(likelihood)
  lik <- if (likelihood == "normal") {
    "  for (i in 1:N) { y[i] ~ dnorm(mu[k[i]], tau[k[i]]) }\n"
  } else {
    paste0("  for (i in 1:N) { y[i] ~ dt(mu[k[i]], tau[k[i]], nu) }\n",
           "  nu ~ dexp(0.1) T(2, 100)\n")
  }
  paste0(
    "model {\n", lik,
    "  for (j in 1:4) {\n",
    "    tau[j] ~ dgamma(a0, b0)\n",
    "    mu[j]  ~ dnorm(m0, k0 * tau[j])\n",   # <- conditional precision
    "  }\n",
    "  for (m in 1:M) {\n",
    "    w_acac[m] <- cgr[m] * (1 - r)\n",
    "    w_acpl[m] <- (1 - cgr[m]) * s\n",
    "    w_plpl[m] <- cgr[m] * r\n",
    "    w_plac[m] <- (1 - cgr[m]) * (1 - s)\n",
    "    delta[m] <- (w_acac[m]*mu[1] + w_acpl[m]*mu[2]) / (w_acac[m]+w_acpl[m]) -\n",
    "                (w_plpl[m]*mu[4] + w_plac[m]*mu[3]) / (w_plpl[m]+w_plac[m])\n",
    "  }\n}\n")
}

#' @rdname jags_model_string
cgr_jags <- function(df, grid = seq(0, 1, length.out = 101),
                     likelihood = c("normal", "t"),
                     n_iter = 10000, n_burn = 2000, n_chains = 4, seed = 1,
                     prior = list(m0 = 0, k0 = 1e-6, a0 = 1e-3, b0 = 1e-3),
                     direction = 1) {

  likelihood <- match.arg(likelihood)
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("needs rjags + a JAGS install; see https://mcmc-jags.sourceforge.io",
         call. = FALSE)
  }

  st <- cgr_strata(df); rat <- cgr_ratios(st)
  eps <- 1e-6
  grid_safe <- pmin(pmax(grid, eps), 1 - eps)   # c in {0,1} can zero a weight

  y <- unlist(st, use.names = FALSE)
  k <- rep(seq_along(st), lengths(st))

  dat <- c(list(y = y, k = k, N = length(y),
                cgr = grid_safe, M = length(grid_safe),
                r = rat$r, s = rat$s), prior)

  inits <- lapply(seq_len(n_chains), function(i)
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + i))

  con <- textConnection(jags_model_string(likelihood))
  on.exit(close(con), add = TRUE)
  m <- rjags::jags.model(con, data = dat, inits = inits,
                         n.chains = n_chains, quiet = TRUE)
  stats::update(m, n_burn)
  mon <- if (likelihood == "t") c("delta", "nu") else "delta"
  samp <- rjags::coda.samples(m, variable.names = mon, n.iter = n_iter)

  all <- as.matrix(samp)
  dc <- grep("^delta\\[", colnames(all), value = TRUE)
  dc <- dc[order(as.integer(sub("^delta\\[(\\d+)\\]$", "\\1", dc)))]
  dr <- all[, dc, drop = FALSE]

  ess <- rep(NA_real_, length(dc)); rh <- ess
  if (requireNamespace("coda", quietly = TRUE)) {
    e <- try(coda::effectiveSize(samp), silent = TRUE)
    if (!inherits(e, "try-error")) ess <- unname(e[dc])
    if (n_chains > 1) {
      g <- try(coda::gelman.diag(samp, multivariate = FALSE)$psrf[, 1],
               silent = TRUE)
      if (!inherits(g, "try-error")) rh <- unname(g[dc])
    }
  }

  out <- data.frame(
    cgr = grid, method = if (likelihood == "t") "jags-t" else "jags",
    est = unname(colMeans(dr)),
    sd  = unname(apply(dr, 2, stats::sd)),
    lo  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .025)))),
    hi  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .975)))),
    p_fav = unname(apply(dr, 2, function(x) mean(direction * x > 0))),
    ess = ess, rhat = rh, stringsAsFactors = FALSE
  )
  out$mcse <- out$sd / sqrt(out$ess)   # autocorrelation-corrected
  # For the Student-t likelihood, expose the posterior-mean estimated degrees
  # of freedom: small nu (say < 10) means genuinely heavy tails and the robust
  # fit is doing work; large nu means it collapsed back toward the normal.
  if ("nu" %in% colnames(all)) attr(out, "nu") <- mean(all[, "nu"])
  out
}

#' Check that the conjugate and JAGS backends agree
#'
#' Runs both backends over the same grid and reports the largest absolute
#' differences in posterior means, credible limits and `P(favourable)`, the
#' largest difference in means measured in Monte Carlo standard errors
#' (`max_abs_z`), and JAGS convergence diagnostics. The verdict is `PASS` when
#' `max_abs_z < 5`.
#'
#' `cgr_unknown_check_backends()` does the same for the six-stratum UNKNOWN
#' extension.
#'
#' @param df A data frame with columns `condition`, `guess` and `value`, coded
#'   as for [cgr_strata()] (binary) or [cgr_unknown_strata()] (UNKNOWN).
#' @inheritParams cgr_conjugate
#' @param n_draws Conjugate posterior draws per stratum.
#' @param u_target Target UNKNOWN rate. `NULL` uses the observed rate.
#' @param jags_args Named list of extra arguments for [cgr_jags()] or
#'   [cgr_unknown_jags()].
#' @return A list of agreement statistics and a `verdict` string.
cgr_check_backends <- function(df, grid = seq(0, 1, length.out = 101),
                               n_draws = 40000, jags_args = list()) {
  a <- cgr_conjugate(df, grid, n_draws = n_draws)
  b <- do.call(cgr_jags, c(list(df = df, grid = grid), jags_args))
  a <- a[order(a$cgr), ]; b <- b[order(b$cgr), ]
  mcse <- sqrt(a$mcse^2 + b$mcse^2)
  res <- list(
    max_abs_diff_mean = max(abs(a$est   - b$est)),
    max_abs_diff_lo   = max(abs(a$lo    - b$lo)),
    max_abs_diff_hi   = max(abs(a$hi    - b$hi)),
    max_abs_diff_pfav = max(abs(a$p_fav - b$p_fav)),
    max_abs_z         = max(abs((a$est - b$est) / mcse)),
    max_rhat          = max(b$rhat, na.rm = TRUE),
    min_ess           = min(b$ess,  na.rm = TRUE)
  )
  res$verdict <- if (is.finite(res$max_abs_z) && res$max_abs_z < 5) {
    "PASS - differences within Monte Carlo error; identity claim supported"
  } else {
    "FAIL - remove the identity claim; describe as comparable weak-prior models"
  }
  res
}

#' JAGS backend for the UNKNOWN-preserving extension
#'
#' Six-stratum analogue of [cgr_jags()]. The six weights are computed inside
#' the model from the within-class shares `r`, `s` and `t` and the target
#' UNKNOWN rate `u`, all passed as data. Strata are indexed in
#' [UNKNOWN_STRATA] order.
#'
#' `pooling = "none"` (the default) uses independent strata and matches
#' [cgr_unknown_conjugate()]. `pooling = "partial"` is a sensitivity analysis
#' that partially pools the six stratum means toward a global mean with a
#' learned spread. It does not match the conjugate posterior.
#'
#' @param df A data frame with one row per participant and columns
#'   `condition` (`"AC"` or `"PL"`), `guess` (`"AC"`, `"PL"` or `"UNKNOWN"`)
#'   and numeric `value`.
#' @inheritParams jags_model_string
#' @param u_target Target UNKNOWN rate. `NULL` uses the observed rate.
#' @param pooling `"none"` or `"partial"`.
#' @param prior Named list with `m0`, `k0`, `a0` and `b0`. For partial pooling,
#'   optionally also `a_mu` and `b_mu`, the gamma prior on the pooling
#'   precision.
#' @return As for [cgr_jags()], with attributes `"u"` and `"pooling"`.
jags_unknown_model_string <- function(likelihood = c("normal", "t"),
                                      pooling = c("none", "partial")) {
  likelihood <- match.arg(likelihood); pooling <- match.arg(pooling)
  lik <- if (likelihood == "normal") {
    "  for (i in 1:N) { y[i] ~ dnorm(mu[k[i]], tau[k[i]]) }\n"
  } else {
    paste0("  for (i in 1:N) { y[i] ~ dt(mu[k[i]], tau[k[i]], nu) }\n",
           "  nu ~ dexp(0.1) T(2, 100)\n")
  }
  mu_prior <- if (pooling == "partial") {
    paste0(
      "  for (j in 1:6) {\n",
      "    tau[j] ~ dgamma(a0, b0)\n",
      "    mu[j]  ~ dnorm(mu0, prec_mu)\n",       # partial pooling toward mu0
      "  }\n",
      "  mu0 ~ dnorm(m0, 1.0E-6)\n",
      "  prec_mu ~ dgamma(a_mu, b_mu)\n")
  } else {
    paste0(
      "  for (j in 1:6) {\n",
      "    tau[j] ~ dgamma(a0, b0)\n",
      "    mu[j]  ~ dnorm(m0, k0 * tau[j])\n",
      "  }\n")
  }
  paste0(
    "model {\n", lik, mu_prior,
    "  for (m in 1:M) {\n",
    "    w_acac[m] <- (1 - u) * cgr[m] * (1 - r)\n",
    "    w_acpl[m] <- (1 - u) * (1 - cgr[m]) * s\n",
    "    w_acu[m]  <- u * t\n",
    "    w_plac[m] <- (1 - u) * (1 - cgr[m]) * (1 - s)\n",
    "    w_plpl[m] <- (1 - u) * cgr[m] * r\n",
    "    w_plu[m]  <- u * (1 - t)\n",
    "    delta[m] <- (w_acac[m]*mu[1] + w_acpl[m]*mu[2] + w_acu[m]*mu[3]) /\n",
    "                (w_acac[m] + w_acpl[m] + w_acu[m]) -\n",
    "                (w_plac[m]*mu[4] + w_plpl[m]*mu[5] + w_plu[m]*mu[6]) /\n",
    "                (w_plac[m] + w_plpl[m] + w_plu[m])\n",
    "  }\n}\n")
}

#' @rdname jags_unknown_model_string
cgr_unknown_jags <- function(df, grid = seq(0, 1, length.out = 101),
                             u_target = NULL, likelihood = c("normal", "t"),
                             pooling = c("none", "partial"),
                             n_iter = 10000, n_burn = 2000, n_chains = 4, seed = 1,
                             prior = list(m0 = 0, k0 = 1e-6, a0 = 1e-3, b0 = 1e-3),
                             direction = 1) {
  likelihood <- match.arg(likelihood); pooling <- match.arg(pooling)
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("needs rjags + a JAGS install; see https://mcmc-jags.sourceforge.io",
         call. = FALSE)
  }
  st  <- cgr_unknown_strata(df)
  o   <- cgr_unknown_observed(st)
  u   <- if (is.null(u_target)) o$u_obs else u_target
  cgr_unknown_estimable(st, u, grid, warn = FALSE)
  rat <- cgr_unknown_ratios(st)
  z   <- function(v) if (is.na(v)) 0 else v      # empty-class share -> 0 weight
  eps <- 1e-6
  grid_safe <- pmin(pmax(grid, eps), 1 - eps)

  # empty strata carry no data; their mu[j] is a prior draw multiplied by an
  # exactly-zero weight, so it never enters delta.
  y <- unlist(st[UNKNOWN_STRATA], use.names = FALSE)
  k <- rep(seq_along(UNKNOWN_STRATA), lengths(st[UNKNOWN_STRATA]))

  dat <- c(list(y = y, k = k, N = length(y), cgr = grid_safe, M = length(grid_safe),
                r = z(rat$r), s = z(rat$s), t = z(rat$t), u = u), prior)
  if (pooling == "partial") {
    dat[["k0"]] <- NULL                    # unused by the pooled prior; avoids a JAGS note
    dat <- c(dat, list(a_mu = if (is.null(prior$a_mu)) 1e-3 else prior$a_mu,
                       b_mu = if (is.null(prior$b_mu)) 1e-3 else prior$b_mu))
  }
  inits <- lapply(seq_len(n_chains), function(i)
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + i))

  con <- textConnection(jags_unknown_model_string(likelihood, pooling))
  on.exit(close(con), add = TRUE)
  m <- rjags::jags.model(con, data = dat, inits = inits,
                         n.chains = n_chains, quiet = TRUE)
  stats::update(m, n_burn)
  mon  <- if (likelihood == "t") c("delta", "nu") else "delta"
  samp <- rjags::coda.samples(m, variable.names = mon, n.iter = n_iter)

  all <- as.matrix(samp)
  dc <- grep("^delta\\[", colnames(all), value = TRUE)
  dc <- dc[order(as.integer(sub("^delta\\[(\\d+)\\]$", "\\1", dc)))]
  dr <- all[, dc, drop = FALSE]

  ess <- rep(NA_real_, length(dc)); rh <- ess
  if (requireNamespace("coda", quietly = TRUE)) {
    e <- try(coda::effectiveSize(samp), silent = TRUE)
    if (!inherits(e, "try-error")) ess <- unname(e[dc])
    if (n_chains > 1) {
      g <- try(coda::gelman.diag(samp, multivariate = FALSE)$psrf[, 1], silent = TRUE)
      if (!inherits(g, "try-error")) rh <- unname(g[dc])
    }
  }
  meth <- paste0(if (likelihood == "t") "jags-t" else "jags",
                 if (pooling == "partial") "-pooled" else "")
  out <- data.frame(
    cgr = grid, method = meth,
    est = unname(colMeans(dr)), sd = unname(apply(dr, 2, stats::sd)),
    lo  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .025)))),
    hi  = unname(apply(dr, 2, function(x) unname(stats::quantile(x, .975)))),
    p_fav = unname(apply(dr, 2, function(x) mean(direction * x > 0))),
    ess = ess, rhat = rh, stringsAsFactors = FALSE)
  out$mcse <- out$sd / sqrt(out$ess)
  attr(out, "u") <- u
  attr(out, "pooling") <- pooling
  if ("nu" %in% colnames(all)) attr(out, "nu") <- mean(all[, "nu"])
  out
}

#' @rdname cgr_check_backends
cgr_unknown_check_backends <- function(df, grid = seq(0, 1, length.out = 101),
                                       u_target = NULL, n_draws = 40000,
                                       jags_args = list()) {
  a <- cgr_unknown_conjugate(df, grid, u_target = u_target, n_draws = n_draws)
  b <- do.call(cgr_unknown_jags,
               c(list(df = df, grid = grid, u_target = u_target), jags_args))
  a <- a[order(a$cgr), ]; b <- b[order(b$cgr), ]
  mcse <- sqrt(a$mcse^2 + b$mcse^2)
  res <- list(
    max_abs_diff_mean = max(abs(a$est   - b$est)),
    max_abs_diff_lo   = max(abs(a$lo    - b$lo)),
    max_abs_diff_hi   = max(abs(a$hi    - b$hi)),
    max_abs_diff_pfav = max(abs(a$p_fav - b$p_fav)),
    max_abs_z         = max(abs((a$est - b$est) / mcse)),
    max_rhat          = max(b$rhat, na.rm = TRUE),
    min_ess           = min(b$ess,  na.rm = TRUE))
  res$verdict <- if (is.finite(res$max_abs_z) && res$max_abs_z < 5) {
    "PASS - UNKNOWN-extension backends agree within Monte Carlo error"
  } else {
    "FAIL - do not claim backend identity for the UNKNOWN extension"
  }
  res
}
