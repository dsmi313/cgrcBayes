# Support functions for the Shiny design app (inst/app/app.R). These are UI,
# interpolation and input-normalisation helpers ONLY. Every statistical quantity
# still comes from the package's cgr_* functions - the app never carries its own
# copy of the estimand.

#' Precomputed design lookups
#'
#' The app reads operating characteristics from grids precomputed with
#' [cgr_operating()] and [cgr_unknown_operating()] (built by the scripts in
#' `data-raw/`), so it can respond instantly.
#'
#' * `cgrc_lookup()` and `cgrc_unknown_lookup()` load the binary and UNKNOWN
#'   grids. The UNKNOWN accessor returns `NULL` if its grid is not installed.
#' * `cgrc_op_at()` and `cgrc_unknown_op_at()` return the operating
#'   characteristics for one DTE by AEB scenario, interpolated bilinearly in
#'   `n` and `p_cg`. `true_effect`, `mu_aeb` and `u` snap to the nearest grid
#'   level. Flags report whether the point was interpolated, or clamped to the
#'   edge of the grid.
#' * `cgrc_power_curve()` and `cgrc_unknown_power_curve()` return the power of
#'   the adjusted analysis over `n`, from the DTE-on, AEB-off scenario.
#' * `cgrc_verdict()` and `cgrc_unknown_verdict()` write a plain-language
#'   readout of the trade-off between power and false treatment attribution,
#'   using `P(favourable) > 0.95` for both the adjusted and the unadjusted
#'   analysis. They do not declare a design safe.
#'
#' @param lut A lookup table from `cgrc_lookup()` or `cgrc_unknown_lookup()`.
#' @param n Planned number of participants.
#' @param p_cg Expected (directional) correct-guess rate.
#' @param true_effect True direct treatment effect, in outcome units.
#' @param u Expected UNKNOWN-response rate.
#' @param dte,aeb Scenario: `1` or `0` for the direct effect and for the
#'   expectancy effect.
#' @param mu_aeb Size of the expectancy effect, in outcome units.
#' @return `cgrc_lookup()` and `cgrc_unknown_lookup()`: data frames.
#'   `cgrc_op_at()` and `cgrc_unknown_op_at()`: a one-row data frame. The power
#'   curves: a data frame with `n` and `power`. The verdicts: a character
#'   string.
#' @examples
#' lut <- cgrc_lookup()
#' cgrc_op_at(lut, n = 120, p_cg = 0.8, true_effect = 3, dte = 1, aeb = 1)
#' cgrc_verdict(lut, n = 120, p_cg = 0.8, true_effect = 3)
cgrc_lookup <- function() {
  f <- system.file("extdata", "cgrc_lookup.rds", package = "cgrc.bayes")
  if (!nzchar(f) || !file.exists(f)) {
    stop("lookup table not found; run data-raw/build_lookup.R then reinstall.",
         call. = FALSE)
  }
  readRDS(f)
}

#' Launch the interactive design and analysis app
#'
#' Opens the Shiny app. Panel A (Design) answers "how reliable is CGR
#' adjustment for my trial design?" from a precomputed simulation grid. Panel B
#' (Analyse) takes a CSV of your own trial, audits it, and runs the adjustment.
#'
#' Requires \pkg{shiny}, \pkg{ggplot2} and \pkg{scales}.
#'
#' @param ... Passed to [shiny::runApp()].
#' @return Called for its side effect. Does not return while the app runs.
cgrc_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("the app needs the 'shiny' package: install.packages('shiny')",
         call. = FALSE)
  }
  app_dir <- system.file("app", package = "cgrc.bayes")
  if (!nzchar(app_dir)) stop("app directory not found in the install", call. = FALSE)
  shiny::runApp(app_dir, ...)
}

# Bilinear interpolation of a single column over two grid axes.
.cgrc_bilin <- function(df, xcol, ycol, vcol, x, y) {
  xs <- sort(unique(df[[xcol]])); ys <- sort(unique(df[[ycol]]))
  x <- min(max(x, min(xs)), max(xs)); y <- min(max(y, min(ys)), max(ys))
  x0 <- max(xs[xs <= x]); x1 <- min(xs[xs >= x])
  y0 <- max(ys[ys <= y]); y1 <- min(ys[ys >= y])
  g  <- function(xx, yy) df[[vcol]][df[[xcol]] == xx & df[[ycol]] == yy][1]
  tx <- if (x1 == x0) 0 else (x - x0) / (x1 - x0)
  ty <- if (y1 == y0) 0 else (y - y0) / (y1 - y0)
  (1 - tx) * (1 - ty) * g(x0, y0) + tx * (1 - ty) * g(x1, y0) +
    (1 - tx) * ty * g(x0, y1) + tx * ty * g(x1, y1)
}

#' @rdname cgrc_lookup
cgrc_op_at <- function(lut, n, p_cg, true_effect, dte, aeb, mu_aeb = 7.7) {
  sub <- lut[lut$DTE == dte & lut$AEB == aeb, ]
  ma  <- NA_real_
  if ("mu_aeb" %in% names(sub)) {
    mas <- sort(unique(sub$mu_aeb)); ma <- mas[which.min(abs(mas - mu_aeb))]
    sub <- sub[sub$mu_aeb == ma, ]
  }
  effs <- sort(unique(sub$true_effect))
  eff  <- effs[which.min(abs(effs - true_effect))]
  sub  <- sub[sub$true_effect == eff, ]
  if (!nrow(sub)) stop("no lookup rows for that scenario", call. = FALSE)
  # newer lookups carry the unadjusted Bayesian flag and the (research-only) 0.975
  # level; include whichever are present so older tables still interpolate.
  metrics <- c("unadj_bias", "adj_bias", "adj_rmse", "coverage95",
               "p_fav_gt_95", "freq_sig", "empty_stratum_rate")
  for (extra in c("unadj_p_fav_gt_95", "p_fav_gt_975"))
    if (extra %in% names(sub)) metrics <- c(metrics, extra)
  vals <- vapply(metrics, function(m) .cgrc_bilin(sub, "n", "p_cg", m, n, p_cg),
                 numeric(1))
  on_grid <- n %in% sub$n && p_cg %in% sub$p_cg
  # clamped: the requested n or p_cg is OUTSIDE this level's grid, so an edge
  # value is being reported under the requested label (e.g. n=60 is not
  # simulated at mu_aeb != 7.7, so it reads the n=80 row). Flag it, do not hide.
  clamped <- n < min(sub$n) || n > max(sub$n) ||
             p_cg < min(sub$p_cg) || p_cg > max(sub$p_cg)
  cbind(data.frame(n = n, p_cg = p_cg, true_effect = eff, mu_aeb = ma,
                   DTE = dte, AEB = aeb,
                   interpolated = !on_grid, clamped = clamped),
        as.data.frame(as.list(vals)))
}

#' @rdname cgrc_lookup
cgrc_power_curve <- function(lut, p_cg, true_effect) {
  ns <- sort(unique(lut$n))
  data.frame(n = ns,
             power = vapply(ns, function(nn)
               cgrc_op_at(lut, nn, p_cg, true_effect, 1, 0)$p_fav_gt_95,
               numeric(1)))
}

#' @rdname cgrc_lookup
cgrc_unknown_lookup <- function() {
  f <- system.file("extdata", "cgrc_unknown_lookup.rds", package = "cgrc.bayes")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  readRDS(f)
}

#' @rdname cgrc_lookup
cgrc_unknown_op_at <- function(lut, n, p_cg, true_effect, u, dte, aeb) {
  sub <- lut[lut$DTE == dte & lut$AEB == aeb, ]
  us  <- sort(unique(sub$u)); uu <- us[which.min(abs(us - u))]
  sub <- sub[sub$u == uu, ]
  effs <- sort(unique(sub$true_effect)); eff <- effs[which.min(abs(effs - true_effect))]
  sub  <- sub[sub$true_effect == eff, ]
  if (!nrow(sub)) stop("no UNKNOWN-lookup rows for that scenario", call. = FALSE)
  metrics <- c("unadj_bias", "adj_bias", "adj_rmse", "coverage95",
               "p_fav_gt_95", "unadj_p_fav_gt_95", "p_fav_gt_975", "freq_sig",
               "empty_stratum_rate")
  metrics <- metrics[metrics %in% names(sub)]
  vals <- vapply(metrics, function(m) .cgrc_bilin(sub, "n", "p_cg", m, n, p_cg),
                 numeric(1))
  on_grid <- n %in% sub$n && p_cg %in% sub$p_cg
  clamped <- n < min(sub$n) || n > max(sub$n) ||
             p_cg < min(sub$p_cg) || p_cg > max(sub$p_cg)
  cbind(data.frame(n = n, p_cg = p_cg, true_effect = eff, u = uu,
                   DTE = dte, AEB = aeb,
                   interpolated = !on_grid, clamped = clamped),
        as.data.frame(as.list(vals)))
}

#' @rdname cgrc_lookup
cgrc_unknown_power_curve <- function(lut, p_cg, true_effect, u) {
  ns <- sort(unique(lut$n))
  data.frame(n = ns,
             power = vapply(ns, function(nn)
               cgrc_unknown_op_at(lut, nn, p_cg, true_effect, u, 1, 0)$p_fav_gt_95,
               numeric(1)))
}

#' @rdname cgrc_lookup
cgrc_unknown_verdict <- function(lut, n, p_cg, true_effect, u) {
  minstr <- cgr_unknown_min_stratum(n, p_cg, u)
  fp  <- cgrc_unknown_op_at(lut, n, p_cg, true_effect, u, 0, 1)   # pure expectancy
  feas <- if (minstr < 15) sprintf(
    paste0(" (The smallest of the six strata averages ~%.0f participants here, so ",
           "the estimand can be undefined for some trials - see feasibility.)"),
    minstr) else ""
  ub <- function(row) if (!is.null(row$unadj_p_fav_gt_95)) row$unadj_p_fav_gt_95 else row$freq_sig
  if (true_effect == 0) {
    sprintf(paste0("UNKNOWN-preserving design at n=%d, directional CGR %.0f%%, ",
      "UNKNOWN rate %.0f%%, no true effect, using the standard threshold posterior ",
      "P(favourable)>0.95 throughout: when the apparent signal is pure expectancy, ",
      "an unadjusted analysis flags a false effect in %.0f%% of simulated trials, ",
      "while the UNKNOWN-adjusted analysis flags it in only %.0f%%.%s"),
      n, 100 * p_cg, 100 * u, 100 * ub(fp), 100 * fp$p_fav_gt_95, feas)
  } else {
    pw <- cgrc_unknown_op_at(lut, n, p_cg, true_effect, u, 1, 0)   # clean power
    sprintf(paste0("UNKNOWN-preserving design at n=%d, directional CGR %.0f%%, ",
      "UNKNOWN rate %.0f%%, using the standard threshold posterior ",
      "P(favourable)>0.95 throughout: for a real %.1f-point effect the ",
      "UNKNOWN-adjusted analysis flags the effect in %.0f%% of simulated trials ",
      "(its power; the unadjusted analysis %.0f%%). When the apparent effect is ",
      "instead pure expectancy, an unadjusted analysis flags a false effect in ",
      "%.0f%% and the adjusted analysis in only %.0f%%.%s"),
      n, 100 * p_cg, 100 * u, true_effect, 100 * pw$p_fav_gt_95, 100 * ub(pw),
      100 * ub(fp), 100 * fp$p_fav_gt_95, feas)
  }
}

#' Normalise arm and guess codings
#'
#' Map common codings to `"AC"` and `"PL"` (and `"UNKNOWN"` for guesses),
#' erroring on values they do not recognise rather than guessing.
#'
#' Active codes: AC, MD, active, drug, microdose, treatment, true, t, 1, yes,
#' y. Placebo codes: PL, placebo, control, sham, false, f, 0, no, n. Matching
#' ignores case and surrounding whitespace.
#'
#' `cgrc_normalise_guess()` keeps `NA` and blank values as missing. With
#' `allow_unknown = TRUE` it also maps "I don't know" answers (unknown, unsure,
#' not sure, don't know, DK, IDK and similar, plus any `unknown_labels`) to
#' `"UNKNOWN"`. With `allow_unknown = FALSE` those answers are rejected, so an
#' UNKNOWN response cannot enter a binary analysis unnoticed.
#'
#' @param x A vector of codes.
#' @param what Description of `x` for error messages.
#' @param allow_unknown Recognise UNKNOWN responses.
#' @param unknown_labels Extra labels to treat as UNKNOWN.
#' @return A character vector.
#' @examples
#' cgrc_normalise_arm(c("drug", "placebo", "1"))
#' cgrc_normalise_guess(c("active", "I don't know", ""), allow_unknown = TRUE)
cgrc_normalise_arm <- function(x, what = "value") {
  s <- toupper(trimws(as.character(x)))
  active  <- c("AC", "MD", "ACTIVE", "DRUG", "MICRODOSE", "TREATMENT",
               "TRUE", "T", "1", "YES", "Y")
  placebo <- c("PL", "PLACEBO", "CONTROL", "SHAM", "FALSE", "F", "0", "NO", "N")
  out <- ifelse(s %in% active, "AC", ifelse(s %in% placebo, "PL", NA_character_))
  bad <- unique(as.character(x)[is.na(out)])
  if (length(bad)) {
    stop(sprintf(paste0("could not map %s level(s) %s to active/placebo; ",
                        "recode to AC/PL (active/placebo)."),
                 what, paste(shQuote(bad), collapse = ", ")), call. = FALSE)
  }
  out
}

#' @rdname cgrc_normalise_arm
cgrc_normalise_guess <- function(x, allow_unknown = FALSE, unknown_labels = NULL,
                                 what = "treatment guessed") {
  raw <- as.character(x)
  s <- toupper(trimws(raw))
  # normalise apostrophe variants to a plain ' with fixed byte replacement so this
  # works in a C locale too, where a multibyte regex character class fails. \u
  # escapes keep the source pure-ASCII; useBytes avoids any locale translation.
  for (ap in c("\u2019", "\u02bc", "\u2018", "`"))
    s <- gsub(ap, "'", s, fixed = TRUE, useBytes = TRUE)
  s <- gsub("[[:space:]]+", " ", s)             # collapse internal whitespace
  active  <- c("AC", "MD", "ACTIVE", "DRUG", "MICRODOSE", "TREATMENT",
               "TRUE", "T", "1", "YES", "Y")
  placebo <- c("PL", "PLACEBO", "CONTROL", "SHAM", "FALSE", "F", "0", "NO", "N")
  unknown <- c("UNKNOWN", "UNSURE", "UNCERTAIN", "NOT SURE", "DO NOT KNOW",
               "DON'T KNOW", "DONT KNOW", "DID NOT KNOW", "DIDN'T KNOW",
               "DIDNT KNOW", "I DO NOT KNOW", "I DON'T KNOW", "I DONT KNOW",
               "DK", "IDK", "NO IDEA", "CAN'T TELL", "CANT TELL",
               toupper(trimws(as.character(unknown_labels))))
  missing <- is.na(raw) | s == ""               # NA / blank stay missing, never UNKNOWN
  out <- rep(NA_character_, length(s))
  out[s %in% active]  <- "AC"
  out[s %in% placebo] <- "PL"
  if (allow_unknown) out[s %in% unknown] <- "UNKNOWN"
  out[missing] <- NA_character_
  bad <- unique(raw[!missing & is.na(out)])      # present, but mapped to nothing
  if (length(bad)) {
    hint <- if (allow_unknown) "AC/PL or an UNKNOWN response" else "AC/PL (active/placebo)"
    stop(sprintf("could not map %s level(s) %s to %s; recode them.",
                 what, paste(shQuote(bad), collapse = ", "), hint), call. = FALSE)
  }
  out
}

#' Audit an uploaded trial
#'
#' Classifies every row instead of silently dropping incomplete ones. Rows
#' missing a condition, guess or outcome, or with a non-numeric outcome, are
#' excluded and logged with the reason. An observed UNKNOWN guess is not
#' missing data and is kept. Codes that cannot be mapped at all are errors (see
#' [cgrc_normalise_arm()]).
#'
#' @param condition,guess,value Columns of the uploaded trial.
#' @param unknown_level Extra label to treat as an UNKNOWN guess.
#' @return A list with the `clean` analysis data frame, an exclusion `log`
#'   (row number and reason), a named vector of counts (`summary`), and
#'   `has_unknown`.
cgrc_input_audit <- function(condition, guess, value, unknown_level = "UNKNOWN") {
  n  <- length(condition)
  cr <- trimws(as.character(condition)); gr <- trimws(as.character(guess))
  vr <- as.character(value)
  miss_c <- is.na(condition) | cr == ""
  miss_g <- is.na(guess)     | gr == ""
  miss_v <- is.na(value)     | trimws(vr) == ""
  vnum   <- suppressWarnings(as.numeric(vr))
  nonnum <- !miss_v & is.na(vnum)

  cond <- rep(NA_character_, n); g <- rep(NA_character_, n)
  if (any(!miss_c))
    cond[!miss_c] <- cgrc_normalise_arm(condition[!miss_c], "treatment received")
  if (any(!miss_g))
    g[!miss_g] <- cgrc_normalise_guess(guess[!miss_g], allow_unknown = TRUE,
                                       unknown_labels = unknown_level)
  is_unknown <- !miss_g & !is.na(g) & g == "UNKNOWN"

  # exclusion reason, by precedence (a missing condition wins over a bad outcome)
  reason <- rep("ok", n)
  reason[nonnum] <- "non-numeric outcome"
  reason[miss_v] <- "missing outcome"
  reason[miss_g] <- "missing guess"
  reason[miss_c] <- "missing condition"
  excluded <- reason != "ok"

  summary <- c(n_input = n,
               missing_condition  = sum(miss_c),
               missing_guess      = sum(miss_g),
               observed_unknown   = sum(is_unknown),
               missing_outcome    = sum(miss_v),
               nonnumeric_outcome = sum(nonnum),
               excluded           = sum(excluded),
               n_analysis         = sum(!excluded))
  list(
    clean = data.frame(condition = cond[!excluded], guess = g[!excluded],
                       value = vnum[!excluded], stringsAsFactors = FALSE),
    log = data.frame(row = which(excluded), reason = reason[excluded],
                     stringsAsFactors = FALSE),
    summary = summary,
    has_unknown = sum(is_unknown) > 0)
}

#' Analysis reports for the app
#'
#' Build the app's downloadable outputs from its assembled analysis object.
#' They live in the package, not the app, so they are tested and the app
#' carries no statistics of its own.
#'
#' * `cgrc_analysis_summary()`: one row per quantity, shared by the on-screen
#'   summary and the summary CSV.
#' * `cgrc_build_report()`: a Markdown report with the input audit, stratum
#'   counts, estimates, the identity check, ROPE results, the seed and the
#'   package version.
#' * `cgrc_build_html_report()`: a compact, self-contained HTML report.
#'
#' @param f The analysis object assembled by the app (`inst/app/app.R`): a
#'   list with `mode` (`"binary"` or `"unknown"`), the cleaned `trial`, the
#'   `audit` from [cgrc_input_audit()], `delta`, `dir`, `seed`, and the fitted
#'   results (`fit`, `head` and `rope` in binary mode; `ufit`, `uhead` and
#'   `urope` in UNKNOWN mode).
#' @return `cgrc_analysis_summary()`: a data frame with columns `quantity` and
#'   `value`. `cgrc_build_report()`: a character vector of Markdown lines.
#'   `cgrc_build_html_report()`: a single HTML string.
cgrc_build_report <- function(f) {
  ver <- tryCatch(as.character(utils::packageVersion("cgrc.bayes")), error = function(e) "NA")
  s <- f$audit$summary
  L <- c()
  add <- function(...) L <<- c(L, sprintf(...))
  add("# CGRC analysis report")
  add("")
  add("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add("- Package: cgrc.bayes %s", ver)
  add("- Random seed: %s", if (is.null(f$seed)) "not set" else as.character(f$seed))
  meth <- if (f$mode == "unknown")
    paste("Proposed UNKNOWN-preserving CGRC extension (six strata; observed UNKNOWN rate",
          "held fixed while the directional correct-guess rate is varied). This",
          "is an extension implemented by cgrc.bayes, not the original Szigeti",
          "estimand.")
  else "Binary four-stratum CGRC (Normal-Inverse-Gamma conjugate posterior)."
  add("- Method: %s", meth)
  add("")
  add("## Input audit")
  add("")
  add("| Quantity | n |")
  add("|---|---|")
  add("| Rows uploaded | %d |", s[["n_input"]])
  add("| Analysed | %d |", nrow(f$trial))
  add("| Missing condition | %d |", s[["missing_condition"]])
  add("| Missing guess | %d |", s[["missing_guess"]])
  add("| Missing outcome | %d |", s[["missing_outcome"]])
  add("| Non-numeric outcome | %d |", s[["nonnumeric_outcome"]])
  add("| Observed UNKNOWN responses | %d |", s[["observed_unknown"]])
  if (isTRUE(f$n_excl_unknown > 0))
    add("| UNKNOWN excluded (complete-case) | %d |", f$n_excl_unknown)
  add("")

  if (f$mode == "unknown") {
    st <- cgr_unknown_strata(f$trial); n <- lengths(st)
    add("## Stratum counts (UNKNOWN preserved)")
    add("")
    add("| Stratum | n |")
    add("|---|---|")
    for (k in UNKNOWN_STRATA) add("| %s | %d |", k, n[[k]])
    add("")
    add("- Observed directional CGR: %.4f", f$ufit$observed_directional_cgr)
    add("- Observed UNKNOWN rate: %.4f (held fixed at %.4f)",
        f$ufit$observed_unknown_rate, f$ufit$target_unknown_rate)
    add("- n total / directional / UNKNOWN: %d / %d / %d",
        f$ufit$n_total, f$ufit$n_directional, f$ufit$n_unknown)
    dg <- cgr_unknown_estimable(st, f$ufit$target_unknown_rate,
                                c(0.5, f$ufit$observed_directional_cgr),
                                warn = FALSE)
    if (dg$boundary_share) add(paste0(
      "- **Fragility warning:** observed zero cell(s) ",
      paste(dg$boundary_cells, collapse = ", "),
      " imply a plug-in within-class arm share of 0 or 1."))
    sm <- f$ufit$summary; z <- cgr_unknown_reference_line_test(f$trial)
    h <- f$uhead
    obs_est <- sm$post_mean[1]; adj_est <- sm$post_mean[2]
  } else {
    st <- cgr_strata(f$trial)
    add("## Stratum counts")
    add("")
    add("| Stratum | n |")
    add("|---|---|")
    for (k in STRATA) add("| %s | %d |", k, length(st[[k]]))
    add("")
    add("- Observed CGR: %.4f", f$fit$observed_cgr)
    sm <- f$fit$summary; z <- cgr_reference_line_test(f$trial, f$fit$observed_cgr)
    h <- f$head
    obs_est <- sm$post_mean[1]; adj_est <- sm$post_mean[2]
  }
  add("")
  add("## Estimates")
  add("")
  add("| | posterior mean | 95%% CrI | P(favourable) |")
  add("|---|---|---|---|")
  add("| Raw (observed CGR) | %.3f | %.3f to %.3f | %.3f |",
      sm$post_mean[1], sm$cri_lo[1], sm$cri_hi[1], sm$p_favourable[1])
  add("| Adjusted (CGR 0.50) | %.3f | %.3f to %.3f | %.3f |",
      sm$post_mean[2], sm$cri_lo[2], sm$cri_hi[2], sm$p_favourable[2])
  add("")
  add("- Meaningful-difference threshold delta: %.3f outcome units", f$delta)
  add("- P(meaningful) raw / adjusted: %.3f / %.3f",
      h$p_meaningful_obs, h$p_meaningful_blind)
  add("")
  add("## Identity check")
  add("")
  add("At the observed CGR the reweighting is a no-op: Delta = raw arm-mean")
  add("difference to %.1e (should be ~0).", abs(z$D_at_obs - z$raw_mean_diff))
  add("")
  add("## Region of practical equivalence (at CGR 0.50)")
  add("")
  rope <- if (f$mode == "unknown") f$urope else f$rope
  i <- which.min(abs(rope$cgr - 0.5))
  add("- P(meaningful benefit): %.3f", rope$p_benefit[i])
  add("- P(practically negligible): %.3f", rope$p_negligible[i])
  add("- P(meaningful harm): %.3f", rope$p_harm[i])
  add("")
  add("_Report text: %s_", h$text)
  L
}

#' @rdname design_helpers
cgr_expected_strata <- function(n, p_cg) {
  c(ACAC = 0.5 * n * p_cg,       ACPL = 0.5 * n * (1 - p_cg),
    PLAC = 0.5 * n * (1 - p_cg), PLPL = 0.5 * n * p_cg)
}

#' Correct-guess rates by arm
#'
#' The overall correct-guess rate can hide an asymmetry between arms. This
#' returns the active-arm rate, P(guess active | received active), the
#' placebo-arm rate, P(guess placebo | received placebo), the overall rate,
#' allocation and guess counts, and the UNKNOWN-response rate in each arm
#' (`unknown_active`, `unknown_placebo`). An UNKNOWN guess is never counted as
#' correct. Unlike [cgr_strata()], it works when a stratum is empty.
#'
#' @param trial A data frame with `condition` (`"AC"` or `"PL"`) and `guess`
#'   (`"AC"` or `"PL"`, optionally `"UNKNOWN"`).
#' @return A named list.
cgr_guess_rates <- function(trial) {
  cond <- as.character(trial$condition); guess <- as.character(trial$guess)
  n_ac <- sum(cond == "AC"); n_pl <- sum(cond == "PL")
  correct <- cond == guess & guess %in% c("AC", "PL")
  list(
    overall  = if (length(cond)) mean(correct) else NA_real_,
    active   = if (n_ac) sum(cond == "AC" & guess == "AC") / n_ac else NA_real_,
    placebo  = if (n_pl) sum(cond == "PL" & guess == "PL") / n_pl else NA_real_,
    n_active = n_ac, n_placebo = n_pl,
    guess_active  = sum(guess == "AC"),
    guess_placebo = sum(guess == "PL"),
    guess_unknown = sum(guess == "UNKNOWN"),
    unknown_active  = if (n_ac) sum(cond == "AC" & guess == "UNKNOWN") / n_ac else NA_real_,
    unknown_placebo = if (n_pl) sum(cond == "PL" & guess == "UNKNOWN") / n_pl else NA_real_)
}

#' Feasibility category for a design
#'
#' A four-level category based only on the expected smallest stratum and the
#' simulated empty-stratum rate: "Adjustment frequently undefined", "Fragile
#' design", "Use with caution", or "Feasibility looks good under simulated
#' conditions". It does not account for bias, coverage or power, which are
#' reported separately.
#'
#' @param min_stratum Expected smallest stratum, for example from
#'   [cgr_min_stratum()].
#' @param empty_stratum_rate Simulated rate of trials with an empty stratum.
#' @param thin Smallest-stratum size below which a design is fragile. Below
#'   `2 * thin` it calls for caution.
#' @param degen_warn,degen_bad Empty-stratum rates above which a design is
#'   fragile, or frequently undefined.
#' @return A list with `category` and the CSS `class` the app uses for the
#'   badge.
cgrc_reliability <- function(min_stratum, empty_stratum_rate = NA_real_,
                             thin = 15, degen_warn = 0.02, degen_bad = 0.10) {
  degen <- if (is.na(empty_stratum_rate)) 0 else empty_stratum_rate
  if (degen > degen_bad) {
    list(category = "Adjustment frequently undefined", class = "warn")
  } else if (min_stratum < thin || degen > degen_warn) {
    list(category = "Fragile design", class = "warn")
  } else if (min_stratum < 2 * thin) {
    list(category = "Use with caution", class = "caution")
  } else {
    list(category = "Feasibility looks good under simulated conditions", class = "ok")
  }
}

#' @rdname cgrc_build_report
cgrc_analysis_summary <- function(f) {
  if (f$mode == "unknown") {
    sm <- f$ufit$summary; h <- f$uhead
    obs_cgr <- f$ufit$observed_directional_cgr
    st <- cgr_unknown_strata(f$trial); nn <- lengths(st); smallest <- min(nn[nn > 0])
  } else {
    sm <- f$fit$summary; h <- f$head
    obs_cgr <- f$fit$observed_cgr
    st <- cgr_strata(f$trial); smallest <- min(lengths(st)[STRATA])
  }
  data.frame(
    quantity = c(
      "Raw randomized treatment contrast (at observed CGR)",
      "CGRC-adjusted contrast (target CGR 0.50)",
      "Adjusted 95% CrI lower", "Adjusted 95% CrI upper",
      "P(favourable) raw", "P(favourable) adjusted",
      sprintf("P(exceeds %.3g-unit threshold) raw", f$delta),
      sprintf("P(exceeds %.3g-unit threshold) adjusted", f$delta),
      "Observed correct-guess rate", "Smallest stratum size",
      "Meaningful-difference threshold (outcome units)", "Random seed"),
    value = c(
      round(sm$post_mean[1], 3), round(sm$post_mean[2], 3),
      round(sm$cri_lo[2], 3), round(sm$cri_hi[2], 3),
      round(sm$p_favourable[1], 3), round(sm$p_favourable[2], 3),
      round(h$p_meaningful_obs, 3), round(h$p_meaningful_blind, 3),
      round(obs_cgr, 4), smallest, round(f$delta, 3),
      if (is.null(f$seed)) NA_real_ else f$seed),
    stringsAsFactors = FALSE)
}

#' @rdname cgrc_build_report
cgrc_build_html_report <- function(f) {
  ver <- tryCatch(as.character(utils::packageVersion("cgrc.bayes")),
                  error = function(e) "NA")
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  summ <- cgrc_analysis_summary(f)
  gr <- tryCatch(cgr_guess_rates(f$trial), error = function(e) NULL)
  method <- if (f$mode == "unknown")
    "Proposed UNKNOWN-preserving CGRC extension (six strata; observed UNKNOWN rate held fixed)."
    else "Binary four-stratum CGRC (Normal-Inverse-Gamma conjugate posterior)."
  s <- f$audit$summary
  warns <- character(0)
  smallest <- summ$value[summ$quantity == "Smallest stratum size"]
  if (length(smallest) && is.finite(smallest) && smallest < 15)
    warns <- c(warns, sprintf("Smallest stratum has %d participants (&lt; 15): the reweighted estimate is fragile.", as.integer(smallest)))
  if (!is.null(gr) && is.finite(gr$active) && is.finite(gr$placebo) &&
      abs(gr$active - gr$placebo) > 0.2)
    warns <- c(warns, sprintf("Arm asymmetry: active-arm correct-guess rate %.0f%% vs placebo-arm %.0f%%.",
                              100 * gr$active, 100 * gr$placebo))
  if (f$mode == "unknown") {
    st <- cgr_unknown_strata(f$trial)
    dg <- cgr_unknown_estimable(st, f$ufit$target_unknown_rate,
                                c(0.5, f$ufit$observed_directional_cgr),
                                warn = FALSE)
    if (dg$boundary_share) warns <- c(warns, paste0(
      "Observed zero cell(s) ", paste(dg$boundary_cells, collapse = ", "),
      " imply a plug-in within-class arm share of 0 or 1; the estimate is fragile."))
  }
  rows <- paste(sprintf("<tr><td>%s</td><td style='text-align:right'>%s</td></tr>",
                        esc(summ$quantity), ifelse(is.na(summ$value), "NA",
                          formatC(summ$value, format = "g", digits = 6))),
                collapse = "\n")
  warn_html <- if (length(warns))
    paste0("<div class='warn'><b>Warnings</b><ul>",
           paste(sprintf("<li>%s</li>", warns), collapse = ""), "</ul></div>")
    else "<p class='muted'>No feasibility warnings.</p>"
  paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>CGRC analysis report</title><style>",
    "body{font-family:system-ui,Arial,sans-serif;font-size:16px;max-width:820px;",
    "margin:24px auto;padding:0 16px;color:#222;line-height:1.5;}",
    "h1{font-size:24px;} h2{font-size:20px;margin-top:1.4em;}",
    "table{border-collapse:collapse;width:100%;} td,th{border:1px solid #ccc;padding:6px 8px;}",
    ".muted{color:#555;} .warn{background:#FDEDEC;border-left:5px solid #C0392B;",
    "padding:8px 12px;border-radius:4px;}</style></head><body>",
    "<h1>CGRC analysis report</h1>",
    "<p class='muted'>Generated ", esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    " &middot; cgrc.bayes ", esc(ver),
    " &middot; seed ", if (is.null(f$seed)) "not set" else esc(as.character(f$seed)),
    "</p>",
    "<h2>Settings</h2><ul>",
    "<li>Method: ", esc(method), "</li>",
    "<li>Favourable direction: ", if (isTRUE(f$dir < 0)) "lower is better" else "higher is better", "</li>",
    "<li>Meaningful-difference threshold: ", formatC(f$delta, format = "g", digits = 4), " outcome units</li>",
    "<li>Rows uploaded: ", s[["n_input"]], "; analysed: ", nrow(f$trial), "</li>",
    "</ul>",
    "<h2>Results</h2><table>", rows, "</table>",
    "<p class='muted'>The adjusted contrast is the CGRC reweighted estimand at ",
    "a 50% correct-guess rate, ",
    "not automatically the true pharmacological effect. Probabilities are ",
    "continuous &mdash; no significant/not-significant cut-off is imposed.</p>",
    "<h2>Warnings</h2>", warn_html,
    "</body></html>")
}

#' Should percentage attenuation be reported?
#'
#' Returns `FALSE` when the unadjusted 95% credible interval includes zero, or
#' when the adjusted and unadjusted estimates have opposite signs. In both
#' cases the ratio is not meaningful.
#'
#' @param unadj_est,unadj_lo,unadj_hi Unadjusted estimate and 95% credible
#'   limits.
#' @param adj_est Adjusted estimate.
#' @return `TRUE` or `FALSE`.
cgrc_pct_ok <- function(unadj_est, unadj_lo, unadj_hi, adj_est) {
  distinct  <- !(unadj_lo <= 0 && unadj_hi >= 0)
  same_sign <- sign(unadj_est) == sign(adj_est)
  distinct && same_sign
}

#' @rdname cgrc_lookup
cgrc_verdict <- function(lut, n, p_cg, true_effect, mu_aeb = 7.7) {
  minstr <- cgr_min_stratum(n, p_cg)
  fp  <- cgrc_op_at(lut, n, p_cg, true_effect, 0, 1, mu_aeb)   # pure expectancy
  feas <- if (minstr < 15) sprintf(
    paste0(" (The smallest stratum averages ~%.0f participants here, so the ",
           "estimand can be undefined for some trials - see the feasibility ",
           "readout.)"), minstr) else ""
  # unadjusted comparator is the Bayesian raw-contrast flag at the SAME P>0.95
  # threshold; on an older lookup without it, fall back to the t-test rate.
  ub <- function(row) if (!is.null(row$unadj_p_fav_gt_95)) row$unadj_p_fav_gt_95 else row$freq_sig
  if (true_effect == 0) {
    sprintf(paste0("At n=%d with %.0f%% correct guessing and no true effect, using ",
      "the standard threshold posterior P(favourable)>0.95 throughout: when the ",
      "apparent signal is pure expectancy, an unadjusted analysis flags a false ",
      "effect in %.0f%% of simulated trials, while the CGR-adjusted analysis flags ",
      "it in only %.0f%% - the false treatment attribution the adjustment removes.%s"),
      n, 100 * p_cg, 100 * ub(fp), 100 * fp$p_fav_gt_95, feas)
  } else {
    pw <- cgrc_op_at(lut, n, p_cg, true_effect, 1, 0, mu_aeb)  # clean power
    sprintf(paste0("At n=%d with %.0f%% correct guessing, using the standard ",
      "threshold posterior P(favourable)>0.95 throughout: for a real %.1f-point ",
      "effect the CGR-adjusted analysis flags the effect in %.0f%% of simulated ",
      "trials (its power; the unadjusted analysis %.0f%%). When the apparent effect ",
      "is instead pure expectancy, an unadjusted analysis flags a false effect in ",
      "%.0f%% and the adjusted analysis in only %.0f%%. Weigh that trade-off against ",
      "your own priorities.%s"),
      n, 100 * p_cg, true_effect, 100 * pw$p_fav_gt_95, 100 * ub(pw),
      100 * ub(fp), 100 * fp$p_fav_gt_95, feas)
  }
}
