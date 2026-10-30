# Open questions

Discrepancies and untested assumptions, ranked by how much they could change a
conclusion. Items resolved during the reproduction are listed at the end with
the evidence that resolved them. Section numbers refer to the
[detailed overview](detailed_overview.md); S1–S4 refer to [sources](sources.md).

## Open

### U2. n = 232 versus the published n = 233

The paper states the analysis rule verbatim: "we only used data from the first
week of the experiment, thus, each datapoint is independent... In the current
analysis n = 233 datapoints were included" (S1). That is the `tp == "w1s1"`
filter used here, so the rule is confirmed. The public file nonetheless gives
232 on every acute scale, with no duplicate trial IDs and no missing outcomes.
The papers do not account for the extra record. The most likely explanations
are a record withheld from the public release or a typo. Needs the original
author.

### U11. The published CGR-adjusted rejection rates are not reproduced

The unadjusted Table 1 rates reproduce, but the CGR-adjusted ones do not: the
KDE procedure flags the effect at CGR 0.50 far less often than published. The
obvious suspect is Eq. 4, which does not say whether the DTE and AEB noise terms
apply to the whole sample (`"all"`) or only to the affected subgroup (`"arm"`).
Trying every combination (300 trials per cell, n = 230, CGR 0.7, 32
source-faithful resamples, flag = averaged p < 0.05):

| DTE noise | AEB noise | Unadjusted, 4 scenarios | KDE adjusted, DTE only | KDE adjusted, DTE + AEB |
|---|---|---|---|---|
| all | all | 0.05 / 0.88 / 0.80 / 0.99 | 0.63 | 0.36 |
| all | arm | 0.04 / 0.88 / 0.93 / 1.00 | 0.64 | 0.43 |
| arm | all | 0.06 / 0.96 / 0.78 / 1.00 | 0.81 | 0.46 |
| arm | arm | 0.04 / 0.95 / 0.93 / 1.00 | 0.80 | 0.57 |
| **Published** | | **0.05 / 0.86 / 0.78 / 0.99** | **0.84** | **0.82** |

No reading reproduces all six numbers. Gating the DTE noise to the active arm
matches the adjusted DTE-only rate but overshoots the unadjusted one (0.96
against 0.86). The published adjusted rate with expectancy on (0.82) is out of
reach under every reading: adjusting for expectancy should leave power at or
below the DTE-only rate, since the expectancy term adds noise. Reducing the
expectancy noise does not help either, because it pushes the unadjusted
expectancy-only rate to 0.96 or above against the published 0.78. The gap
therefore lies outside the Eq. 4 noise scope: in the simulation parameters, the
adjustment step, or the decision rule. Needs the author's simulation code.

### U9. Guess rates run opposite to the pattern the author reports

Szigeti and Heifets (2024, S3) state that correct-guess rates are "generally
higher in the active arms". In this dataset the pattern is reversed: the
placebo arm guesses correctly at 0.723 and the microdose arm at 0.528
(Section 10). It is driven by a response bias toward guessing "placebo" (9,508
placebo versus 6,115 microdose guesses across the dataset). This does not break
the estimand, which conditions on the observed strata through r and s, but it
changes the intuition for which arm expectancy inflates. Worth raising with the
author.

### U6. Benign versus malicious unblinding cannot be tested with these data

If participants guess correctly *because* they improved, CGR adjustment removes
a real effect. The paper's argument for malicious unblinding is specific: 55% of
participants named body or perceptual sensations as their main cue (muscle
tension 58%, stomach discomfort 27%) against 23% naming mental or psychological
benefits, and the placebo-microdose PANAS differences (2.1 positive, 0.8
negative) are far smaller than natural day-to-day variability (about 10 and 6)
(S1). The 2024 review adds that microdose effects fall below the 0.5-SMD
minimal important difference (S3). That is an argument from evidence, not a
demonstration. No analysis of these data can settle it, so this stays open on
principle.

### U10. Assumptions behind the UNKNOWN-preserving extension

The six-stratum estimator is an addition by this package, not part of the
Szigeti estimand. Two interpretive assumptions remain open by nature:

- **The preserved UNKNOWN arm share t is an assumption.** Holding
  `t = ACU / (ACU + PLU)` fixed preserves the observed within-UNKNOWN arm ratio,
  exactly as the CGRC holds r and s fixed. It does not show that assignment is
  independent of the UNKNOWN response. `cgr_unknown_independent()` targets a
  different assumption and is offered as a contrast, not a resolution.
- **c = 0.50 is not perfect blinding.** It is directional guessing at chance
  with the UNKNOWN rate held fixed. Whether an UNKNOWN answer signals good
  blinding or disengagement cannot be decided from the counts.

Operating characteristics have been simulated with `sim_aeb_unknown()` and
`cgr_unknown_operating()` under two explicit assumptions: (A1) the UNKNOWN rate
is equal in both arms, and (A2) an UNKNOWN responder carries no expectancy.
With n = 300, u = 0.2 and 500 trials per scenario (`seed = 1`), the estimator
of Δ(0.50, u_obs) has bias of at most 0.08 against a true effect of 3, 95%
coverage of 0.950–0.960 in all four scenarios, and a pure-expectancy
false-favourable rate of 0.04, against 0.73 for a two-sided t-test. A precomputed lookup over n, p_cg, effect size and
u drives the app's design panel.

Still conditional: this is validation under A1 and A2 only. Differential
UNKNOWN rates (A1 relaxed) and UNKNOWN responders who carry expectancy (A2
relaxed) are not yet characterised. Six strata are thinner than four, so
`empty_stratum_rate` climbs faster at small n or high guess or UNKNOWN rates.
The UNKNOWN lookup fixes the expectancy effect at 7.7; the adjusted estimator is
nearly insensitive to it under A2, but the unadjusted comparator is not.

### U7. The KDE bandwidth was never a deliberate choice

The original uses scikit-learn's default fixed bandwidth of 1.0, not a
data-adaptive rule. On a 0–100 VAS that is very narrow. It does not matter here,
because the estimand uses only means and KDE smoothing preserves them, but it
would matter for any quantile- or tail-based estimand.

## Resolved

### U1. Agreement between the conjugate and JAGS backends

The JAGS prior was rewritten to the Normal-Inverse-Gamma form
`mu[j] ~ dnorm(m0, k0 * tau[j])`, so both backends target the same posterior.
`cgr_check_backends()` on the PANAS data (full 101-point grid; conjugate
backend at 40,000 draws; four JAGS chains at 10,000 iterations each; JAGS 4.3.2;
`set.seed(1)`) gives posterior means that agree within 0.008, 95% credible
limits within 0.035 (lower) and 0.047 (upper), and posterior probabilities
within 0.006. The largest mean difference is 0.74 Monte Carlo standard errors,
and convergence is clean (maximum R-hat 1.00, minimum ESS 39,664). Exact figures
vary slightly between runs. See Section 8.9. If the check
ever fails, the identity claim should be dropped rather than the prior weakened
to force agreement.

### U3. The quoted CGR of 0.72

The author's code (S4) sets `trial_cgrs = {'sbmd': 0.72}` in `src/config.py`, a
hardcoded constant drawn as the Figure 4 reference line. The same code computes
the trial CGR from the data as `(n_plpl + n_acac) / n` = 0.647 but does not use it
for the line. 0.72 matches the placebo-arm correct-guess rate (0.7234). The
reference-line test (Section 8.6) confirms that the curve equals the reported
unadjusted values at 0.647, not at 0.72. Whether the constant was meant as the
placebo-arm figure is a question only the author can answer; it does not affect
any estimate.

### U4. The author's source code

The code is named in the 2023 paper's data-availability statement:
`github.com/szb37/CorrectGuessRateCurve` (S4). Its `get_strata_ratio` and
`get_strata_sample_sizes` form `r = PLPL / (PLPL + ACAC)` and
`s = ACPL / (ACPL + PLAC)`, identical to this implementation, which independently
confirms the estimand. Reading it also showed that the resample count in the
repository's configuration is 32 over a 13-point grid, not 100, and that
`get_strata_ratio` rounds proportions to two decimals; the source-faithful KDE
mode reproduces that rounding. The Eq. 4 noise scope (`noise = "all"`) has not
been confirmed from the simulation module; it is used because it most closely
reproduces the published Table 1 and Hedges' g.

### U5. Head-to-head comparison of the KDE and Bayesian procedures

Section 8.8 compares the two procedures over 300 simulated trials with a true
effect of 3. They agree on the estimand: per-trial point estimates correlate at
about 0.996 with the same small bias. They differ on the evidence: at matched
thresholds the posterior flags the effect in about 74% of trials against about
54% for the averaged p-value. Setting the KDE bandwidth to zero leaves the flag
rate unchanged, so the gap comes from averaging p-values, not from the kernel.

### U8. Stratum normality

Skewness and excess kurtosis are reported per stratum in Section 11. A robust
Student-t likelihood (`cgr_jags(likelihood = "t")`) estimates ν ≈ 18 for PANAS,
so the t collapses toward the normal and the Gaussian conclusion holds despite
the skew in the ACAC stratum. Running it requires JAGS.
