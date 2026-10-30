# Assumptions

Every assumption the analysis rests on, grouped by where it enters. "Testable
here" means testable with these data. Section numbers refer to the
[detailed overview](detailed_overview.md); U-numbers refer to
[open questions](open-questions.md).

## A. Assumptions of the CGRC estimand (inherited, not introduced)

| # | Assumption | Consequence if false | Testable here |
|---|---|---|---|
| A1 | Unblinding is malicious: perceived treatment causes expectancy, which causes the outcome. | If unblinding is benign (improvement causes the correct guess), adjustment removes a real effect: about a third of it when a quarter of guesses are benign, close to two-thirds from half onward ([stress test](benign-unblinding.md)). The original paper warns about this explicitly. | No (U6) |
| A2 | The four observed strata are exchangeable with the corresponding strata in a hypothetical well-blinded trial. | Reweighting transports the wrong population. Someone who guesses right under poor blinding may differ from someone who guesses right by luck. | No |
| A3 | The within-class arm shares r and s are the right quantities to hold fixed. | Overall target arm mass changes with the class masses, and each arm mean is renormalised. Whether the preserved shares are scientifically appropriate is an assumption. | Partly (Section 11) |
| A4 | r and s are known constants, not estimates. | Uncertainty is understated. Treating them as Beta-distributed widens the interval (Section 11). | Yes |
| A5 | Only stratum means matter; the estimand is linear in them. | Questions about tails or quantiles need a different estimator. | By construction |
| A6 | Treatment belief is captured by one binary guess. | A confident correct guess and a lucky one are treated the same. Confidence ratings would separate them. | No |
| A7 | Expectancy adds linearly to the direct treatment effect. | The AEB model assumes additivity; the paper notes that effects may not be additive. | No |
| A8 | CGR = 0.50 is the right target. | Szigeti et al. motivate chance-level guessing as an approximation to effective blinding; the reweighting alone does not identify a causal, perfectly blinded effect. | No |

## B. Assumptions introduced by the Bayesian implementation

| # | Assumption | Consequence if false | Testable here |
|---|---|---|---|
| B1 | Outcomes within a stratum are normal. | Posterior means are robust; intervals less so. Skew and kurtosis are reported, and a Student-t likelihood gives ν ≈ 18 for PANAS (Section 11, U8). | Yes |
| B2 | Observations are independent within a stratum. | Restricting to week 1 secures this. Pooling weeks would break it through repeated measures. | By design |
| B3 | Strata have independent parameters (no pooling). | No borrowing of strength between strata, which is conservative for small strata. A partial-pooling sensitivity exists for the UNKNOWN extension (`cgr_unknown_jags(pooling = "partial")`) but not for the binary model. | Yes, not run for the binary model |
| B4 | The priors are vague enough not to matter. | Vague and weakly informative priors agree; a strongly informative prior shifts the estimate, as it should (Section 11). | Yes |
| B5 | Each stratum has its own variance. | No equal-variance assumption is imposed, which is weaker than the pooled-variance t-test in the original. | By construction |
| B6 | `n_draws` is large enough that Monte Carlo error in the summaries is negligible. | The Monte Carlo standard error is reported with every estimate. | Yes |
| B7 | The conjugate and JAGS posteriors coincide. | Verified: the backends agree within Monte Carlo error across the full grid (Section 8.9, U1). | Yes |

## C. Assumptions about the data

| # | Assumption | Status |
|---|---|---|
| C1 | Week 1 only, so each data point is independent. | Follows the paper. |
| C2 | The raw value is the outcome, with no baseline adjustment. | Forced: the public data have no baseline. This reproduces the published curves. |
| C3 | Condition MD/PL maps to active/placebo. | Verified in the data audit. |
| C4 | Guess MD/PL maps to guessed active/guessed placebo. | Verified. |
| C5 | The public file matches the analysed file. | Questionable: n = 232 against the published 233 (U2). |
| C6 | The observed CGR is 0.647, not 0.72. | Verified from the data. 0.72 is a hardcoded constant in the author's code that matches the placebo-arm rate (U3, [sources](sources.md) S4). |
| C7 | No missing outcome values in week 1. | Verified: no blanks on the four scales analysed. |
| C8 | No duplicate participants in week 1. | Verified: 232 rows, 232 unique trial IDs. |

## D. Assumptions in the simulation study

| # | Assumption | Status |
|---|---|---|
| D1 | The AEB model is the right data-generating process. | It is the paper's model. Operating characteristics are conditional on it and would differ under a misspecified model. |
| D2 | `noise = "all"` is the intended reading of Eq. 4. | Not confirmed from source code. It is used because it most closely reproduces Hedges' g and Table 1 (U4). |
| D3 | The true Δ(0.5) equals the direct treatment effect. | Follows from the model: at CGR 0.5 the guess distribution is the same in both arms, so the expectancy term cancels. |
| D4 | 500 trials per scenario are enough. | The Monte Carlo SE of a coverage estimate is about sqrt(0.05 × 0.95 / 500) = 0.010, so 0.936 against 0.95 is within noise. |

## E. Claims deliberately not made

- That the Bayesian method gives a different or better estimate. Both
  procedures target the same estimand (Section 8).
- That the Energy VAS result is a direct drug effect.
- That attenuation measures how much of an effect was caused by expectancy.
- That reproducing the published curve validates any causal assumption.
- That every result is robust to every modelling choice. The Student-t check
  supports the normal likelihood for PANAS only.
- That the AEB simulation reproduces the author's simulation module line by
  line. The r and s formulas are confirmed from source; the Eq. 4 noise scope is
  not.
