# Stress test of the "I don't know" extension

The six-stratum estimator keeps UNKNOWN guesses as a third category. Its
operating characteristics were first simulated under two assumptions:

- **A1:** the UNKNOWN rate is the same in both arms.
- **A2:** an UNKNOWN responder carries no expectancy.

This page checks what happens when each assumption fails, and compares the
extension with the obvious alternative: dropping UNKNOWN responders and running
the binary CGRC on the rest (complete case).

## Design

`cgr_unknown_stress()` simulates trials with `sim_aeb_unknown()` and two extra
parameters:

- `u_diff`, the active-minus-placebo difference in UNKNOWN rates (A1 holds at
  0). The arms answer UNKNOWN at `u + u_diff / 2` and `u - u_diff / 2`.
- `aeb_unknown`, the share of the expectancy effect an UNKNOWN responder still
  carries through the belief they did not state (A2 holds at 0).

Settings: n = 300, directional correct-guess rate 0.7, average UNKNOWN rate
0.2, expectancy effect 7.7, 500 trials per cell. Each cell is run with no direct
effect (so the rate of `P(favourable) > 0.95` is a false-positive rate) and with
a direct effect of 3 (so it is power). Both methods are fitted to the same
posterior draws in every trial. Bias is averaged over the two scenarios.

Reproduce with `Rscript data-raw/stress_unknown.R`; the full results are in
[`results/unknown-stress.csv`](results/unknown-stress.csv).

## Results

### A2 fails: UNKNOWN responders carry expectancy

| Share of expectancy carried | Bias, six-stratum | Bias, complete case | False positives, six-stratum | False positives, complete case | Power, six-stratum | Power, complete case |
|---|---|---|---|---|---|---|
| 0 | 0.04 | 0.05 | 0.06 | 0.06 | 0.84 | 0.73 |
| 0.25 | 0.11 | −0.05 | 0.06 | 0.04 | 0.81 | 0.68 |
| 0.5 | 0.30 | −0.01 | 0.11 | 0.04 | 0.87 | 0.70 |
| 1 | 0.62 | −0.00 | 0.18 | 0.07 | 0.88 | 0.69 |

### A1 fails: UNKNOWN rates differ between arms

| Difference in UNKNOWN rates (active − placebo) | Bias, six-stratum | Bias, complete case | False positives, six-stratum | False positives, complete case | Power, six-stratum | Power, complete case |
|---|---|---|---|---|---|---|
| −0.2 | 0.77 | 0.00 | 0.21 | 0.04 | 0.93 | 0.69 |
| −0.1 | 0.40 | 0.01 | 0.12 | 0.06 | 0.89 | 0.72 |
| 0 | 0.04 | 0.05 | 0.06 | 0.06 | 0.84 | 0.73 |
| 0.1 | −0.46 | −0.05 | 0.02 | 0.05 | 0.69 | 0.69 |
| 0.2 | −0.75 | 0.05 | 0.01 | 0.04 | 0.61 | 0.72 |

These agree with the separate A1 check in Section 6.2 of the
[summary paper](summary.md) (−0.82 at a 20-point gap, −0.36 at a 10-point gap,
at n = 230).

Across the full grid, which also crosses the two violations, 95% coverage
ranges from 0.86 to 0.97 for the six-stratum estimator and from 0.93 to 0.97
for complete case. The largest complete-case bias is 0.13, and its
false-positive rate never exceeds 0.07.

## What this shows

**When A1 and A2 hold, the extension earns its place.** Both methods are
unbiased with the same false-positive rate, but keeping UNKNOWN responders
raises power from 0.73 to 0.84.

**When either fails, the extension is biased and complete case is not.** The
six-stratum bias grows in proportion to the violation. With UNKNOWN responders
carrying all their expectancy, false positives triple (0.18). When the placebo
arm answers UNKNOWN 20 points more often, they reach 0.21. When the active arm
does, the bias runs the other way and power falls to 0.61.

**The bias follows a simple formula.** Under the AEB model,

  bias ≈ u · λ · μ_AEB · (2 p_cg − 1) − (1 − λ) · u_diff · μ_AEB / 2,

where λ is the share of expectancy UNKNOWN responders carry. The first term is
expectancy hidden in the UNKNOWN class: unstated beliefs still track the arm, so
the UNKNOWN cells differ by arm even though the method holds them fixed. The
second is the expectancy missing from whichever arm answers UNKNOWN more often.
Across the 20 cells of this grid the formula is within 0.08 of the simulated
bias, and it holds to within about 0.1 at a second design (n = 400,
p_cg = 0.8, u = 0.35). It is a heuristic approximation, checked at these two
designs only. `cgr_unknown_bias_approx()` computes it.

**Complete case is robust here for a specific reason.** In this model, whether
someone answers UNKNOWN is random within each arm, so dropping them leaves an
unbiased sample of directional responders. That would fail if answering UNKNOWN
depended on the belief itself or on the outcome, for example if people who
would have guessed wrong were more likely to say "I don't know". That case is
not tested here, and Section 6.2 of the summary paper shows that an
outcome-correlated UNKNOWN response damages the six-stratum estimator badly
(bias −2.02).

## Recommendation

- Report the six-stratum and complete-case estimates side by side. When they
  agree, the extra power of the six-stratum estimate can be used; when they
  disagree, the difference is evidence that A1 or A2 is failing.
- Check A1 directly: compare the UNKNOWN rate in each arm, which
  `cgr_guess_rates()` reports as `unknown_active` and `unknown_placebo`.
- A2 cannot be checked from the data. Use `cgr_unknown_bias_approx()` with the
  observed `u_diff` and a range of plausible λ to see how large a violation
  would have to be to change the conclusion.
