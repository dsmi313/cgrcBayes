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
unbiased sample of directional responders. The next section tests what happens
when answering UNKNOWN depends on the unstated belief or on the outcome.

## When the UNKNOWN response is selective

Two further departures make *who* answers UNKNOWN depend on something other
than arm, with UNKNOWN rates equal across arms on average:

- **Belief-dependent** (`u_wrong`): participants who would have guessed wrong
  answer UNKNOWN more often than those who would have guessed right, by the
  stated amount, with each arm's average UNKNOWN rate held at 0.2.
- **Outcome-dependent** (`u_outcome`): the probability of answering UNKNOWN
  changes by the stated amount per standard deviation of the outcome before
  expectancy. At ±0.2 per SD, someone one SD above average answers UNKNOWN
  with probability 0.4 (or 0) instead of 0.2.

Each is run with UNKNOWN responders carrying none or half of the expectancy
effect, at the same settings as above (500 trials per cell). Reproduce with
`Rscript data-raw/stress_unknown_selection.R`; full results are in
[`results/unknown-stress-selection.csv`](results/unknown-stress-selection.csv).
Each cell shows the six-stratum value first and the complete-case value second.

### Belief-dependent UNKNOWN responses

| Extra UNKNOWN rate, would-be-wrong guessers | Expectancy share | Bias, six-stratum (no effect / effect) | Bias, complete case (no effect / effect) | False positives, six / CC | Power, six / CC |
|---|---|---|---|---|---|
| 0 | 0 | 0.03 / −0.02 | 0.02 / 0.03 | 0.06 / 0.06 | 0.80 / 0.73 |
| 0 | 0.5 | 0.31 / 0.39 | 0.01 / 0.10 | 0.11 / 0.05 | 0.89 / 0.71 |
| 0.1 | 0 | 0.01 / −0.03 | 0.01 / −0.05 | 0.07 / 0.06 | 0.78 / 0.67 |
| 0.1 | 0.5 | 0.20 / 0.14 | 0.10 / 0.00 | 0.06 / 0.06 | 0.79 / 0.68 |
| 0.2 | 0 | 0.02 / −0.03 | 0.04 / −0.04 | 0.05 / 0.05 | 0.73 / 0.63 |
| 0.2 | 0.5 | −0.03 / −0.03 | −0.00 / 0.01 | 0.05 / 0.04 | 0.76 / 0.65 |
| 0.3 | 0 | 0.01 / 0.07 | 0.00 / 0.12 | 0.04 / 0.04 | 0.74 / 0.68 |
| 0.3 | 0.5 | −0.14 / −0.23 | 0.06 / −0.04 | 0.04 / 0.05 | 0.66 / 0.62 |

### Outcome-dependent UNKNOWN responses

| Change in UNKNOWN rate per outcome SD | Expectancy share | Bias, six-stratum (no effect / effect) | Bias, complete case (no effect / effect) | False positives, six / CC | Power, six / CC |
|---|---|---|---|---|---|
| −0.2 | 0 | −0.04 / 0.16 | −0.06 / −0.35 | 0.05 / 0.04 | 0.83 / 0.63 |
| −0.2 | 0.5 | 0.41 / 0.46 | 0.06 / −0.28 | 0.13 / 0.05 | 0.89 / 0.63 |
| −0.1 | 0 | 0.02 / 0.19 | 0.02 / −0.04 | 0.05 / 0.05 | 0.80 / 0.69 |
| −0.1 | 0.5 | 0.30 / 0.39 | −0.04 / −0.09 | 0.10 / 0.05 | 0.85 / 0.69 |
| 0.1 | 0 | 0.05 / −0.19 | 0.05 / −0.07 | 0.05 / 0.06 | 0.75 / 0.70 |
| 0.1 | 0.5 | 0.38 / 0.27 | 0.07 / 0.01 | 0.12 / 0.06 | 0.87 / 0.74 |
| 0.2 | 0 | −0.01 / −0.18 | −0.02 / −0.24 | 0.05 / 0.05 | 0.79 / 0.66 |
| 0.2 | 0.5 | 0.35 / 0.17 | 0.05 / −0.35 | 0.12 / 0.07 | 0.86 / 0.62 |

Across both tests, 95% coverage ranges from 0.91 to 0.97 for the six-stratum
estimator and from 0.93 to 0.98 for complete case, and the complete-case
false-positive rate never exceeds 0.07.

### What this adds

**Belief-dependent UNKNOWN responses do not break either method on their
own.** With no expectancy among UNKNOWN responders, both stay unbiased
(|bias| ≤ 0.12) at every level tested. Dropping UNKNOWN responders when the
dropout depends only on whether someone would have guessed right removes people
at random *within* each of the four binary strata, so the binary CGRC is
unaffected. When UNKNOWN responders do carry expectancy, belief-dependence
partly cancels the six-stratum bias, because would-be-wrong guessers lean
toward the arm they did not receive; at 0.3 it overshoots to −0.14 to −0.23.

**Outcome-dependent UNKNOWN responses bias both methods when there is a real
effect.** With a true effect of 3, complete case underestimates it by about
0.24–0.35 at ±0.2 per SD, and the six-stratum estimator is off by 0.16–0.19 in
the absence of expectancy among UNKNOWN responders. With no true effect,
neither method's false-positive rate rises above its level without selection:
about 0.05 for complete case, about 0.1–0.13 for the six-stratum estimator when
UNKNOWN responders carry expectancy. The mechanism is that a real effect raises
outcomes in the active arm, so outcome-dependent answering becomes
arm-dependent answering: at +0.2 per SD the active arm answers UNKNOWN at about
0.25 and the placebo arm at about 0.19.

**Taken together**, complete case is the safer default: across every departure
tested, its false-positive rate stays at about 0.05, and its worst bias is
about 0.35 points (roughly 12% of a 3-point effect), toward the null. The
six-stratum estimator is more powerful in every cell of these two tests, but
when UNKNOWN responders carry expectancy its false-positive rate can roughly
double (0.10–0.13 here, 0.18 at full expectancy in the first test). The size of
that cost is set mainly by `aeb_unknown`, which the data cannot reveal.

## Recommendation

- Report the six-stratum and complete-case estimates side by side. Treat the
  complete-case estimate as the primary analysis, since it kept a false-positive
  rate of about 0.05 under every departure tested. When the two agree, the
  extra power of the six-stratum estimate can be used; when they disagree, the
  difference is evidence that the extension's assumptions are failing.
- Check A1 directly: compare the UNKNOWN rate in each arm, which
  `cgr_guess_rates()` reports as `unknown_active` and `unknown_placebo`.
- A2 cannot be checked from the data. Use `cgr_unknown_bias_approx()` with the
  observed `u_diff` and a range of plausible λ to see how large a violation
  would have to be to change the conclusion.
