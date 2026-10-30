# Stress test of benign unblinding

The CGRC assumes unblinding is *malicious*: participants work out their arm
from side effects, and the belief then changes the outcome through expectancy.
If unblinding is *benign*, participants work out their arm from how they are
doing. A correct guess is then a consequence of the treatment effect, and
reweighting away from correct guessers removes part of a real effect. Szigeti et
al. (2023) flag this risk, and it cannot be tested from trial data (see U6 in
[open questions](open-questions.md)). This page puts numbers on it.

## Design

`sim_aeb()` takes a `benign` argument: the share of participants who infer their
arm from their outcome before expectancy, guessing active when it is above the
trial median. Their guess carries no expectancy. Everyone else guesses from side
effects as in the original model, correctly with probability `p_cg`.

`cgr_operating()` passes `benign` through, so the usual operating
characteristics can be read off for each share: n = 230, `p_cg` = 0.7, direct
effect 3, expectancy effect 7.7, 500 trials per scenario. The true effect is the
direct effect, because benign guessers carry no expectancy. Reproduce with
`Rscript data-raw/stress_benign.R`; full results are in
[`results/benign-stress.csv`](results/benign-stress.csv).

## Results

Bias of the CGR-adjusted estimate at CGR 0.50, and its power to detect the real
effect (`P(favourable) > 0.95`):

| Share of benign guessers | Observed CGR, real effect | Adjusted bias, real effect only | Adjusted bias, real effect + expectancy | 95% coverage, real effect only | Power, real effect only | Power, real effect + expectancy |
|---|---|---|---|---|---|---|
| 0 | 0.70 | 0.06 | −0.00 | 0.94 | 0.85 | 0.69 |
| 0.25 | 0.67 | −1.00 | −0.75 | 0.83 | 0.61 | 0.52 |
| 0.5 | 0.64 | −1.73 | −1.23 | 0.53 | 0.38 | 0.39 |
| 0.75 | 0.61 | −1.96 | −1.70 | 0.32 | 0.36 | 0.28 |
| 1 | 0.58 | −1.90 | −1.87 | 0.12 | 0.55 | 0.26 |

With no effect and no expectancy:

| Share of benign guessers | Observed CGR | Adjusted bias | 95% coverage | False positives |
|---|---|---|---|---|
| 0 | 0.70 | 0.03 | 0.95 | 0.05 |
| 0.25 | 0.65 | −0.52 | 0.85 | 0.00 |
| 0.5 | 0.60 | −0.68 | 0.71 | 0.00 |
| 0.75 | 0.55 | −0.47 | 0.82 | 0.01 |
| 1 | 0.50 | 0.00 | 0.96 | 0.05 |

## What this shows

**Adjustment removes a large part of a real effect.** If a quarter of
participants infer their arm from how they are doing, the adjusted estimate
loses about a third of a 3-point effect (bias −1.00). From half onward it loses
roughly 1.7–2.0 points, close to two-thirds, and power falls from 0.85 to
0.36–0.55.

**It can manufacture an effect in the wrong direction.** With no effect at all,
outcome-driven guessing sorts participants who are doing well into the
"guessed active" strata in both arms. Moving weight toward incorrect guessers
then pulls the active arm down and the placebo arm up, giving a spurious
negative estimate of about −0.5 to −0.7 and coverage as low as 0.71. The error
points toward harm, so false positives for a favourable effect do not rise.

**The data give no warning.** Benign guessing lowers the observed correct-guess
rate toward 0.50 (from 0.70 to 0.58 here), which looks like *better* blinding.
Both mechanisms also make participants who believe they received the active
treatment do better, so the stratum means cannot tell them apart.

**Malicious expectancy is still removed.** In the expectancy-only rows, the
adjusted bias stays within 0.17 at every share, so the method still does its
intended job for the participants whose unblinding is malicious. The damage
comes entirely from the benign guessers.

## Implications

- The CGRC-adjusted estimate recovers the direct effect only if unblinding is
  entirely malicious. Under any benign unblinding it can understate the effect
  substantially.
- A trial report should present the adjusted estimate as one end of a range,
  alongside the unadjusted contrast, not as a corrected effect.
- Evidence about *how* participants guessed is the only safeguard. In the
  microdosing trial, 55% named bodily sensations and 23% named mental or
  psychological benefits as their main cue (Szigeti et al. 2023). That supports
  a mostly malicious mechanism, but a benign share of even a quarter would
  matter. Collecting the basis for each guess would let an analysis bound
  `benign` instead of assuming it is zero.
