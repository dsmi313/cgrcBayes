# cgrc.bayes

[![R-CMD-check](https://github.com/dsmi313/cgrcBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dsmi313/cgrcBayes/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Bayesian implementation of the **Correct Guess Rate Curve** (CGRC) of
[Szigeti et al. (2023)](https://doi.org/10.1038/s41598-023-34938-7), for
placebo-controlled trials where participants can tell which arm they are in.

When participants work out their allocation, expectancy can inflate the
treatment contrast. The CGRC reweights the participants a trial already has to
show what the contrast would be at a different correct-guess rate, usually 0.50
(guessing at chance). This package:

- **reproduces** the CGRC estimand exactly, and reports its uncertainty as a
  **posterior distribution** instead of an average of resampled p-values;
- **checks a design before analysis** by simulating the trial under the
  activated-expectancy-bias model, as the original paper recommends;
- **keeps "I don't know" guesses** with a proposed six-stratum extension that
  reduces exactly to the original when there are none;
- ships an **interactive Shiny app** so trialists can use all of this without
  writing code.

The adjusted estimate is a sensitivity analysis for unblinding, not a causal
decomposition: shrinkage under adjustment does not measure how much of an
effect was expectancy. If some participants work out their arm from how they are
doing rather than from side effects, adjustment also removes part of a real
effect; the [benign-unblinding stress test](docs/benign-unblinding.md) puts
numbers on this.

## Install

```r
# install.packages("remotes")
remotes::install_github("dsmi313/cgrcBayes")
```

Checked on the current R release on every push, and tested locally on R 4.3. The app also needs `shiny`, `ggplot2` and `scales`. The
optional JAGS backend needs `rjags` and a system
[JAGS](https://mcmc-jags.sourceforge.io) installation.

## Quick start

```r
library(cgrc.bayes)

set.seed(1)
trial <- sim_aeb(n = 230, p_cg = 0.7, dte_on = TRUE, aeb_on = TRUE)  # simulated trial

fit <- cgrc(trial)       # posterior curve, summarised at the observed CGR and at 0.50
fit
plot(fit)                # effect with 95% credible band, and P(favourable)

cgrc_headline(trial)     # two plain probabilities, before and after adjustment
```

### Your own trial

`cgrc()` needs one row per participant with three columns: `condition`
(`"AC"`/`"PL"`), `guess` (`"AC"`/`"PL"`) and `value`.

```r
d <- read.csv("my_trial.csv")
my_trial <- data.frame(
  condition = ifelse(d$arm   == "drug", "AC", "PL"),
  guess     = ifelse(d$guess == "drug", "AC", "PL"),
  value     = d$outcome
)
cgrc(my_trial)                    # use direction = -1 if lower scores are better
```

### Check the design first

Szigeti et al. advise that researchers "should first run simulations to
determine whether CGR produces acceptable error rates for the parameters of
their data". `cgr_operating()` does that at your sample size and guess rate:

```r
cgr_operating(n = 120, p_cg = 0.85, n_trials = 500)
```

It reports bias, RMSE, 95% coverage, power and false-positive rates, and the
share of simulated trials in which a stratum comes up empty. At a high guess
rate with a small sample, the wrong-guess strata can empty out and the
adjustment is undefined, not merely imprecise.

### Interactive app

```r
cgrc_app()
```

**Design** reads a precomputed simulation grid and shows, for your planned
design, power, the false-attribution trade-off with and without adjustment, and
expected stratum sizes. **Analyse** takes a CSV of your trial, audits every row,
runs the adjustment (showing the complete-case and six-stratum estimates side
by side when there are "I don't know" answers), and exports a summary CSV, the curve, a plot and an HTML
report.

### "I don't know" guesses

```r
trial_u <- sim_aeb_unknown(n = 300, u = 0.2, dte_on = TRUE)  # guesses include "UNKNOWN"
cgrc_unknown(trial_u)
```

UNKNOWN answers are kept as a third category instead of being dropped or scored
as wrong. The directional correct-guess rate is reweighted while the UNKNOWN
rate is held fixed. This extension is **not part of the original CGRC**. It is
unbiased and more powerful than dropping UNKNOWN responders when UNKNOWN rates
are equal across arms and UNKNOWN responders carry no expectancy, and biased when
either fails. The [stress test](docs/unknown-stress-test.md) quantifies this;
`cgr_unknown_stress()` and `cgr_unknown_bias_approx()` let you check your own
design. Report it alongside the complete-case analysis, which kept a
false-positive rate of about 5% under every departure tested.

## What the reproduction shows

On the public microdosing data (Szigeti et al. 2021; week 1, n = 232):

| Outcome | Unadjusted | Adjusted (CGR 0.50) | Published |
|---|---|---|---|
| PANAS | 3.16 | 1.08 (95% CrI −1.58 to 3.79) | 3.2 → 1.1 |
| Mood VAS | 6.34 | 2.52 (95% CrI −3.16 to 8.31) | 6.4 → 2.7 |
| Energy VAS | 11.38 | 7.11 (95% CrI 2.18 to 12.07) | 11.5 → 6.8 |

As published, only Energy VAS survives adjustment. At the observed CGR the
reweighted estimate equals the raw arm difference exactly, and the conjugate and
JAGS backends agree within Monte Carlo error across the whole curve.

Over 500 simulated trials per scenario (n = 230, CGR 0.7), the adjusted
estimator has bias of at most 0.03 and 95% credible-interval coverage of
0.936–0.964. When the apparent effect is pure expectancy, it flags a favourable
effect in 5% of trials against 80% for an unadjusted t-test. The cost is power:
with a real effect plus expectancy, it flags the effect in 69% of trials.

## Documentation

| | |
|---|---|
| [Summary paper](docs/summary.md) | The method, the "I don't know" extension, and three external trials |
| [Detailed overview](docs/detailed_overview.md) | Full reproduction, worked example, derivations, simulations and sensitivity analyses |
| [Plain-language guide](docs/plain-language.md) | The method without equations |
| [Assumptions](docs/assumptions.md) | Every assumption, and what follows if it fails |
| [Benign unblinding](docs/benign-unblinding.md) | How much real effect adjustment removes when participants infer their arm from how they are doing |
| [Stress test of the extension](docs/unknown-stress-test.md) | How the "I don't know" extension behaves when its assumptions fail |
| [Open questions](docs/open-questions.md) | Discrepancies with the original paper, ranked |
| [Sources](docs/sources.md) | Papers and code, and which facts come from each |

Every exported function has a help page: start with `?cgrc`.

## Repository layout

| Path | Contents |
|---|---|
| `R/` | The package: estimand, conjugate and JAGS posteriors, KDE reimplementation, simulation, plotting, ROPE, the UNKNOWN extension and app helpers |
| `inst/app/` | The Shiny app |
| `inst/extdata/` | Precomputed design grids and three external-trial datasets with provenance |
| `tests/testthat/` | Test suite |
| `data/` | The pinned public dataset, with its SHA-256 in `PROVENANCE.txt` |
| `data-raw/` | Scripts that fetch the data and build the design grids |
| `docs/` | The writeups, their figures, and the R Markdown source of the detailed overview |

## Reproducing

```sh
Rscript data-raw/download_data.R   # re-fetch data/pacutes.csv and verify its checksum
Rscript tests/testthat.R           # run the test suite from source
Rscript docs/render.R              # regenerate the detailed overview (needs JAGS)
R CMD build . && R CMD check cgrc.bayes_0.1.0.tar.gz
```

## Citation

Run `citation("cgrc.bayes")`, or use the "Cite this repository" button on
GitHub. Please also cite the original method:

> Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D (2023). The difference
> between 'placebo group' and 'placebo control': a case study in psychedelic
> microdosing. *Scientific Reports* 13:12107.
> <https://doi.org/10.1038/s41598-023-34938-7>

This project reproduces and re-analyses published work. It is independent of,
and not endorsed by, the original authors.

## License

MIT © Dave Smith
