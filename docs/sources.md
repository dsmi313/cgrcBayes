# Sources

The literature and code this project reproduces and cites. The journal articles
are copyrighted and are not included in the repository; cite them by DOI. The
labels S1–S4 are stable identifiers used in [open questions](open-questions.md),
not a ranking.

| Role | Source |
|---|---|
| Primary methodological source | S1, Szigeti et al. 2023, *Scientific Reports* |
| Empirical dataset | Szigeti et al. 2021, *eLife* |
| Supporting literature | S3, Szigeti & Heifets 2024, *Biol Psychiatry CNNI* |
| Related expectancy literature | S2, Szigeti et al. 2024, *Psychological Medicine* |
| Original analysis code | S4, `szb37/CorrectGuessRateCurve` |

## S1. The CGRC method paper (primary methodological source)

Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D. **"The difference between
'placebo group' and 'placebo control': a case study in psychedelic
microdosing."** *Scientific Reports* 13:12107 (2023).
<https://doi.org/10.1038/s41598-023-34938-7>

The source of everything reproduced here: the activated expectancy bias (AEB)
model (Fig 1, Eq 1–4), the Correct Guess Rate Curve (Fig 2), the simulation
operating characteristics (Table 1), the microdose re-analysis (Table 2), and the
empirical CGR curves (Fig 4). Facts checked against it:

- **Analysis rule (U2).** "We only used data from the first week... each
  datapoint is independent... In the current analysis n = 233 datapoints were
  included." This confirms the `tp == "w1s1"` filter; the public file gives 232.
- **The 0.72 CGR (U3).** The Fig 4 caption: "vertical green dashed line
  corresponds to the trial's original CGR (= 0.72)". The source code (S4) shows
  this is a hardcoded constant, not the data's CGR of 0.647. That it comes from
  the placebo-arm rate (0.7234) is a hypothesis.
- **`noise = "all"` (U4).** Table 1 rates are 0.05 / 0.86 / 0.78 / 0.99, with a
  DTE of Hedges' g = 0.4, in the stated regime "n ~ 200, CGR ~ 0.7, effect ~ 0.4
  Hedges' g". `noise = "all"` reproduces these most closely, though not exactly.
  The exact Eq. 4 noise scope could not be confirmed from the available source.
- **The malicious-unblinding argument (U6).** 55% named body or perceptual cues
  (muscle tension 58%, stomach discomfort 27%) against 23% naming mental or
  psychological benefits; placebo-microdose PANAS differences of 2.1 and 0.8
  against about 10 and 6 of natural within-subject variability.
- **Table 2.** PANAS 3.2 → 1.1, Mood VAS 6.4 → 2.7, Energy VAS 11.5 → 6.8
  (g 0.58 → 0.34), with only Energy VAS surviving CGR adjustment. All reproduced
  from the public data.

## The empirical dataset

Szigeti B, et al. **"Self-blinding citizen science to explore psychedelic
microdosing."** *eLife* 10:e62878 (2021).
<https://doi.org/10.7554/eLife.62878>

The source study for the public data re-analysed here. The package pins a
checksum-verified copy as `data/pacutes.csv` (mirrored at
`szb37/mcrds_public`); the URL, access date and SHA-256 are in
`data/PROVENANCE.txt`.

## S2. Related expectancy literature (context only, different trial)

Szigeti B, Weiss B, Rosas FE, Erritzoe D, Nutt D, Carhart-Harris R. **"Assessing
expectancy and suggestibility in a trial of escitalopram v. psilocybin for
depression."** *Psychological Medicine* 54:1717–1724 (2024).
<https://doi.org/10.1017/S0033291723003653>

A different trial, and not a source of the CGRC method or the data used here.
Cited because it develops the same theme of expectancy and blinding failure (for
example, it notes 94% correct guessing in another psilocybin trial).

## S3. Supporting literature (context, not a methodological source)

Szigeti B, Heifets BD. **"Expectancy Effects in Psychedelic Trials."**
*Biological Psychiatry: Cognitive Neuroscience and Neuroimaging* 9:512–521
(2024). <https://doi.org/10.1016/j.bpsc.2024.02.004>

A review by the same author. It supports the interpretation of expectancy
effects and blinding failure but does not provide the CGRC method. Two facts are
used:

- The microdose correct-guess rate is "only ~65% to 70%", against about 50% in a
  truly blind trial. This is consistent with the public data's 0.647, not with
  the 0.72 line in the 2023 Fig 4 (U3).
- Microdose effects fall below the 0.5-SMD minimal important difference, "too
  small to be noticeable" (U6).

## S4. The original analysis code

Szigeti B. **CorrectGuessRateCurve**, analysis code for S1. MIT licence,
DrugNerdsLab, 2022. <https://github.com/szb37/CorrectGuessRateCurve>

Named in S1's data-availability statement. Reading it settled two questions and
qualified two descriptions of the method:

- **0.72 is hardcoded (U3).** `src/config.py` sets `trial_cgrs = {'sbmd': 0.72}`,
  drawn as the Figure 4 reference line. The same code computes the trial CGR from
  the data as `(n_plpl + n_acac) / n` = 0.647 (`src/cgrc/core.py`, around line
  181), but the plotted line does not use it.
- **The estimand is confirmed (U4).** `get_strata_ratio` and
  `get_strata_sample_sizes` in `src/cgrc/core.py` form `r = PLPL / (PLPL + ACAC)`
  and `s = ACPL / (ACPL + PLAC)`, identical to this implementation. The default
  `strata_sampling = 'all_prop'` controls pseudo-sample allocation; it says
  nothing about the scope of the Eq. 4 noise.
- **Resample count.** The Figure 4 configuration `cgrC_low` in `config.py`
  specifies `n_cgrc_trials = 32` over `np.linspace(0, 1, 13)` (the options are 32,
  64 and 96, never 100). Two caveats: the repository's README says it reproduces
  the preprint, which had a different title, so the configuration may predate the
  published version; and the Figure 4 block in `run.py` sits inside `if False:`.
  This is the repository's stated configuration, not proof of what generated the
  published figure. If the true count is 32, the Monte Carlo error is about 1.8
  times that of 100 resamples.
- **Rounding.** `get_strata_ratio` applies `round(x, 2)` to every stratum
  proportion, so two-decimal proportions are part of the source-faithful
  allocation; integer rounding of KDE scores is also required.
