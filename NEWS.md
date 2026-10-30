# cgrc.bayes 0.1.0

First release.

* `cgrc()` and `cgr_conjugate()`: the four-stratum CGRC estimand with a
  conjugate Normal-Inverse-Gamma posterior. Reproduces Table 2 of Szigeti et
  al. (2023) from the public data.
* `cgr_jags()`: a JAGS backend whose prior matches the conjugate model, with
  `cgr_check_backends()` to confirm agreement and a Student-t sensitivity
  analysis.
* `cgr_kde()`, `cgr_kde_curve()` and `szigeti_panel()`: a source-faithful
  reimplementation of the original kernel density resampling procedure and
  figure.
* `sim_aeb()` and `cgr_operating()`: the activated-expectancy-bias model and
  repeated-simulation operating characteristics for checking a design.
* `cgr_rope()` and `cgrc_headline()`: a region of practical equivalence
  decomposition and a two-probability summary.
* `cgrc_unknown()`: a proposed extension that keeps "I don't know" guesses as a
  third response category. It is not part of the original method.
* `cgrc_app()`: a Shiny app for design and analysis.
