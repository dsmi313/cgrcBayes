#' cgrc.bayes: Bayesian Implementation of the Correct Guess Rate Curve
#'
#' When participants in a placebo-controlled trial can tell which arm they are
#' in, expectancy can inflate the treatment contrast. The Correct Guess Rate
#' Curve (CGRC) of Szigeti et al. (2023) reweights the participants a trial
#' already has to show what the contrast would look like at a different
#' correct-guess rate, usually 0.50 (guessing at chance). This package
#' reproduces that estimand and reports its uncertainty as a posterior
#' distribution.
#'
#' @section Main functions:
#' * [cgrc()] adjusts a trial; [cgrc_headline()] summarises it in two
#'   probabilities.
#' * [cgr_operating()] checks, by simulation, whether CGR adjustment is reliable
#'   for a trial design.
#' * [cgrc_unknown()] analyses trials whose guess question allowed
#'   "I don't know".
#' * [cgr_rope()] splits the adjusted posterior into harm, negligible and
#'   benefit.
#' * [cgrc_app()] opens the interactive design and analysis app.
#'
#' @references Szigeti B, Nutt D, Carhart-Harris R, Erritzoe D (2023). The
#'   difference between 'placebo group' and 'placebo control': a case study in
#'   psychedelic microdosing. *Scientific Reports* 13:12107.
#'   \doi{10.1038/s41598-023-34938-7}
"_PACKAGE"

# Column names used inside ggplot2::aes() in cgr_plot() and cgr_unknown_plot().
utils::globalVariables(c("cgr", "est", "hi", "lo", "method", "y"))
