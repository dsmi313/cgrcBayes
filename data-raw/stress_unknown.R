# Stress test of the UNKNOWN-preserving estimator when its simulation
# assumptions fail (see docs/unknown-stress-test.md). Writes
# docs/results/unknown-stress.csv. Run from the repository root; takes about
# two minutes.

for (f in list.files("R", full.names = TRUE)) source(f)

res <- cgr_unknown_stress(u_diff = c(-0.2, -0.1, 0, 0.1, 0.2),
                          aeb_unknown = c(0, 0.25, 0.5, 1),
                          n_trials = 500, n = 300, p_cg = 0.7, u = 0.2,
                          mu_dte = 3, mu_aeb = 7.7, n_draws = 2000, seed = 1)

# Closed-form approximation to the six-stratum bias when only A2 fails
# (u_diff = 0): u * aeb_unknown * mu_aeb * (2 * p_cg - 1).
res$approx_bias <- ifelse(res$method == "six-stratum" & res$u_diff == 0,
                          0.2 * res$aeb_unknown * 7.7 * (2 * 0.7 - 1), NA_real_)

dir.create("docs/results", showWarnings = FALSE)
utils::write.csv(res, "docs/results/unknown-stress.csv", row.names = FALSE)
