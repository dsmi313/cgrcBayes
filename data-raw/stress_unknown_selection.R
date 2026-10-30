# Second stress test of the UNKNOWN-preserving estimator: the UNKNOWN response
# depends on the unstated belief (u_wrong) or on the outcome (u_outcome), with
# and without expectancy among UNKNOWN responders. Writes
# docs/results/unknown-stress-selection.csv. Run from the repository root.

for (f in list.files("R", full.names = TRUE)) source(f)

belief <- cgr_unknown_stress(u_diff = 0, aeb_unknown = c(0, 0.5),
                             u_wrong = c(0, 0.1, 0.2, 0.3), u_outcome = 0,
                             n_trials = 500, n = 300, p_cg = 0.7, u = 0.2,
                             n_draws = 2000, seed = 11)
outcome <- cgr_unknown_stress(u_diff = 0, aeb_unknown = c(0, 0.5),
                              u_wrong = 0, u_outcome = c(-0.2, -0.1, 0.1, 0.2),
                              n_trials = 500, n = 300, p_cg = 0.7, u = 0.2,
                              n_draws = 2000, seed = 12)
res <- rbind(cbind(test = "belief", belief), cbind(test = "outcome", outcome))
utils::write.csv(res, "docs/results/unknown-stress-selection.csv", row.names = FALSE)
