# Stress test of the binary CGRC under benign unblinding: a share of
# participants infer their arm from how they are doing, so their correct guess
# follows the effect instead of causing it (see docs/benign-unblinding.md).
# Writes docs/results/benign-stress.csv. Run from the repository root.

for (f in list.files("R", full.names = TRUE)) source(f)

res <- do.call(rbind, lapply(c(0, 0.25, 0.5, 0.75, 1), function(b) {
  op <- cgr_operating(n_trials = 500, n = 230, p_cg = 0.7, n_draws = 2000,
                      seed = 21, benign = b)
  set.seed(22)
  cgr <- sapply(seq_len(4), function(k) {
    z <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))[[k]]
    mean(replicate(500, { d <- sim_aeb(230, 0.7, z[1] == 1, z[2] == 1, benign = b)
                          mean(d$condition == d$guess) }))
  })
  cbind(benign = b, mean_observed_cgr = cgr, op)
}))
dir.create("docs/results", showWarnings = FALSE)
utils::write.csv(res, "docs/results/benign-stress.csv", row.names = FALSE)
