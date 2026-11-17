# reproduce_numbers.R -- prints the numerical values quoted in Sections 3-6 of the thesis, section by section
# (the validation statistics of Appendix B are printed by validate_flexible.R and validate_direction.R).
# Run:  Rscript reproduce_numbers.R
source("flexible_learning.R")
fine <- seq(0, 10, length.out = 100001)
window <- function(qG, qB, p = 0.5) { po <- posteriors(qG, qB, p); c(alpha_at(po$pi_b), alpha_at(po$pi_g)) }

cat("=== Section 3: fixed symmetric signal ===\n")
cat(sprintf("risk-neutral threshold (s-L)/(H-L) = %.4f\n", (S - L) / (H - L)))
w <- window(0.75, 0.75); ap <- alpha_indifference(0.5)
cat(sprintf("q = 0.75: valuable for alpha in (%.3f, %.3f); alpha_p = %.3f; VOI(alpha_p) = %.4f\n", w[1], w[2], ap, voi(ap, .75, .75)))
cat(sprintf("q = 0.95: unfavourable posterior %.2f < risk-neutral threshold %.4f, so Regime I is empty\n", 0.05, (S - L) / (H - L)))

cat("\n=== Section 4: three fixed signals (Table 1) ===\n")
for (sg in list(list("symmetric", .75, .75), list("alarm", .90, .55), list("all-clear", .55, .90))) {
  qG <- sg[[2]]; qB <- sg[[3]]; po <- posteriors(qG, qB); w <- window(qG, qB); v <- voi(fine, qG, qB)
  cat(sprintf("%-10s P(b|G)=%.2f P(g|B)=%.2f pi_b=%.3f pi_g=%.3f window=(%.3f, %.3f) max VOI=%.4f at alpha=%.3f I=%.4f\n",
              sg[[1]], 1 - qG, 1 - qB, po$pi_b, po$pi_g, w[1], w[2], max(v), fine[which.max(v)], mutual_info(qG, qB)))
}
d_al <- function(a) voi(a, .90, .55) - voi(a, .75, .75)
d_ac <- function(a) voi(a, .55, .90) - voi(a, .75, .75)
cat(sprintf("alarm beats symmetric below alpha = %.3f; all-clear beats symmetric above alpha = %.3f\n",
            bisect(d_al, 0.5, 1.3, 100), bisect(d_ac, 1.4, 2.0, 100)))
cat(sprintf("at alpha = 0.25: VOI symmetric = %.4f, alarm = %.4f, all-clear = %.4f\n",
            voi(.25, .75, .75), voi(.25, .9, .55), voi(.25, .55, .9)))
cat(sprintf("example (q_G, q_B) = (0.4, 0.9): P(g|G)/P(g|B) = %.1f\n", 0.4 / 0.1))

cat("\n=== Section 5: flexible learning ===\n")
cat(sprintf("K0 = %.4f nats = %.1f%% of ln 2\n", K0, 100 * K0 / log(2)))
for (a in c(0.25, 3)) {
  s <- optimal_signal(0.5, pstar(a), K0)
  cat(sprintf(paste0("alpha = %g: p* = %.3f; (q_G*, q_B*) = (%.3f, %.3f); P(g) = %.3f, P(b) = %.3f; posteriors (%.3f, %.3f); ",
                     "missed-opportunity rate %.3f, false-green-light rate %.3f; distances to threshold: below %.2f, above %.2f\n"),
              a, pstar(a), s$qG, s$qB, s$P_g, 1 - s$P_g, s$pi_b, s$pi_g, 1 - s$qG, 1 - s$qB, pstar(a) - s$pi_b, s$pi_g - pstar(a)))
}
mn <- function(a) { s <- optimal_signal(0.5, pstar(a), K0); min(s$qG, s$qB) - 0.5 }
cat(sprintf("both accuracies at least 1/2 exactly for alpha in (%.3f, %.3f)\n", bisect(mn, 0.3, 1.3, 100), bisect(mn, 1.4, 3.0, 100)))

cat("\n=== Section 6: numerical results ===\n")
for (a in c(0.25, ap, 3))
  cat(sprintf("alpha = %.3f: VOI flexible = %.4f, VOI symmetric = %.4f\n", a, voi_flexible(a, 0.5, K0)$voi, voi(a, .75, .75)))
cat(sprintf("information of alarm / all-clear = %.4f < K0 = %.4f\n", mutual_info(.9, .55), K0))
for (m in c(0.5, 1, 2)) {
  s <- optimal_signal(0.5, pstar(1), m * K0); t <- optimal_signal(0.5, pstar(0.25), m * K0)
  cat(sprintf("K = %g K0: posteriors at alpha = 1: (%.3f, %.3f); q_G* - q_B* at alpha = 0.25: %.2f\n",
              m, s$pi_b, s$pi_g, t$qG - t$qB))
}
for (p in c(0.3, 0.7)) cat(sprintf("p = %g: alpha_p = %.3f\n", p, alpha_indifference(p)))

cat("\n=== Appendix B: limit of double precision (baseline, K0) ===\n")
ok <- function(a) tryCatch({ optimal_signal(0.5, pstar(a), K0); -1 }, error = function(e) 1)
a_lim <- bisect(ok, 5, 10, 60)
s <- optimal_signal(0.5, pstar(a_lim - 1e-3), K0)
cat(sprintf("the solver stops above alpha = %.2f (there 1 - pi_g* = %.1e, p* = %.3f)\n", a_lim, 1 - s$pi_g, pstar(a_lim)))
for (m in c(0.5, 1, 2)) {                     # largest risk aversion plotted (alpha = 5), all budgets of Figure 5
  s <- optimal_signal(0.5, pstar(5), m * K0)
  cat(sprintf("K = %g K0, alpha = 5: distance of the closest optimal posterior to certainty = %.1e\n", m, min(s$pi_b, 1 - s$pi_g)))
}
