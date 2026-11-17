# fig_asymmetric.R -- Figure 2 (Section 4.3): symmetric, alarm and all-clear signals, baseline calibration.
# (a) threshold p*(alpha) and the posteriors of the alarm and all-clear signals; (b) their VOI.
# Run:  Rscript fig_asymmetric.R   ->  figures/fig_asymmetric.(pdf|png)
source("flexible_learning.R")
P <- 0.5
sig <- list(list(q = c(0.75, 0.75), col = "black", lty = 1, lab = "symmetric (0.75, 0.75)"),
            list(q = c(0.90, 0.55), col = ALARM,   lty = 2, lab = "alarm (0.90, 0.55)"),
            list(q = c(0.55, 0.90), col = CLEAR,   lty = 4, lab = "all-clear (0.55, 0.90)"))
a <- seq(0, 4, length.out = 1601)
a_p <- alpha_indifference(P)
save_fig("fig_asymmetric", 9, 2.9, function() {
  par(mfrow = c(1, 2), mar = c(3.4, 3.6, 1.8, 0.6), mgp = c(2.2, 0.7, 0), cex = 0.8, bty = "l")
  plot(a, pstar(a), type = "l", lwd = 1.8, xlim = c(0, 4), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = expression("risk aversion " * alpha), ylab = "belief")
  title(expression("(a) Threshold " * p^"*" * (alpha) * " and posterior beliefs"), font.main = 1, cex.main = 1)
  for (s in sig[2:3]) { po <- posteriors(s$q[1], s$q[2], P); abline(h = c(po$pi_g, po$pi_b), col = s$col, lty = s$lty) }
  text(3.95, 0.846 + 0.04, expression("all-clear " * pi[g] == 0.846), adj = 1, cex = 0.8, col = CLEAR)
  text(3.95, 0.667 - 0.05, expression("alarm " * pi[g] == 0.667), adj = 1, cex = 0.8, col = ALARM)
  text(3.95, 0.333 + 0.04, expression("all-clear " * pi[b] == 0.333), adj = 1, cex = 0.8, col = CLEAR)
  text(3.95, 0.154 - 0.05, expression("alarm " * pi[b] == 0.154), adj = 1, cex = 0.8, col = ALARM)
  plot(NA, xlim = c(0, 4), ylim = c(0, 0.23), xaxs = "i", yaxs = "i", xlab = expression("risk aversion " * alpha),
       ylab = "value of information")
  title("(b) Value of the three signals", font.main = 1, cex.main = 1)
  for (s in sig) lines(a, voi(a, s$q[1], s$q[2], P), col = s$col, lty = s$lty, lwd = 1.8)
  abline(v = a_p, col = "gray", lty = 3); text(a_p + 0.05, 0.008, expression(alpha[p]), adj = 0, cex = 0.85)
  legend("topright", bty = "n", cex = 0.8, col = sapply(sig, `[[`, "col"), lty = sapply(sig, `[[`, "lty"), lwd = 1.8,
         legend = sapply(sig, function(s) sprintf("%s, I = %.3f", s$lab, mutual_info(s$q[1], s$q[2], P))))
})
cat("Figure 2 written.\n")
