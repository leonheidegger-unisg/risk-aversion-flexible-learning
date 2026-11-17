# fig_benchmark.R -- Figure 1 (Section 3.5): fixed symmetric signal q = 0.75, baseline calibration.
# (a) threshold p*(alpha) and the posteriors 1-q, q; (b) VOI(alpha, q) with the three regimes shaded.
# alpha = 0 is evaluated from the analytical risk-neutral limit.
# Run:  Rscript fig_benchmark.R   ->  figures/fig_benchmark.(pdf|png)
source("flexible_learning.R")
Q <- 0.75; P <- 0.5
a_lo <- alpha_at(1 - Q); a_hi <- alpha_at(Q); a_p <- alpha_indifference(P)
a <- seq(0, 4, length.out = 801)
v <- voi(a, Q, Q, P)
regions <- list(c(0, a_lo), c(a_lo, a_hi), c(a_hi, 4))
cols <- c("#e8eef6", "#fdf0e3", "#eaf4ea")
labs <- c("I: always risky", "II: signal-dependent", "III: always safe")
shade <- function(ytop) for (k in 1:3) rect(regions[[k]][1], 0, regions[[k]][2], ytop, col = cols[k], border = NA)
save_fig("fig_benchmark", 9, 2.9, function() {
  par(mfrow = c(1, 2), mar = c(3.4, 3.6, 1.8, 0.6), mgp = c(2.2, 0.7, 0), cex = 0.8, bty = "l")
  plot(NA, xlim = c(0, 4), ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = expression("risk aversion " * alpha), ylab = "belief")
  title(expression("(a) Threshold " * p^"*" * (alpha) * " and posterior beliefs"), font.main = 1, cex.main = 1)
  shade(1); box(bty = "l")
  lines(a, pstar(a), lwd = 1.8)
  abline(h = c(Q, 1 - Q), col = "gray", lty = 2)
  text(3.95, Q + 0.04, expression(pi[g] == q * " = 0.75"), adj = 1, cex = 0.85)
  text(3.95, 1 - Q + 0.04, expression(pi[b] == 1 - q * " = 0.25"), adj = 1, cex = 0.85)
  plot(NA, xlim = c(0, 4), ylim = c(0, 0.235), xaxs = "i", yaxs = "i", xlab = expression("risk aversion " * alpha),
       ylab = expression(VOI(alpha, q)))
  title("(b) Value of the symmetric signal", font.main = 1, cex.main = 1)
  shade(0.235); box(bty = "l")
  lines(a, v, lwd = 1.8)
  abline(v = a_p, col = "gray", lty = 3); text(a_p + 0.05, 0.008, expression(alpha[p]), adj = 0, cex = 0.85)
  text(0.03, 0.222, labs[1], adj = 0, cex = 0.7)
  for (k in 2:3) text(mean(regions[[k]]), 0.222, labs[k], cex = 0.7)
})
cat(sprintf("Figure 1: alpha_low=%.3f, alpha_high=%.3f, alpha_p=%.3f, max VOI on grid=%.4f\n", a_lo, a_hi, a_p, max(v)))
