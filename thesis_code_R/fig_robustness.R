# fig_robustness.R -- Figure 5 (Section 6.3): VOI of the optimal flexible signal for the budgets
# K in {0.5, 1, 2} K0, baseline calibration (H, L, s, p) = (4, 1, 1.5, 0.5).
# (The peak check for the priors p = 0.3 and 0.7 is in validate_flexible.R, Appendix B.2.)
# Run:  Rscript fig_robustness.R   ->  figures/fig_robustness.(pdf|png)
source("flexible_learning.R")
a <- seq(0.02, 5, length.out = 300)
curve_K <- function(K) sapply(a, function(x) voi_flexible(x, 0.5, K)$voi)
mults <- c(0.5, 1, 2); ltys <- c(3, 1, 2)
curves <- lapply(mults, function(m) curve_K(m * K0))
ap <- alpha_indifference(0.5)
save_fig("fig_robustness", 6.2, 2.9, function() {
  par(mar = c(3.4, 3.6, 0.8, 0.6), mgp = c(2.2, 0.7, 0), cex = 0.8, bty = "l")
  plot(NA, xlim = c(0, 5), ylim = c(0, max(unlist(curves)) * 1.05), xaxs = "i", yaxs = "i",
       xlab = expression("risk aversion " * alpha), ylab = "value of information")
  for (k in 1:3) lines(a, curves[[k]], col = FLEX, lty = ltys[k], lwd = 1.8)
  abline(v = ap, col = "gray", lty = 3); text(ap + 0.05, 0.01, expression(alpha[p]), adj = 0, cex = 0.85)
  legend("topright", bty = "n", cex = 0.85, col = FLEX, lty = ltys, lwd = 1.8,
         legend = c(expression(K == 0.5 * K[0]), expression(K == K[0]), expression(K == 2 * K[0])))
})
cat("Figure 5 written.\n")
