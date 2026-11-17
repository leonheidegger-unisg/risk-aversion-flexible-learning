# fig_flex_vs_fixed.R -- Figure 4 (Section 6.2): VOI of the optimal flexible signal (budget K0)
# versus the fixed symmetric signal q = 0.75, which uses exactly K0 nats at p = 1/2.
# Run:  Rscript fig_flex_vs_fixed.R   ->  figures/fig_flex_vs_fixed.(pdf|png)
source("flexible_learning.R")
a <- seq(0.02, 4, length.out = 400)
flex <- sapply(a, function(x) voi_flexible(x, 0.5, K0)$voi)
sym <- voi(a, 0.75, 0.75, 0.5)
ap <- alpha_indifference(0.5)
save_fig("fig_flex_vs_fixed", 5.8, 3.0, function() {
  par(mar = c(3.4, 3.6, 0.8, 0.6), mgp = c(2.2, 0.7, 0), cex = 0.8, bty = "l")
  plot(NA, xlim = c(0, 4), ylim = c(0, 0.23), xaxs = "i", yaxs = "i", xlab = expression("risk aversion " * alpha),
       ylab = "value of information")
  polygon(c(a, rev(a)), c(flex, rev(sym)), col = adjustcolor(FLEX, alpha.f = 0.12), border = NA)
  lines(a, flex, col = FLEX, lwd = 2); lines(a, sym, lwd = 1.6, lty = 2)
  abline(v = ap, col = "gray", lty = 3); text(ap + 0.05, 0.008, expression(alpha[p]), adj = 0, cex = 0.85)
  legend("topright", bty = "n", cex = 0.85, col = c(FLEX, "black"), lty = c(1, 2), lwd = c(2, 1.6),
         legend = c(expression("optimal flexible signal (budget " * K[0] * ")"), expression("fixed symmetric signal " * q == 0.75)))
})
cat("Figure 4 written.\n")
