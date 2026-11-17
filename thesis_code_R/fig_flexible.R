# fig_flexible.R -- Figure 3 (Section 5.5): optimal signal under the budget K0, baseline calibration.
# (a) optimal posteriors and the threshold p*(alpha); (b) optimal accuracies q_G*, q_B*.
# Run:  Rscript fig_flexible.R   ->  figures/fig_flexible.(pdf|png)
source("flexible_learning.R")
P <- 0.5
a <- seq(0.02, 4, length.out = 400)
sol <- lapply(a, function(x) optimal_signal(P, pstar(x), K0))
pib <- sapply(sol, `[[`, "pi_b"); pig <- sapply(sol, `[[`, "pi_g")
qG <- sapply(sol, `[[`, "qG"); qB <- sapply(sol, `[[`, "qB")
a_p <- alpha_indifference(P)
save_fig("fig_flexible", 9, 2.9, function() {
  par(mfrow = c(1, 2), mar = c(3.4, 3.6, 1.8, 0.6), mgp = c(2.2, 0.7, 0), cex = 0.8, bty = "l")
  plot(a, pig, type = "l", col = CLEAR, lwd = 1.8, xlim = c(0, 4), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = expression("risk aversion " * alpha), ylab = "belief")
  title("(a) Optimal posterior beliefs", font.main = 1, cex.main = 1)
  lines(a, pib, col = ALARM, lwd = 1.8, lty = 2); lines(a, pstar(a), lwd = 1.4)
  abline(h = P, col = "gray", lty = 3)
  legend("bottomright", bty = "n", cex = 0.8, col = c(CLEAR, ALARM, "black"), lty = c(1, 2, 1), lwd = c(1.8, 1.8, 1.4),
         legend = c(expression(pi[g]^"*" * " (favourable report)"), expression(pi[b]^"*" * " (unfavourable report)"),
                    expression("threshold " * p^"*" * (alpha))))
  plot(a, qG, type = "l", col = CLEAR, lwd = 1.8, xlim = c(0, 4), ylim = c(0.25, 1.02), xaxs = "i", yaxs = "i",
       xlab = expression("risk aversion " * alpha), ylab = "accuracy")
  title("(b) Optimal signal accuracies", font.main = 1, cex.main = 1)
  lines(a, qB, col = ALARM, lwd = 1.8, lty = 2)
  abline(v = a_p, col = "gray", lty = 3); text(a_p + 0.05, 0.28, expression(alpha[p]), adj = 0, cex = 0.85)
  abline(h = 0.5, col = "gray", lty = 3)
  text(0.45, 0.36, "alarm type", cex = 0.8); text(0.45, 0.30, expression(q[G]^"*" > q[B]^"*"), cex = 0.8)
  text(3.2, 0.36, "all-clear type", cex = 0.8); text(3.2, 0.30, expression(q[G]^"*" < q[B]^"*"), cex = 0.8)
  legend("right", bty = "n", cex = 0.8, col = c(CLEAR, ALARM), lty = c(1, 2), lwd = 1.8,
         legend = c(expression(q[G]^"*" == P(g ~ "|" ~ G)), expression(q[B]^"*" == P(b ~ "|" ~ B))))
})
cat("Figure 3 written.\n")
