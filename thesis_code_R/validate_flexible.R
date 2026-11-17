# validate_flexible.R -- checks of the solver reported in Appendix B (B.1 accuracy, B.2 validation).
#  1. residuals of system (27) over a wide range of priors, thresholds and budgets;
#  2. comparison with a brute-force grid search over (q_G, q_B) for several priors, risk aversions and budgets;
#  3. shadow value (25) against a central finite difference of J*(K);
#  4. location of the maximum of the monetary VOI (Appendix B.2, "Location of the peak");
#  5. rise of the monetary VOI below alpha_p and fall above it on the grid of Figure 5 (Section 6.2);
#  6. ordering of the degree of asymmetry |q_G* - q_B*| across the budgets of Figure 5 (Section 6.3).
# Run:  Rscript validate_flexible.R   (takes a few minutes)
source("flexible_learning.R")
KAPPA <- K0 / log(2)

# 1. residuals ---------------------------------------------------------------------------------
P_GRID <- round(seq(0.1, 0.9, by = 0.1), 2)
PS_GRID <- round(seq(0.05, 0.95, by = 0.05), 2)
SHARES <- c(0.02, 0.1, 0.3, 0.6, 0.9)
eq <- 0; inf <- 0; n_ok <- 0; n_fail <- 0
for (p in P_GRID) for (ps in PS_GRID) for (sh in SHARES) {
  s <- tryCatch(optimal_signal(p, ps, sh * hbin(p)), error = function(e) NULL)
  if (is.null(s)) { n_fail <- n_fail + 1; next }
  n_ok <- n_ok + 1; eq <- max(eq, s$eq_residual); inf <- max(inf, s$info_rel_dev)
}
cat(sprintf("[1] residuals: %d of %d combinations solved (%d beyond double precision); max |F(z_b)-F(z_g)| = %.1e; max relative deviation of I from K = %.1e\n",
            n_ok, n_ok + n_fail, n_fail, eq, inf))

# 2. brute-force grid --------------------------------------------------------------------------
diffs <- c()
for (p in c(0.3, 0.5, 0.7)) for (alpha in c(0.3, 1, 2, 3)) for (sh in c(0.1, 0.2, 0.5)) {
  K <- sh * hbin(p)
  s <- optimal_signal(p, pstar(alpha), K)
  diffs <- c(diffs, W1(alpha, s$qG, s$qB, p) - grid_optimum(alpha, p, K, n = 801)[["W"]])
}
cat(sprintf("[2] grid comparison (801 x 801): %d cases; solver at least as good in %d; differences between %.1e and %.1e\n",
            length(diffs), sum(diffs >= 0), min(diffs), max(diffs)))

# 3. shadow value ------------------------------------------------------------------------------
EPS <- 1e-6; rel <- c()
for (p in c(0.3, 0.5, 0.7, 0.85)) for (ps in c(0.2, 0.4, 0.5, 0.6, 0.8)) for (sh in c(0.1, 0.4, 0.8)) {
  K <- sh * hbin(p)
  lam <- optimal_signal(p, ps, K)$lam
  fd <- (J_star(p, ps, K + EPS) - J_star(p, ps, K - EPS)) / (2 * EPS)
  rel <- c(rel, abs(fd - lam) / lam)
}
cat(sprintf("[3] shadow value: %d cases; largest relative deviation = %.1e\n", length(rel), max(rel)))

# 4. location of the maximum of the monetary VOI -----------------------------------------------
cat("[4] Location of the peak: VOI at alpha_p vs. maximum on 300 values of alpha in [0.02, 5]\n",
    "   (the three budgets of Figure 5 and the priors 0.3 and 0.7 of Section 6.3)\n", sep = "")
grid5 <- seq(0.02, 5, length.out = 300)
settings <- list(list(0.5, 0.5 * K0, "p = 0.5, K = 0.5 K0"), list(0.5, K0, "p = 0.5, K = K0"),
                 list(0.5, 2 * K0, "p = 0.5, K = 2 K0"), list(0.3, KAPPA * hbin(0.3), "p = 0.3, K = kappa0 h(0.3)"),
                 list(0.7, KAPPA * hbin(0.7), "p = 0.7, K = kappa0 h(0.7)"))
vf <- function(x, p, K) voi_flexible(x, p, K)$voi
worst <- -Inf; curves <- list()
for (st in settings) {
  p <- st[[1]]; K <- st[[2]]
  ap <- alpha_indifference(p); v_ap <- vf(ap, p, K)
  v <- sapply(grid5, vf, p = p, K = K); curves[[st[[3]]]] <- list(ap = ap, v = v)
  gmax <- max(v)
  fine <- c(seq(ap - 0.05, ap + 0.05, length.out = 2001), seq(0.02, 5, length.out = 1001)); fine <- fine[fine > 0]
  excess <- max(sapply(fine, vf, p = p, K = K)) - v_ap
  worst <- max(worst, excess)
  cat(sprintf("    %-28s alpha_p = %.3f  VOI(alpha_p) = %.4f  grid max = %.4f  largest excess on fine grids = %.1e\n",
              st[[3]], ap, v_ap, gmax, excess))
}
cat(sprintf("    no sampled alpha exceeds VOI(alpha_p) by more than %.1e\n", max(worst, 0)))

# 5. rise and fall of the monetary VOI around alpha_p (Section 6.2) ------------------------------
cat("[5] Monotonicity: VOI on the 300 values of alpha in [0.02, 5] of Figure 5;\n",
    "    consecutive grid points both below alpha_p must increase, both above alpha_p must decrease\n", sep = "")
for (lab in names(curves)) {
  ap <- curves[[lab]]$ap; d <- diff(curves[[lab]]$v)
  below <- grid5[-1] < ap                              # both points of the step lie below alpha_p
  above <- grid5[-length(grid5)] > ap                  # both points of the step lie above alpha_p
  cat(sprintf("    %-28s rises in %d of %d steps below alpha_p, falls in %d of %d steps above; smallest |step| = %.1e\n",
              lab, sum(d[below] > 0), sum(below), sum(d[above] < 0), sum(above), min(abs(d[below | above]))))
}

# 6. degree of the asymmetry across budgets (Section 6.3) --------------------------------------
cat("[6] Degree of asymmetry: |q_G* - q_B*| for K in {0.5, 1, 2} K0 at p = 1/2 on the grid of Figure 5\n")
ap <- alpha_indifference(0.5); n_ord <- 0; n_tot <- 0; gap <- Inf
for (x in grid5) {
  if (abs(x - ap) < 1e-9) next
  asym <- sapply(c(0.5, 1, 2), function(m) { s <- optimal_signal(0.5, pstar(x), m * K0); abs(s$qG - s$qB) })
  n_tot <- n_tot + 1
  n_ord <- n_ord + (asym[1] > asym[2] && asym[2] > asym[3])
  gap <- min(gap, asym[1] - asym[2], asym[2] - asym[3])
}
cat(sprintf("    |q_G* - q_B*| strictly largest with 0.5 K0 and smallest with 2 K0 at %d of %d values of alpha; smallest difference between adjacent budgets = %.1e\n",
            n_ord, n_tot, gap))
