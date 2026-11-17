# validate_direction.R -- direction of the asymmetry, Appendix B.2.
# Checks (i) the criterion q_G* > q_B* <=> l(pi_g*) + l(pi_b*) < 2 l(p) (Section 5.5; (A.5) for p = 1/2)
# and (ii) the sign rules reported in Appendix B.2.
# Grid: p and p* each in {0.10, 0.15, ..., 0.90}, five budgets K/h(p).
# In the budget-dependent region (p* strictly between p and 1/2) the critical budget, at which the optimal
# signal is symmetric (the zero of G below), is computed and the sign at each sampled budget is compared
# with the side of that critical budget. The script also reports the range of the critical budgets.
# Run:  Rscript validate_direction.R
source("flexible_learning.R")
GRID <- round(seq(0.10, 0.90, by = 0.05), 2)
SHARES <- c(0.02, 0.1, 0.3, 0.6, 0.9)
TOL <- 1e-9

sign_rule <- function(p, ps) {
  # sign of q_G* - q_B* for every budget (Corollary 2 at p = 1/2; Appendix B.2 otherwise), or NA in the budget-dependent region
  if (p == 0.5) return(sign(0.5 - ps))
  if (p > 0.5) { if (ps >= p) return(-1); if (ps <= 0.5) return(1); return(NA) }   # NA: 1/2 < p* < p
  if (ps <= p) return(1); if (ps >= 0.5) return(-1); NA                           # NA: p < p* < 1/2
}

critical_share <- function(p, ps) {
  # K/h(p) at which the optimal signal is symmetric, for the budget-dependent region: zero of
  # G(d) = ln(1 + e^(m + d)) - ln(1 + e^(m - d)) - 2 p* d, m = l(p), where the posteriors are symmetric in log-odds around m.
  flip <- p < 0.5                                    # relabelling of the states (Appendix A.2)
  if (flip) { p <- 1 - p; ps <- 1 - ps }
  m <- logit(p)
  G <- function(d) softplus(m + d) - softplus(m - d) - 2 * ps * d
  d0 <- m - partner(m, ps)                           # limit K -> 0
  d_star <- bisect(G, d0 + 1e-9, 1e4, n = 200)
  zb <- m - d_star
  list(share = info_pair(p, zb, partner(zb, ps))[["I"]] / hbin(p), flip = flip)
}

n <- 0; skipped <- 0; crit_ok <- 0; rule_n <- 0; rule_ok <- 0; open_n <- 0; open_ok <- 0
c_all <- c()                                          # critical budgets K/h(p) in the budget-dependent region
for (p in GRID) for (ps in GRID) {
  rule <- sign_rule(p, ps)
  if (is.na(rule)) { cs <- critical_share(p, ps); c_all <- c(c_all, cs$share) }
  for (sh in SHARES) {
    s <- tryCatch(optimal_signal(p, ps, sh * hbin(p)), error = function(e) NULL)
    if (is.null(s)) { skipped <- skipped + 1; next }    # optimal posterior beyond double precision
    n <- n + 1
    d <- s$qG - s$qB
    sgn <- if (abs(d) < TOL) 0 else sign(d)
    crit <- 2 * logit(p) - logit(s$pi_g) - logit(s$pi_b)
    crit_ok <- crit_ok + (sgn == (if (abs(crit) < TOL) 0 else sign(crit)))
    if (!is.na(rule)) { rule_n <- rule_n + 1; rule_ok <- rule_ok + (sgn == rule) }
    else {
      pred <- if (sh < cs$share) 1 else -1             # q_G* > q_B* below the critical budget (p > 1/2)
      open_n <- open_n + 1; open_ok <- open_ok + (sgn == (if (cs$flip) -pred else pred))
    }
  }
}
cat(sprintf("%d of %d cases solved (%d beyond double precision)\n", n, n + skipped, skipped))
cat(sprintf("log-odds criterion holds in %d of %d cases\n", crit_ok, n))
cat(sprintf("sign rule holds in %d of %d cases in the regions where the sign does not depend on K\n", rule_ok, rule_n))
cat(sprintf("budget-dependent region: sign agrees with the side of the critical budget in %d of %d cases\n", open_ok, open_n))
cat(sprintf("budget-dependent region: %d pairs (p, p*); critical budget K/h(p) between %.3f and %.8f; below h(p) in %d of %d pairs; above the largest sampled share %g in %d pairs\n",
            length(c_all), min(c_all), max(c_all), sum(c_all < 1), length(c_all), max(SHARES), sum(c_all > max(SHARES))))
