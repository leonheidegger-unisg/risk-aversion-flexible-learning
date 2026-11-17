# flexible_learning.R
# Core model for the thesis "Risk Aversion and Flexible Learning: Optimal Asymmetric Signals under Limited Attention".
#
# Investor: CARA utility u(x) = -exp(-alpha x); two states G/B with prior p = P(G);
# actions S (safe, pays s) and R (risky, pays H in G and L in B), H > s > L.
# Information: binary signal z in {g, b} with accuracies q_G = P(g|G), q_B = P(b|B);
# flexible learning under the budget I(theta; z) <= K (nats).
#
# Equation numbers refer to the thesis. The solver of Appendix B.1 works in log-odds,
# uses nested bisection with 80 halvings, and evaluates every logarithm with overflow-safe formulas.
# Base R only (no packages needed).

# ----- baseline calibration (Section 3) -----
H <- 4.0; L <- 1.0; S <- 1.5
P_PRIOR <- 0.5

# ----- numerical settings (Appendix B.1) -----
Z_MAX <- 1e5          # largest log-odds used in the brackets
DELTA <- 1e-6         # mirror-point rule near the threshold
N_HALVINGS <- 80      # halvings per bisection
EPS_CERTAIN <- 1e-15  # posteriors closer than this to 0 or 1 are treated as indistinguishable from certainty

# ================================================================ overflow-safe helpers
softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))          # ln(1 + e^z)
sigmoid  <- function(z) exp(-softplus(-z))                          # 1 / (1 + e^{-z})
logit    <- function(x) log(x) - log1p(-x)                          # log-odds l(x)
hbin <- function(x) {                                               # binary entropy h(x) in nats, eq. (13)
  v <- suppressWarnings(-x * log(x) - (1 - x) * log1p(-x))
  v[x <= 0 | x >= 1] <- 0
  v
}
h_logodds <- function(z) sigmoid(z) * softplus(-z) + sigmoid(-z) * softplus(z)

# ================================================================ preferences and fixed signals (Sections 3-4)
u <- function(alpha, x) -exp(-alpha * x)

pstar <- function(alpha, H = get("H", envir = globalenv()), L = get("L", envir = globalenv()),
                  S = get("S", envir = globalenv())) {
  # action threshold (3): choose R iff the posterior is at least p*(alpha); at alpha = 0: (s-L)/(H-L)
  small <- alpha < 1e-12
  a <- ifelse(small, 1, alpha)
  ifelse(small, (S - L) / (H - L), expm1(-a * (S - L)) / expm1(-a * (H - L)))
}

mutual_info <- function(qG, qB, p = P_PRIOR) {                       # I(theta; z), eq. (14)
  pg <- p * qG + (1 - p) * (1 - qB)
  hbin(pg) - p * hbin(qG) - (1 - p) * hbin(qB)
}

posteriors <- function(qG, qB, p = P_PRIOR) {                        # P(g), pi_g, pi_b from Bayes' rule, eq. (11)
  pg <- p * qG + (1 - p) * (1 - qB)
  list(P_g = pg, pi_g = p * qG / pg, pi_b = p * (1 - qG) / (1 - pg))
}

W0 <- function(alpha, p = P_PRIOR, H = 4, L = 1, S = 1.5)            # expected utility without information, eq. (7)
  pmax(u(alpha, S), p * u(alpha, H) + (1 - p) * u(alpha, L))

W1 <- function(alpha, qG, qB, p = P_PRIOR, H = 4, L = 1, S = 1.5) {  # with the signal, eq. (8)
  uS <- u(alpha, S); uH <- u(alpha, H); uL <- u(alpha, L)
  Pgg <- p * qG; Pgb <- (1 - p) * (1 - qB); Pbg <- p * (1 - qG); Pbb <- (1 - p) * qB
  pmax(uS * (Pgg + Pgb), Pgg * uH + Pgb * uL) + pmax(uS * (Pbg + Pbb), Pbg * uH + Pbb * uL)
}

voi_risk_neutral <- function(qG, qB, p, H, L, S) {                   # limit of the monetary VOI as alpha -> 0
  Pgg <- p * qG; Pgb <- (1 - p) * (1 - qB); Pbg <- p * (1 - qG); Pbb <- (1 - p) * qB
  pmax(S * (Pgg + Pgb), Pgg * H + Pgb * L) + pmax(S * (Pbg + Pbb), Pbg * H + Pbb * L) -
    max(S, p * H + (1 - p) * L)
}

voi <- function(alpha, qG, qB, p = P_PRIOR, H = 4, L = 1, S = 1.5) { # monetary VOI = CE(W1) - CE(W0), eq. (9)
  small <- alpha < 1e-12
  a <- ifelse(small, 1, alpha)
  cara <- (-log(-W1(a, qG, qB, p, H, L, S)) + log(-W0(a, p, H, L, S))) / a
  ifelse(small, voi_risk_neutral(qG, qB, p, H, L, S), cara)
}

alpha_at <- function(target, H = 4, L = 1, S = 1.5) {               # alpha with p*(alpha) = target
  if (target <= pstar(0, H, L, S)) return(0)
  lo <- 1e-9; hi <- 200
  for (i in 1:200) {
    mid <- 0.5 * (lo + hi)
    if (pstar(mid, H, L, S) < target) lo <- mid else hi <- mid
  }
  0.5 * (lo + hi)
}

alpha_indifference <- function(p = P_PRIOR, H = 4, L = 1, S = 1.5) { # alpha_p with p*(alpha_p) = p (Section 3.5)
  if (p <= (S - L) / (H - L)) stop("no alpha_p > 0: the prior does not exceed the risk-neutral threshold")
  alpha_at(p, H, L, S)
}

# ================================================================ flexible learning: solver of Appendix B.1
F_eq <- function(z, ps) ps * z - softplus(z)   # f(sigma(z)) = p* z - ln(1 + e^z); condition (24) reads F(z_b) = F(z_g)

bisect <- function(fn, lo, hi, n = N_HALVINGS) {
  flo <- fn(lo); fhi <- fn(hi)
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)
  if ((flo > 0) == (fhi > 0)) stop("root not bracketed: parameters too extreme for double precision")
  for (i in seq_len(n)) {
    mid <- 0.5 * (lo + hi); fm <- fn(mid)
    if ((fm > 0) == (flo > 0)) { lo <- mid; flo <- fm } else hi <- mid
  }
  0.5 * (lo + hi)
}

partner <- function(z, ps) {
  # inner step: the log-odds on the other side of l(p*) with the same value of F (equidistance (24)).
  # Same bisection as bisect(), written out as a loop for speed.
  zs <- logit(ps)
  if (abs(z - zs) < DELTA) return(2 * zs - z)   # ill-conditioned near the threshold: mirror point, error O(DELTA^2)
  Fz <- function(y) ps * y - (max(y, 0) + log1p(exp(-abs(y))))
  target <- Fz(z)
  if (z < zs) { lo <- zs; hi <- Z_MAX } else { lo <- -Z_MAX; hi <- zs }   # F is strictly decreasing above l(p*)
  flo <- Fz(lo) - target; fhi <- Fz(hi) - target
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)
  if ((flo > 0) == (fhi > 0)) stop("root not bracketed: parameters too extreme for double precision")
  for (i in seq_len(N_HALVINGS)) {
    mid <- 0.5 * (lo + hi); fm <- ps * mid - (max(mid, 0) + log1p(exp(-abs(mid)))) - target
    if ((fm > 0) == (flo > 0)) { lo <- mid; flo <- fm } else hi <- mid
  }
  0.5 * (lo + hi)
}

info_pair <- function(p, zb, zg) {             # information used, eq. (14), and P(g) for the pair (zb, zg)
  pib <- sigmoid(zb); pig <- sigmoid(zg)
  if (pig - pib <= 0) return(c(I = 0, P_g = 1)) # degenerate pair at the uninformative end of the bracket
  Pg <- (p - pib) / (pig - pib)
  c(I = hbin(p) - Pg * h_logodds(zg) - (1 - Pg) * h_logodds(zb), P_g = Pg)
}

optimal_signal <- function(p, ps, K) {
  # Solve system (27): F(z_b) = F(z_g) and I = K, by nested bisection in log-odds (Appendix B.1).
  # Valid for any prior p and threshold ps in (0, 1) and 0 < K < h(p).
  if (!(K > 0 && K < hbin(p))) stop("need 0 < K < h(p)")
  zs <- logit(ps); m <- logit(p)
  z_up <- if (p <= ps) m else partner(m, ps)                  # upper end of the z_b bracket: uninformative point
  F_top <- F_eq(0.99 * Z_MAX, ps)                              # lower end: most extreme representable z_b
  z_lo <- if (F_eq(-0.99 * Z_MAX, ps) >= F_top) -0.99 * Z_MAX else
    bisect(function(z) F_eq(z, ps) - F_top, -Z_MAX, zs)
  zb <- bisect(function(z) info_pair(p, z, partner(z, ps))[["I"]] - K, z_lo, z_up)
  zg <- partner(zb, ps)
  pib <- sigmoid(zb); pig <- sigmoid(zg)
  if (min(sigmoid(-abs(zb)), sigmoid(-abs(zg))) < EPS_CERTAIN)
    stop("an optimal posterior lies within 1e-15 of certainty: beyond double precision")
  Pg <- info_pair(p, zb, zg)[["P_g"]]
  qG <- Pg * pig / p
  qB <- 1 - Pg * sigmoid(-zg) / (1 - p)
  I_channel <- mutual_info(qG, qB, p)                          # independent recomputation from the accuracies
  list(pi_b = pib, pi_g = pig, z_b = zb, z_g = zg, P_g = Pg, qG = qG, qB = qB,
       lam = 1 / (zg - zb),                                    # shadow value (25)
       I = I_channel, eq_residual = abs(F_eq(zb, ps) - F_eq(zg, ps)),
       info_rel_dev = abs(I_channel - K) / K)
}

J_star <- function(p, ps, K) {                 # maximal normalised gain J*(K), eq. (21)
  s <- optimal_signal(p, ps, K)
  s$P_g * (s$pi_g - ps) - max(p - ps, 0)
}

voi_flexible <- function(alpha, p, K, H = 4, L = 1, S = 1.5) {
  sol <- optimal_signal(p, pstar(alpha, H, L, S), K)
  list(voi = voi(alpha, sol$qG, sol$qB, p, H, L, S), sol = sol)
}

# ================================================================ brute-force check (independent of the first-order conditions)
grid_optimum <- function(alpha, p, K, n = 801, H = 4, L = 1, S = 1.5) {
  g <- seq(0, 1, length.out = n)
  QG <- matrix(g, n, n); QB <- matrix(g, n, n, byrow = TRUE)
  feas <- (mutual_info(QG, QB, p) <= K) & (QG + QB >= 1)
  val <- ifelse(feas, W1(alpha, QG, QB, p, H, L, S), -Inf)
  i <- which.max(val)
  c(qG = QG[i], qB = QB[i], W = val[i])
}

K0 <- log(2) - hbin(0.75)                       # baseline budget (16)

# shared plotting helpers -------------------------------------------------------------------
ALARM <- "#c8501d"; CLEAR <- "#2a8a4a"; FLEX <- "#1b6ca8"
save_fig <- function(name, width, height, draw) {
  dir.create("figures", showWarnings = FALSE)
  pdf(file.path("figures", paste0(name, ".pdf")), width = width, height = height, family = "serif"); draw(); dev.off()
  png(file.path("figures", paste0(name, ".png")), width = width, height = height, units = "in", res = 150, family = "serif")
  draw(); invisible(dev.off())
}
