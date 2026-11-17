# Risk Aversion and Flexible Learning — R code

R code for the bachelor thesis *Risk Aversion and Flexible Learning: Optimal Asymmetric Signals under Limited Attention*
(Leonard Heidegger, University of St. Gallen, 2026). The scripts reproduce every figure and every number reported in the thesis.

## Requirements

R (tested with R 4.3.3). Only base R is used, so no packages need to be installed.

## Files

| File | Content | Thesis |
|---|---|---|
| `flexible_learning.R` | model, value of information and the solver of system (27); used by all other scripts | Sections 3–6, Appendix B.1 |
| `fig_benchmark.R` | Figure 1: threshold and value of information of the fixed symmetric signal | Section 3.5 |
| `fig_asymmetric.R` | Figure 2: symmetric, alarm and all-clear signals | Section 4.3 |
| `fig_flexible.R` | Figure 3: optimal signal under the budget K₀ | Section 5.5 |
| `fig_flex_vs_fixed.R` | Figure 4: optimal flexible signal versus the fixed symmetric signal | Section 6.2 |
| `fig_robustness.R` | Figure 5: value of the optimal signal for three budgets | Section 6.3 |
| `reproduce_numbers.R` | prints the numerical values quoted in Sections 3–6 and the precision limits of Appendix B.1 | Sections 3–6, Appendix B.1 |
| `validate_flexible.R` | prints the solver checks: residuals, brute-force grid comparison, shadow value, location of the peak, rise and fall of the value around α_p, degree of asymmetry across budgets | Appendix B.1–B.2 |
| `validate_direction.R` | prints the checks on the direction of the asymmetry: log-odds criterion, sign rules, critical budgets | Appendix B.2 |
| `figures/` | Figures 1–5 as PDF and PNG files, written by the figure scripts | Figures 1–5 |

## Usage

Run each script from this folder, for example in a terminal:

```
Rscript reproduce_numbers.R
Rscript fig_robustness.R
```

or in RStudio after setting the working directory to this folder (`setwd(...)`) and calling `source("reproduce_numbers.R")`.
The figure scripts write PDF and PNG files to the folder `figures`. The figure scripts and `reproduce_numbers.R` take a few
seconds each, `validate_direction.R` under half a minute and `validate_flexible.R` a few minutes.

## Model and calibration

CARA utility u(x) = −e^(−αx), baseline calibration (H, L, s, p) = (4, 1, 1.5, 0.5) and baseline information budget
K₀ = ln 2 − h(0.75) ≈ 0.131 nats. All logarithms are natural logarithms, so information is measured in nats.

## Numerical method

Proposition 2 of the thesis reduces the investor's problem to two equations in the two posterior beliefs: the
Kullback–Leibler equidistance condition (24) and the binding information budget. Both are solved by nested bisection
in log-odds, with 80 halvings per bisection, and all logarithms are evaluated with overflow-safe formulas. Within 10⁻⁶
of the threshold, the inner step uses the mirror point in log-odds. When an optimal posterior lies within 10⁻¹⁵ of
certainty, the solver stops with an error instead of returning an inaccurate value; in the baseline calibration this
happens for risk aversion above about 7.3. Details are given in Appendix B of the thesis.
