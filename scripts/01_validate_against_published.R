## =============================================================================
## scripts/01_validate_against_published.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184.
##
## Checks this R translation against the numbers published in the paper's
## Supplementary Appendix:
##
##   Part A. Given the authors' optimal groupings (hard-coded in their MATLAB
##           files), reproduce Table S.XI: objective, BIC, coefficients,
##           large-T standard errors, and fixed-T standard errors, G = 1..15;
##           the fixed-effects row; and the "True" values of Table S.III.
##   Part B. Starting from random values, check that our implementations of
##           Algorithm 1 and Algorithm 2 (VNS) find the same optima
##           (Tables S.I and S.II), and recover the authors' groupings.
##   Part C. Table S.XIV (GFE in deviations from country means) and
##           Table S.XV (first-differenced IV with estimated groups).
##
## Output: printed comparisons and output/validation_report.csv.
## Published values are reported to 3 decimals. They appear to have been
## rounded from 4-decimal MATLAB output (e.g. 0.406464 -> 0.4065 -> 0.407), so
## "match" means agreement after rounding to 3 decimals, directly or via 4.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

run_search <- TRUE       # Part B: random-start searches (a few minutes)
search_G   <- 2:5        # add 10 to reproduce Table S.II (slower in R)

p   <- load_bm_panel()
pub <- load_bm_published()$groupings
report <- list()
## Round half away from zero on the decimal representation (as when copying
## MATLAB's 4-decimal display into a table by hand).
round_hu <- function(x, d) sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d
add <- function(table, item, published, ours, digits = 3) {
  report[[length(report) + 1]] <<- data.frame(
    table = table, item = item, published = published, this_code = round(ours, 6),
    match = abs(round_hu(ours, digits) - published) < 1e-9 |
      abs(round_hu(round_hu(ours, digits + 1), digits) - published) < 1e-9)
}
cum <- function(th) th[2] / (1 - th[1])
cum_grad <- function(th) c(th[2] / (1 - th[1])^2, 1 / (1 - th[1]))

## -----------------------------------------------------------------------------
## Part A. Table S.XI given the published groupings
## -----------------------------------------------------------------------------
S11 <- data.frame(
  G   = 1:15,
  obj = c(24.301, 19.847, 16.599, 14.319, 12.593, 11.132, 10.059, 9.251, 8.426, 7.749, 7.218, 6.809, 6.391, 5.996, 5.664),
  bic = c(0.052, 0.046, 0.042, 0.039, 0.037, 0.036, 0.035, 0.035, 0.034, 0.034, 0.034, 0.034, 0.035, 0.035, 0.035),
  th1 = c(0.665, 0.601, 0.407, 0.302, 0.255, 0.465, 0.403, 0.333, 0.312, 0.277, 0.293, 0.304, 0.236, 0.237, 0.244),
  th2 = c(0.083, 0.061, 0.089, 0.082, 0.079, 0.064, 0.065, 0.070, 0.069, 0.075, 0.073, 0.074, 0.072, 0.071, 0.071),
  cum = c(0.247, 0.152, 0.151, 0.118, 0.107, 0.119, 0.108, 0.104, 0.101, 0.104, 0.104, 0.107, 0.094, 0.094, 0.094),
  lt1 = c(0.049, 0.041, 0.052, 0.054, 0.050, 0.043, 0.043, 0.044, 0.045, 0.049, 0.042, 0.044, 0.040, 0.042, 0.043),
  lt2 = c(0.014, 0.011, 0.011, 0.009, 0.010, 0.007, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.009, 0.009, 0.009),
  ltc = c(0.018, 0.021, 0.013, 0.011, 0.009, 0.011, 0.011, 0.010, 0.010, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009),
  ft1 = c(NA, 0.061, 0.083, 0.108, 0.088, 0.054, 0.074, 0.085, 0.072, 0.062, 0.062, 0.054, 0.046, 0.047, 0.046),
  ft2 = c(NA, 0.013, 0.015, 0.012, 0.012, 0.008, 0.013, 0.012, 0.010, 0.010, 0.012, 0.009, 0.010, 0.010, 0.010),
  ftc = c(NA, 0.030, 0.022, 0.021, 0.014, 0.014, 0.019, 0.014, 0.011, 0.011, 0.013, 0.010, 0.010, 0.010, 0.010)
)

fits <- lapply(1:15, function(G) gfe(p, G, groups = if (G == 1) NULL else pub[[paste0("G", G)]]))
bic <- gfe_bic(fits)
for (G in 1:15) {
  f <- fits[[G]]
  V <- gfe_vcov_largeT(f); se <- gfe_se_theta(f, V)
  add("S.XI", sprintf("G=%2d objective", G), S11$obj[G], f$objective)
  add("S.XI", sprintf("G=%2d BIC", G), S11$bic[G], bic$BIC[G])
  add("S.XI", sprintf("G=%2d theta1 (lag dem)", G), S11$th1[G], f$theta[1])
  add("S.XI", sprintf("G=%2d theta2 (income)", G), S11$th2[G], f$theta[2])
  add("S.XI", sprintf("G=%2d cumulative income", G), S11$cum[G], cum(f$theta))
  add("S.XI", sprintf("G=%2d large-T se theta1", G), S11$lt1[G], se[1])
  add("S.XI", sprintf("G=%2d large-T se theta2", G), S11$lt2[G], se[2])
  add("S.XI", sprintf("G=%2d large-T se cumulative", G), S11$ltc[G], sqrt(gfe_delta(f, V, cum_grad(f$theta))))
  if (G >= 2) {
    Vf <- gfe_vcov_fixedT(f, gfe_bandwidth(f, snap_to_grid = FALSE)); sf <- gfe_se_theta(f, Vf)
    add("S.XI", sprintf("G=%2d fixed-T se theta1", G), S11$ft1[G], sf[1])
    add("S.XI", sprintf("G=%2d fixed-T se theta2", G), S11$ft2[G], sf[2])
    add("S.XI", sprintf("G=%2d fixed-T se cumulative", G), S11$ftc[G], sqrt(gfe_delta(f, Vf, cum_grad(f$theta))))
  }
}
cat("BIC selects G =", attr(bic, "G_BIC"), "(the paper: G = 10)\n")

## Fixed-effects row of Table S.XI: two-way fixed effects = GFE with G = 1 in
## deviations from country means.
ffe <- gfe(gfe_within(p), 1)
add("S.XI", "FE objective", 17.517, ffe$objective)
add("S.XI", "FE theta1", 0.284, ffe$theta[1])
add("S.XI", "FE theta2", -0.031, ffe$theta[2])
add("S.XI", "FE cumulative income", -0.044, cum(ffe$theta))

## Table S.III, "True" columns: GFE estimates used as the Monte Carlo DGP.
for (G in c(3, 5, 10)) {
  f <- fits[[G]]
  tru <- list(`3` = c(0.407, 0.089, 0.151), `5` = c(0.255, 0.079, 0.107), `10` = c(0.277, 0.075, 0.104))[[as.character(G)]]
  add("S.III", sprintf("G=%2d true theta1", G), tru[1], f$theta[1])
  add("S.III", sprintf("G=%2d true theta2", G), tru[2], f$theta[2])
  add("S.III", sprintf("G=%2d true cumulative", G), tru[3], cum(f$theta))
}

## -----------------------------------------------------------------------------
## Part B. Do our searches find the published optima?
## -----------------------------------------------------------------------------
if (run_search) {
  set.seed(20150501)
  for (G in search_G) {
    t0 <- Sys.time()
    v <- gfe_vns(p, G, n_starts = 10, neighmax = 10, itermax = 10)
    m <- match_groups(v$groups, pub[[paste0("G", G)]], G)
    cat(sprintf("G=%2d VNS(10;10;10): objective %.6f (published grouping: %.6f), %d/10 starts reach it, %d units classified differently, %.0fs\n",
                G, v$objective, fits[[G]]$objective,
                sum(abs(v$objectives_by_start - v$objective) < 1e-8), sum(m$groups != pub[[paste0("G", G)]]),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    tab <- if (G == 10) "S.II" else if (G <= 3) "S.I" else "S.XI"
    add(tab, sprintf("G=%2d objective from VNS search", G), S11$obj[G], v$objective)
    if (G <= 3) {
      a1 <- gfe_alg1(p, G, n_starts = 1000)
      add("S.I", sprintf("G=%2d objective from Algorithm 1 (1,000 starts)", G), S11$obj[G], a1$objective)
    }
  }
}

## -----------------------------------------------------------------------------
## Part C. Extensions: Tables S.XIV and S.XV
## -----------------------------------------------------------------------------
S14 <- data.frame(G = 2:5, obj = c(12.859, 10.400, 9.221, 8.174),
                  th1 = c(0.061, -0.033, -0.072, -0.093), th2 = c(-0.038, -0.035, 0.045, -0.013),
                  se1 = c(0.049, 0.043, 0.046, 0.042), se2 = c(0.027, 0.027, 0.027, 0.026))
pw <- gfe_within(p)
for (k in seq_len(nrow(S14))) {
  G <- S14$G[k]
  f <- gfe(pw, G, groups = pub[[paste0("dev_G", G)]])
  se <- gfe_se_theta(f, gfe_vcov_largeT(f))
  add("S.XIV", sprintf("G=%d objective", G), S14$obj[k], f$objective)
  add("S.XIV", sprintf("G=%d theta1", G), S14$th1[k], f$theta[1])
  add("S.XIV", sprintf("G=%d theta2", G), S14$th2[k], f$theta[2])
  add("S.XIV", sprintf("G=%d large-T se theta1", G), S14$se1[k], se[1])
  add("S.XIV", sprintf("G=%d large-T se theta2", G), S14$se2[k], se[2])
}

S15 <- data.frame(G = 1:5, th1 = c(0.472, 0.338, 0.202, 0.065, 0.085), th2 = c(-0.075, -0.063, -0.064, -0.037, -0.078),
                  se1 = c(0.131, 0.120, 0.094, 0.089, 0.089), se2 = c(0.068, 0.062, 0.055, 0.049, 0.049))
raw <- read.csv("data/ajry_5yearpanel_extract.csv", stringsAsFactors = FALSE)
bal <- attr(p, "data")
for (k in seq_len(nrow(S15))) {
  G <- S15$G[k]
  grp <- if (G == 1) rep(1L, p$N) else pub[[paste0("dev_G", G)]]
  df <- merge(data.frame(code = bal$code, year = bal$year, group = grp[bal$unit]),
              raw[, c("code", "code_numeric", "year", "fhpolrigaug", "lrgdpch")], by = c("code", "year"))
  r <- gfe_fd_iv(df, id = "code_numeric", year = "year", y = "fhpolrigaug", x = "lrgdpch", group = "group")
  add("S.XV", sprintf("G=%d IV theta1", G), S15$th1[k], r$coef[1])
  add("S.XV", sprintf("G=%d IV theta2", G), S15$th2[k], r$coef[2])
  add("S.XV", sprintf("G=%d IV se theta1", G), S15$se1[k], r$se[1])
  add("S.XV", sprintf("G=%d IV se theta2", G), S15$se2[k], r$se[2])
}

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
rep_df <- do.call(rbind, report)
write.csv(rep_df, "output/validation_report.csv", row.names = FALSE)
cat("\nAgreement with the published values, by table:\n")
print(aggregate(match ~ table, data = rep_df, FUN = function(x) sprintf("%d of %d", sum(x), length(x))))
cat("\nItems that differ in the third decimal:\n")
print(rep_df[!rep_df$match, ], row.names = FALSE)
cat("\nNotes:",
    "\n - Fixed-T standard errors use the rule-of-thumb bandwidth without snapping it to the grid",
    "\n   linspace(0.05, 0.15, 10) used in GFE_estimates.m; the unsnapped rule reproduces Table S.XI",
    "\n   to within 0.001. The formula is sensitive to the bandwidth: moving to an adjacent grid",
    "\n   point changes several standard errors by more than that (try gfe_bandwidth(fit, TRUE)).",
    "\n - G = 1 (common time effects only): the published large-T s.e. of theta1 is 0.049; the",
    "\n   formula in GFE_estimates.m gives 0.0483. Stata's regress, vce(cluster) small-sample",
    "\n   factors would give 0.0489, so the published row may come from Stata.\n")
