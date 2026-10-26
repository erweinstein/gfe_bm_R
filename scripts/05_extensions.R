## =============================================================================
## scripts/05_extensions.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184,
## Supplementary Appendix S.7.3 (robustness of the grouped patterns).
##
## Follows (original code by Stephane Bonhomme and Elena Manresa):
##   GFE_estimates.m (continent groups)       Figure S.5
##   GFE_deviations_estimates.m               Table S.XIV, Figure S.6 (top)
##   GFE_replication_code_IV.do               Table S.XV
##   GFE_twolayers_computation.m / _estimates.m   Figure S.6 (bottom)
##   Heterogeneous_coefficientsA.m / B.m      Table S.XVI, Figure S.7
##
## Each exercise asks the same question in a different way: is the division of
## countries into "stable", "early transition", and "late transition" paths a
## robust feature of the data, or an artifact of the baseline specification?
##
## By default the script evaluates the published groupings (and, for the
## heterogeneous-coefficient models, the groupings found by this R code; see
## below), so it runs in seconds. The flags turn on the searches.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

search_deviations <- FALSE   # TRUE: re-estimate the deviation groupings with VNS
search_twolayer   <- FALSE   # TRUE: search the (5, 2, 2) two-layer grouping (slow)
search_hetero     <- FALSE   # TRUE: search the heterogeneous-coefficient groupings (slow)
hetero_starts     <- 20000   # the paper used 1,000,000 random starts

p     <- load_bm_panel()
pub   <- load_bm_published()$groupings
years <- p$time_labels
N <- p$N; T <- p$T
cum <- function(th) th[2] / (1 - th[1])
cum_grad <- function(th) c(th[2] / (1 - th[1])^2, 1 / (1 - th[1]))

## Baseline G = 4 grouping with the paper's labels (see script 02):
## 1 = high democracy, 2 = low democracy, 3 = early transition, 4 = late transition.
base4 <- c(2, 3, 4, 1)[pub$G4]
base_names <- c("high democracy", "low democracy", "early transition", "late transition")

## Panels of three plots (time effects, democracy, income), as in Figure S.4.
three_panels <- function(alpha, groups, file, legend_text, col, lty, alpha_lab = "time effects",
                         alpha_ylim = NULL) {
  pdf(file, width = 12, height = 4.6)
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 1, 1), oma = c(2.5, 0, 0, 0), cex.lab = 1.3)
  G <- ncol(alpha)
  f <- list(panel = p, groups = groups, G = G)
  plot_profiles(alpha, years, ylab = alpha_lab, col = col, lty = lty,
                ylim = if (is.null(alpha_ylim)) range(alpha) else alpha_ylim)
  plot_profiles(group_means(f, p$y), years, ylab = "democracy", ylim = c(0, 1), col = col, lty = lty)
  plot_profiles(group_means(f, p$X[, "lag_income"]), years, ylab = "income", ylim = c(7, 10), col = col, lty = lty)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("bottom", legend = legend_text, col = col, lty = lty, lwd = 3, horiz = TRUE, bty = "n",
         cex = if (G > 5) 0.8 else 1.1)
  par(op); invisible(dev.off())
}

## =============================================================================
## 1. Continent-specific time effects (Figure S.5)
## =============================================================================
## Groups fixed a priori to continents, instead of estimated: a natural
## benchmark for the GFE classification. (Continents are coded 1 Africa,
## 2 Asia incl. Australia and New Zealand, 3 Europe, 4 South America,
## 5 North America incl. Mexico, as in GFE_estimates.m.)
continents <- c("Africa", "Asia", "Europe", "South America", "North America")
fc <- gfe(p, 5, groups = pub$continent5)
cat("Continent-specific time effects: theta =", round(fc$theta, 3),
    "; objective =", round(fc$objective, 3), "(GFE with G = 5:", round(gfe(p, 5, groups = pub$G5)$objective, 3), ")\n")
three_panels(fc$alpha, fc$groups, "output/figureS5_continents.pdf", continents,
             col = gray(c(0, 0.3, 0.5, 0.65, 0.8)), lty = c(1, 2, 4, 1, 2))

## =============================================================================
## 2. Country fixed effects plus grouped patterns (Table S.XIV, Figure S.6 top)
## =============================================================================
## Model (5): y_it = x_it' theta + alpha_{g_i t} + eta_i + v_it.
## GFE applied to deviations from country means (S.21) removes eta_i.
## G = 1 is two-way fixed effects.
pw <- gfe_within(p)
set.seed(5)
dev_fits <- lapply(1:5, function(G) {
  if (G == 1) return(gfe(pw, 1))
  if (search_deviations) gfe(pw, G, algorithm = "vns", n_starts = 10, neighmax = 10, itermax = 10)
  else gfe(pw, G, groups = pub[[paste0("dev_G", G)]])
})
tabS14 <- do.call(rbind, lapply(dev_fits, function(f) {
  V <- gfe_vcov_largeT(f); se <- gfe_se_theta(f, V)
  data.frame(G = f$G, objective = f$objective, theta1 = f$theta[1], se1 = se[1],
             theta2 = f$theta[2], se2 = se[2], cumulative = cum(f$theta),
             se_cum = sqrt(gfe_delta(f, V, cum_grad(f$theta))))
}))
cat("\nTable S.XIV: GFE in deviations from country means\n")
print(format(tabS14, digits = 3), row.names = FALSE)
write.csv(tabS14, "output/tableS14_deviations.csv", row.names = FALSE)
## Rows G = 2..5 match Table S.XIV. For G = 1 (two-way fixed effects) the
## published standard errors (0.058, 0.049, 0.069) are larger than the formula of
## GFE_deviations_estimates.m gives here (0.053, 0.045, 0.064), plausibly
## because that row comes from Stata with country dummies counted in the
## degrees of freedom.

## The paper: with G = 3, "the match with the classification without country
## fixed-effects and G = 4 is perfect for 80 out of the 90 countries", the
## stable group mostly comprising high- and low-democracy countries.
dev3 <- dev_fits[[3]]
xt <- table(deviations_G3 = dev3$groups, baseline_G4 = base_names[base4])
print(xt)
stable <- which.max(xt[, "high democracy"] + xt[, "low democracy"])
early  <- which.max(xt[, "early transition"])
late   <- which.max(xt[, "late transition"])
cat("Countries classified alike:",
    xt[stable, "high democracy"] + xt[stable, "low democracy"] + xt[early, "early transition"] +
      xt[late, "late transition"], "of", N, "(paper: 80)\n")
dev_names <- character(3); dev_names[c(stable, early, late)] <- c("stable", "early transition", "late transition")
three_panels(dev3$alpha, dev3$groups, "output/figureS6_top_deviations_G3.pdf", dev_names,
             col = gray(c(0, 0.4, 0.7))[order(c(stable, early, late))], lty = 1,
             alpha_lab = "time effects (deviations)")

## =============================================================================
## 3. First-differenced IV with estimated groups (Table S.XV)
## =============================================================================
## Delta y_it = theta1 Delta y_{i,t-1} + theta2 Delta x_{i,t-1} + mu_{g_i t} + Delta v_it,
## instrumenting Delta y_{i,t-1} by y_{i,t-2}; groups from GFE in deviations.
## IV corrects the small-T (Nickell) bias of the within estimator of theta1.
raw <- read.csv("data/ajry_5yearpanel_extract.csv", stringsAsFactors = FALSE)
bal <- attr(p, "data")
tabS15 <- do.call(rbind, lapply(1:5, function(G) {
  grp <- dev_fits[[G]]$groups
  df <- merge(data.frame(code = bal$code, year = bal$year, group = grp[bal$unit]),
              raw[, c("code", "code_numeric", "year", "fhpolrigaug", "lrgdpch")], by = c("code", "year"))
  r <- gfe_fd_iv(df, id = "code_numeric", year = "year", y = "fhpolrigaug", x = "lrgdpch", group = "group")
  data.frame(G = G, theta1 = r$coef[1], se1 = r$se[1], theta2 = r$coef[2], se2 = r$se[2],
             cumulative = cum(r$coef), n_obs = r$n)
}))
rownames(tabS15) <- NULL
cat("\nTable S.XV: first-differenced IV with estimated groups\n")
print(format(tabS15, digits = 3), row.names = FALSE)
cat("(The published standard errors of the cumulative effect are not reproduced here: the\n",
    " code that computed them is not part of the replication files.)\n")
write.csv(tabS15, "output/tableS15_iv.csv", row.names = FALSE)

## =============================================================================
## 4. Two layers of groups (Figure S.6, bottom)
## =============================================================================
## Model (S.40): three "layer-1" groups with their own time paths, split into
## (5, 2, 2) subgroups that differ from the layer-1 path by a constant:
##   alpha_{(g,h) t} = a_{g t} + b_{g h}.
H <- c(5, 2, 2)
tl <- if (search_twolayer) gfe_twolayer(p, H, n_starts = 1000) else gfe_twolayer(p, H, groups = pub$twolayer_522)
cat("\nTwo-layer model (5, 2, 2): objective", round(tl$objective, 4), "; theta =", round(tl$theta, 3),
    "; cumulative income effect =", round(cum(tl$theta), 3), "\n")

## Check: the constrained least squares of GFE_twolayers_computation.m equals
## the closed form of GFE_twolayers_estimates.m (two-way demeaning within each
## layer-1 group: y - ybar_{g t} - ybar_{(g,h)} + ybar_g).
cell_mean <- function(v, key) ave(v, key)
g1 <- tl$layer1[p$unit]; sg <- tl$groups[p$unit]
dmean <- function(v) v - cell_mean(v, paste(g1, p$time)) - cell_mean(v, sg) + cell_mean(v, g1)
yd <- dmean(p$y); Xd <- apply(p$X, 2, dmean)
theta_closed <- drop(solve(crossprod(Xd), crossprod(Xd, yd)))
stopifnot(max(abs(theta_closed - tl$theta)) < 1e-8)

## The paper: the two-layer and deviation classifications agree "except in
## five cases".
cat("Two-layer (layer 1) vs deviations (G = 3):\n"); print(table(two_layer = tl$layer1, deviations = dev3$groups))
tl_names <- paste0("(", tl$layer1[match(seq_len(sum(H)), tl$groups)], ",",
                   tl$layer2[match(seq_len(sum(H)), tl$groups)], ")")
## Line styles of GFE_twolayers_estimates.m
three_panels(tl$c, tl$groups, "output/figureS6_bottom_twolayer.pdf", tl_names,
             col = gray(c(0, 0.35, 0.75, 0.65, 0.2, 0.35, 0.35, 0.65, 0.65)),
             lty = c(1, 1, 1, 1, 1, 4, 2, 4, 2), alpha_ylim = c(-0.6, 0.5))

## =============================================================================
## 5. Heterogeneous coefficients (Table S.XVI, Figure S.7)
## =============================================================================
## Model (7): y_it = x_it' theta_{g_i} + alpha_{g_i t} + v_it, G = 4.
##   Version B: income coefficient group-specific, lagged democracy common.
##   Version A: both coefficients group-specific.
## Computation: a counterpart of Algorithm 1 from random initial groupings
## (Heterogeneous_coefficientsA/B.m). The paper used 1,000,000 starts. With
## far fewer starts, a final single-unit local search (polish_top) helps;
## see R/gfe_extensions.R. The authors did not publish these groupings; the
## groupings in data/rcode_hetero_groupings.csv were found by this R code
## (search_hetero = TRUE) and reproduce Table S.XVI as noted below.
found <- read.csv("data/rcode_hetero_groupings.csv", comment.char = "#", stringsAsFactors = FALSE)
het_spec <- list(B = 2, A = 1:2)
het <- lapply(names(het_spec), function(v) {
  if (search_hetero) {
    set.seed(100)
    gfe_hetero(p, 4, hetero = het_spec[[v]], n_starts = hetero_starts, polish_top = 30)
  } else {
    gfe_hetero(p, 4, hetero = het_spec[[v]], start_groups = list(found[[paste0("hetero_", v)]]))
  }
})
names(het) <- names(het_spec)

## Relabel the groups to match the baseline G = 4 classification, so that
## Group 1 is "high democracy", etc., as in Table S.XVI.
tabS16 <- do.call(rbind, lapply(names(het), function(v) {
  h <- het[[v]]
  m <- match_groups(h$groups, base4, 4)
  ord <- order(m$map)                       # estimated label for paper group 1..4
  if (v == "B") {
    th1 <- rep(h$theta_common, 4); s1 <- rep(h$se_theta_common, 4)
    th2 <- h$beta[ord, 1]; s2 <- h$se_beta[ord, 1]
  } else {
    th1 <- h$beta[ord, 1]; s1 <- h$se_beta[ord, 1]
    th2 <- h$beta[ord, 2]; s2 <- h$se_beta[ord, 2]
  }
  data.frame(version = v, group = 1:4, group_name = base_names, objective = h$objective,
             theta1 = th1, se1 = s1, theta2 = th2, se2 = s2, cumulative = th2 / (1 - th1),
             n = tabulate(m$groups, 4), same_as_baseline = tapply(m$groups == base4, m$groups, sum))
}))
rownames(tabS16) <- NULL
cat("\nTable S.XVI: heterogeneous coefficients (G = 4)\n")
print(format(tabS16, digits = 3), row.names = FALSE)
cat("Published, version B: theta1 = 0.288 (0.054); theta2 = 0.103 (0.019), 0.047 (0.014),",
    "0.087 (0.018), 0.082 (0.013).\n")
cat("Published, version A: theta1 = 0.644, 0.319, 0.016, 0.248; theta2 = 0.070, 0.041, 0.122, 0.090.\n")
write.csv(tabS16, "output/tableS16_heterogeneous.csv", row.names = FALSE)

## Figure S.7: democracy and income by group in both versions.
pdf("output/figureS7_heterogeneous.pdf", width = 8, height = 8)
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))
cols <- gray(c(0, 0.6, 0.4, 0.75)); ltys <- c(1, 4, 2, 1)
for (v in c("B", "A")) {
  m <- match_groups(het[[v]]$groups, base4, 4)
  f <- list(panel = p, groups = m$groups, G = 4)
  plot_profiles(group_means(f, p$y), years, ylab = "democracy", ylim = c(0, 1), col = cols, lty = ltys)
  title(if (v == "B") "heterogeneous income coefficient" else "heterogeneous coefficients", cex.main = 1)
  plot_profiles(group_means(f, p$X[, "lag_income"]), years, ylab = "income", ylim = c(7, 10), col = cols, lty = ltys)
}
par(op); invisible(dev.off())

cat("\nFigures and tables written to output/.\n")
