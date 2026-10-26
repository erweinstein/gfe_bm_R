## =============================================================================
## scripts/02_application_income_democracy.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184,
## Section 5 (income and democracy) and Supplementary Appendix S.7.
##
## Follows GFE_estimates.m (original code by Stephane Bonhomme and Elena
## Manresa) and produces:
##
##   output/table_S11_gfe_estimates.csv  Table S.XI: objective, BIC, estimates,
##                                       and three standard errors (large-T,
##                                       fixed-T, bootstrap), G = 1..15
##   output/figure1_estimates_by_G.pdf   Figure 1: estimates of the lagged
##                                       democracy and income coefficients and
##                                       of the cumulative income effect, by G,
##                                       with 95% bootstrap intervals
##   output/figure2_grouped_patterns.pdf Figure 2: G = 4 group-specific time
##                                       effects, democracy, and income by group
##   output/figure2_map.pdf              Figure 2 map (needs the authors'
##                                       world-c.dta / world-d.dta)
##   output/figureS4_profiles_G2356.pdf  Figure S.4: G = 2, 3, 5, 6
##   output/group_membership_G4.csv     group membership by country (G = 4)
##
## Model (equation (13) in the paper):
##   dem_it = theta1 dem_{i,t-1} + theta2 income_{i,t-1} + alpha_{g_i t} + v_it,
## with the cumulative ("long-run") income effect theta2 / (1 - theta1).
##
## By default the script uses the authors' optimal groupings, so it runs in
## seconds. Set `reestimate = TRUE` to recompute them from scratch with VNS
## (Algorithm 2), and `bootstrap_B > 0` to recompute bootstrap standard errors
## (the default uses the authors' bootstrap s.e., std_cluster_b in
## GFE_estimates.m, computed with 100 replications and VNS (5; 10; 5)).
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

reestimate  <- FALSE   # TRUE: recompute the groupings with VNS (10; 10; 10)
bootstrap_B <- 0       # e.g. 100 to recompute bootstrap s.e. (slow: hours for all G)
bm_dir      <- NULL      # folder with the authors' world-c.dta and world-d.dta (for the map)

p   <- load_bm_panel()
pub <- load_bm_published()
years <- p$time_labels
Gs <- 1:15

## -----------------------------------------------------------------------------
## 1. GFE estimates for G = 1..15
## -----------------------------------------------------------------------------
set.seed(1)
fits <- lapply(Gs, function(G) {
  if (G == 1) return(gfe(p, 1))
  if (reestimate) gfe(p, G, algorithm = "vns", n_starts = 10, neighmax = 10, itermax = 10)
  else gfe(p, G, groups = pub$groupings[[paste0("G", G)]])
})
bic <- gfe_bic(fits)

cum <- function(th) th[2] / (1 - th[1])
cum_grad <- function(th) c(th[2] / (1 - th[1])^2, 1 / (1 - th[1]))

tab <- do.call(rbind, lapply(Gs, function(G) {
  f <- fits[[G]]
  V <- gfe_vcov_largeT(f); se <- gfe_se_theta(f, V)
  se_c <- sqrt(gfe_delta(f, V, cum_grad(f$theta)))
  if (G >= 2) {
    Vf <- gfe_vcov_fixedT(f, gfe_bandwidth(f))
    sf <- gfe_se_theta(f, Vf); sf_c <- sqrt(gfe_delta(f, Vf, cum_grad(f$theta)))
  } else { sf <- c(NA, NA); sf_c <- NA }
  b <- pub$bootstrap_se[pub$bootstrap_se$G == G, ]
  if (bootstrap_B > 0 && G >= 2) {
    bs <- gfe_bootstrap(f, B = bootstrap_B, algorithm = "vns", fun = function(th) cum(th),
                        n_starts = 5, neighmax = 10, itermax = 5)
    b <- data.frame(se_lag_dem = bs$se[1], se_lag_income = bs$se[2], se_cumulative_income = bs$se[3])
  }
  data.frame(G = G, objective = f$objective, BIC = bic$BIC[G],
             theta1 = f$theta[1], theta2 = f$theta[2], cumulative = cum(f$theta),
             se1_largeT = se[1], se1_fixedT = sf[1], se1_boot = if (nrow(b)) b$se_lag_dem else NA,
             se2_largeT = se[2], se2_fixedT = sf[2], se2_boot = if (nrow(b)) b$se_lag_income else NA,
             sec_largeT = se_c, sec_fixedT = sf_c, sec_boot = if (nrow(b)) b$se_cumulative_income else NA)
}))
rownames(tab) <- NULL
write.csv(tab, "output/table_S11_gfe_estimates.csv", row.names = FALSE)
print(format(tab[, 1:9], digits = 3), row.names = FALSE)
cat("BIC selects G =", attr(bic, "G_BIC"), "\n")

## -----------------------------------------------------------------------------
## 2. Figure 1: estimates by G with 95% bootstrap intervals
## (G = 1: common time effects only, large-T clustered s.e., as in the paper)
## -----------------------------------------------------------------------------
se_boot <- function(col_boot, col_lt) ifelse(is.na(tab[[col_boot]]), tab[[col_lt]], tab[[col_boot]])
pdf("output/figure1_estimates_by_G.pdf", width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 1, 1), cex.lab = 1.3, cex.axis = 1.1)
plot_by_G(Gs, tab$theta1, se_boot("se1_boot", "se1_largeT"), "lagged democracy coefficient", ylim = c(0, 1))
plot_by_G(Gs, tab$theta2, se_boot("se2_boot", "se2_largeT"), "income coefficient", ylim = c(0, 0.15))
plot_by_G(Gs, tab$cumulative, se_boot("sec_boot", "sec_largeT"), "implied cumulative income effect", ylim = c(0, 0.4))
par(op); invisible(dev.off())

## -----------------------------------------------------------------------------
## 3. Figure 2: grouped patterns of democracy (G = 4)
## -----------------------------------------------------------------------------
## Relabel the published groups as in the paper: 1 = high democracy,
## 2 = low democracy, 3 = early transition, 4 = late transition. (The labels
## in GFE_estimates.m are a permutation of these; the mapping below is checked
## against the "assignment" variable of explaining_the_waves_G4.dta.)
paper_labels <- c(2, 3, 4, 1)          # published label k -> paper label
f4 <- gfe(p, 4, groups = paper_labels[pub$groupings$G4])
group_names <- c("1: high democracy", "2: low democracy", "3: early transition", "4: late transition")

dat <- attr(p, "data")
dem_by_g <- group_means(f4, p$y)
inc_by_g <- group_means(f4, p$X[, "lag_income"])

pdf("output/figure2_grouped_patterns.pdf", width = 12, height = 4.6)
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 1, 1), oma = c(2.5, 0, 0, 0), cex.lab = 1.3)
cols <- gray(c(0, 0.6, 0.4, 0.75)); ltys <- c(1, 4, 2, 1)
plot_profiles(f4$alpha, years, ylab = "time effects", col = cols, lty = ltys)
plot_profiles(dem_by_g, years, ylab = "democracy", ylim = c(0, 1), col = cols, lty = ltys)
plot_profiles(inc_by_g, years, ylab = "income", ylim = c(7, 9.5), col = cols, lty = ltys)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", legend = group_names, col = cols, lty = ltys, lwd = 3, horiz = TRUE, bty = "n", cex = 1.1)
par(op); invisible(dev.off())

membership <- data.frame(code = p$unit_labels, group = f4$groups,
                         group_name = group_names[f4$groups])
write.csv(membership, "output/group_membership_G4.csv", row.names = FALSE)
print(table(membership$group_name))

## Map of the country classification (spmap in GFE_replication_code.do),
## drawn with base graphics from the authors' world-c.dta / world-d.dta.
if (!is.null(bm_dir) && requireNamespace("haven", quietly = TRUE)) {
  wc <- as.data.frame(haven::read_dta(file.path(bm_dir, "world-c.dta")))
  wd <- as.data.frame(haven::read_dta(file.path(bm_dir, "world-d.dta")))
  wd$group <- membership$group[match(wd$GMI_CNTRY, membership$code)]
  wd$group[wd$id == 67] <- NA          # Gaza strip shares Israel's code (as in the .do file)
  fill <- c(gray(0.05), gray(0.4), gray(0.62), gray(0.82))
  pdf("output/figure2_map.pdf", width = 11, height = 6)
  par(mar = c(0, 0, 0, 0))
  plot(range(wc$`_X`, na.rm = TRUE), c(-58, max(wc$`_Y`, na.rm = TRUE)),   # no Antarctica
       type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
  for (id in unique(wc$`_ID`)) {
    pts <- wc[wc$`_ID` == id, ]
    g <- wd$group[match(id, wd$id)]
    if (is.na(g)) polygon(pts$`_X`, pts$`_Y`, col = "white", border = gray(0.6), lwd = 0.3, density = 12, angle = 45)
    else polygon(pts$`_X`, pts$`_Y`, col = fill[g], border = "white", lwd = 0.3)
  }
  legend("bottomleft", legend = c(group_names, "not in sample"), fill = c(fill, "white"),
         density = c(rep(NA, 4), 25), border = gray(0.4), bty = "n")
  invisible(dev.off())
}

## -----------------------------------------------------------------------------
## 4. Figure S.4: profiles for G = 2, 3, 5, 6
## -----------------------------------------------------------------------------
pdf("output/figureS4_profiles_G2356.pdf", width = 12, height = 14)
op <- par(mfrow = c(4, 3), mar = c(4.5, 4.5, 1, 1))
for (G in c(2, 3, 5, 6)) {
  f <- fits[[G]]
  plot_profiles(f$alpha, years, ylab = sprintf("time effects (G = %d)", G))
  plot_profiles(group_means(f, p$y), years, ylab = "democracy", ylim = c(0, 1))
  plot_profiles(group_means(f, p$X[, "lag_income"]), years, ylab = "income", ylim = c(7, 9.5))
}
par(op); invisible(dev.off())

cat("Figures and tables written to output/.\n")
