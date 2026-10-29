## =============================================================================
## scripts/04_montecarlo.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184,
## Supplementary Appendix S.6 (simulation exercises).
##
## Follows the authors' Monte Carlo code (original code by Stephane Bonhomme
## and Elena Manresa):
##   MonteCarlo/DGP_final.m                    data generating process
##   MonteCarlo/Executable Files/MC1.exe       GFE estimation in each sample
##   MonteCarlo/Executable Files/MC_bootstrap.exe  bootstrap in each sample
##   MonteCarlo/Processing_bootstrap_G10_dgp1.m    bootstrap coverage
##   MonteCarlo/BIC_G3_MC.m                    BIC choice of G (Table S.VIII)
##
## Design (Section S.6.1). "We design [the exercise] to mimic the cross-country
## data set that we use in the empirical application": N = 90, T = 7, and
##
##   y_it = theta1 y_{i,t-1} + theta2 x_it + alpha_{g_i t} + v_it ,
##
## with x_it = lagged log income, taken from the data, and (theta, alpha, g_i)
## fixed at their GFE estimates on the data with G0 groups. The outcome is
## generated recursively from each country's observed 1970 lagged democracy.
## Errors:
##   dgp = "iid"   : i.i.d. normal with variance Q / (NT - G T - N - K)
##                   (sigma2_naive in DGP_final.m), Tables S.III-S.IV
##   dgp = "wild"  : unit-level wild bootstrap of the GFE residuals,
##                   v*_it = eta_i * v_hat_it, eta_i ~ N(0, 1), as coded in
##                   DGP_final.m (DGP = 2), Tables S.VI-S.VII
##   dgp = "resample" : unit-specific residual sequences drawn with replacement,
##                   as DESCRIBED in the text of Section S.6.2 ("Nonnormal
##                   design"). The code and the text differ; we offer both.
##
## For each simulated sample we compute (Tables S.III and S.IV):
##   - the GFE estimates (VNS with (5; 10; 5), as in the paper),
##   - the misclassification frequency, after relabeling the estimated groups
##     to best match the true ones. The paper takes the minimum over all G!
##     permutations (over 500,000 random permutations when G = 10); the
##     Hungarian algorithm (match_groups) finds the exact minimum directly,
##   - large-T clustered and fixed-T (Pollard) standard errors and the
##     coverage of the 95% normal intervals,
##   - optionally (bootstrap_B > 0) bootstrap standard errors and the coverage
##     of 95% percentile intervals (Processing_bootstrap_G10_dgp1.m).
##
## Cost. The paper uses 1,000 replications, and 100 bootstrap replications per
## sample with Algorithm 1 and 1,000 starts. The defaults below are much
## smaller so that the script runs in minutes; the results then carry
## simulation noise of the order 1/sqrt(n_rep). Set `cores` to use several
## processors. Each replication sets its own seed (seed + r), so the results
## do not depend on the number of cores.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

G0          <- 3          # true number of groups: 3, 5, or 10 in the paper
dgp         <- "iid"      # "iid", "wild", or "resample" (see above)
n_rep       <- 50         # 1,000 in the paper
vns_tuning  <- list(n_starts = 5, neighmax = 10, itermax = 5)   # (5; 10; 5)
bootstrap_B <- 0          # 100 in the paper (with Algorithm 1, 1,000 starts)
boot_tuning <- list(algorithm = "alg1", n_starts = 100)          # paper: 1,000 starts
run_bic     <- FALSE      # Table S.VIII: BIC choice of G across simulated samples
bic_G       <- if (G0 == 10) 3:15 else 1:7   # G values compared (as in BIC_G3_MC.m, BIC_G10_MC.m)
bic_Gmax    <- 15         # G used for sigma2_hat in the BIC penalty (footnote 12)
## Where the G = 15 grouping for sigma2_hat comes from:
##   "own"   each sample's own G = 15 fit, as (S20) and footnote 12 describe (and as
##           BIC_G10_MC.m does);
##   "g0_10" as BIC_G3_MC.m does for G0 = 3: it loads assignmentG010G15.txt, the
##           G = 15 groupings estimated on the G0 = 10 samples, and applies them to
##           the G0 = 3 samples. They do not fit the G0 = 3 data, so sigma2_hat is
##           about twice as large, and the penalty with it. This reproduces the
##           published Table S.VIII (G = 3 in 98%); "own" chooses G = 5-6.
bic_sigma2  <- "own"
cores       <- 1
seed        <- 20150501

p   <- load_bm_panel()
pub <- load_bm_published()
N <- p$N; T <- p$T; K <- p$K

## -----------------------------------------------------------------------------
## 1. The DGP: parameters fixed at their GFE values on the data
## -----------------------------------------------------------------------------
truth <- gfe(p, G0, groups = pub$groupings[[paste0("G", G0)]])
cum   <- function(th) th[2] / (1 - th[1])
cum_grad <- function(th) c(th[2] / (1 - th[1])^2, 1 / (1 - th[1]))
theta0 <- c(truth$theta, cumulative = cum(truth$theta))
sigma2_naive <- truth$objective / (N * T - G0 * T - N - K)
truth10 <- NULL; sigma2_10 <- NULL
if (run_bic && bic_sigma2 == "g0_10") {
  truth10   <- gfe(p, 10, groups = pub$groupings$G10)
  sigma2_10 <- truth10$objective / (N * T - 10 * T - N - K)
}
cat(sprintf("True values (G0 = %d): theta1 = %.4f, theta2 = %.4f, cumulative = %.4f; sigma = %.4f\n",
            G0, theta0[1], theta0[2], theta0[3], sqrt(sigma2_naive)))

#' Simulate one data set from the DGP of DGP_final.m.
#' Returns a balanced gfe_panel with the same layout as the data.
simulate_panel <- function(truth, dgp = c("iid", "wild", "resample"), sigma2_naive = get("sigma2_naive", envir = globalenv())) {
  dgp <- match.arg(dgp)
  p <- truth$panel; N <- p$N; T <- p$T
  ehat <- matrix(truth$residuals, T, N)                  # T x N GFE residuals
  err <- switch(dgp,
    ## V / std(V) rescales the draws to have sample variance exactly 1
    iid = { V <- rnorm(N * T); matrix(sqrt(sigma2_naive) * V / sd(V), T, N) },
    wild = sweep(ehat, 2, rnorm(N), "*"),
    resample = ehat[, sample.int(N, N, replace = TRUE), drop = FALSE])
  lag_inc <- matrix(p$X[, "lag_income"], T, N)
  delta   <- truth$alpha[, truth$groups, drop = FALSE]   # T x N: alpha_{g_i t}
  dm <- matrix(0, T + 1, N)
  dm[1, ] <- matrix(p$X[, "lag_dem"], T, N)[1, ]          # observed 1970 lagged democracy
  for (t in 1:T)
    dm[t + 1, ] <- truth$theta[1] * dm[t, ] + truth$theta[2] * lag_inc[t, ] + delta[t, ] + err[t, ]
  gfe_panel(as.vector(dm[-1, ]),
            cbind(lag_dem = as.vector(dm[-(T + 1), ]), lag_income = as.vector(lag_inc)),
            N = N, T = T, unit_labels = p$unit_labels, time_labels = p$time_labels)
}

## -----------------------------------------------------------------------------
## 2. One replication
## -----------------------------------------------------------------------------
one_replication <- function(r) {
  set.seed(seed + r)
  ps  <- simulate_panel(truth, dgp)
  fit <- do.call(gfe, c(list(panel = ps, G = G0, algorithm = "vns"), vns_tuning))
  est <- c(fit$theta, cum(fit$theta))

  ## Misclassification, after the best relabeling of the estimated groups.
  mg <- match_groups(fit$groups, truth$groups, G0)
  alpha_rel <- matrix(NA, T, G0); alpha_rel[, mg$map] <- fit$alpha

  ## Standard errors: large-T clustered (S.3) and fixed-T (Pollard).
  se_of <- function(V) c(gfe_se_theta(fit, V), sqrt(gfe_delta(fit, V, cum_grad(fit$theta))))
  se_lt <- se_of(gfe_vcov_largeT(fit))
  se_ft <- se_of(gfe_vcov_fixedT(fit, gfe_bandwidth(fit)))
  cover <- function(se) as.numeric(abs(est - theta0) <= qnorm(0.975) * se)

  out <- list(est = est, misclass = mg$misclassification, alpha = alpha_rel,
              se_largeT = se_lt, se_fixedT = se_ft,
              cov_largeT = cover(se_lt), cov_fixedT = cover(se_ft))

  if (bootstrap_B > 0) {
    bs <- do.call(gfe_bootstrap, c(list(fit = fit, B = bootstrap_B, fun = cum), boot_tuning))
    q <- bs$percentile
    out$se_boot  <- bs$se
    out$cov_boot <- as.numeric(q[1, ] < theta0 & theta0 < q[2, ])
  }

  if (run_bic) {
    Q <- sapply(bic_G, function(G) if (G == 1) gfe(ps, 1)$objective else
      do.call(gfe, c(list(panel = ps, G = G, algorithm = "vns"), vns_tuning))$objective)
    if (bic_sigma2 == "g0_10") {
      ## Emulates BIC_G3_MC.m: a G = 15 grouping estimated on an independent G0 = 10
      ## sample, applied to this G0 = 3 sample.
      p10 <- simulate_panel(truth10, dgp, sigma2_10)
      g15 <- do.call(gfe, c(list(panel = p10, G = bic_Gmax, algorithm = "vns"), vns_tuning))$groups
      Qmax <- gfe(ps, bic_Gmax, groups = g15)$objective
    } else {
      Qmax <- do.call(gfe, c(list(panel = ps, G = bic_Gmax, algorithm = "vns"), vns_tuning))$objective
    }
    out$bic_Q <- Q
    out$bic_sigma2 <- Qmax / (N * T - bic_Gmax * T - N - K)
  }
  out
}

## -----------------------------------------------------------------------------
## 3. Run the replications (optionally in parallel; works on Windows)
## -----------------------------------------------------------------------------
t0 <- Sys.time()
if (cores > 1) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterCall(cl, function(files) { for (f in files) source(f); NULL }, .gfe_source_files())
  parallel::clusterExport(cl, c("truth", "theta0", "sigma2_naive", "simulate_panel", "cum", "cum_grad",
                                "G0", "dgp", "vns_tuning", "bootstrap_B", "boot_tuning", "run_bic",
                                "bic_G", "bic_Gmax", "bic_sigma2", "truth10", "sigma2_10", "seed", "N", "T", "K"))
  sims <- parallel::parLapply(cl, seq_len(n_rep), one_replication)
  parallel::stopCluster(cl)
} else {
  sims <- lapply(seq_len(n_rep), function(r) {
    if (r %% 10 == 0) message("replication ", r, " of ", n_rep)
    one_replication(r)
  })
}
cat("Elapsed:", format(Sys.time() - t0, digits = 3), "\n")

## -----------------------------------------------------------------------------
## 4. Tables S.III / S.VI (bias) and S.IV / S.VII (inference)
## -----------------------------------------------------------------------------
stack <- function(name) do.call(rbind, lapply(sims, `[[`, name))
pars  <- c("theta1", "theta2", "cumulative")
est   <- stack("est"); colnames(est) <- pars

bias_tab <- data.frame(parameter = pars, true = theta0, gfe_mean = colMeans(est),
                       gfe_median = apply(est, 2, median), mc_sd = apply(est, 2, sd))
cat(sprintf("\nBias of the GFE estimator (G0 = %d, dgp = %s, %d replications)\n", G0, dgp, n_rep))
print(format(bias_tab, digits = 3), row.names = FALSE)
cat(sprintf("Average misclassification frequency: %.2f%%\n", 100 * mean(sapply(sims, `[[`, "misclass"))))

inf_tab <- data.frame(parameter = pars,
  se_largeT_median = apply(stack("se_largeT"), 2, median),
  se_fixedT_median = apply(stack("se_fixedT"), 2, median),
  se_boot_median   = if (bootstrap_B > 0) apply(stack("se_boot"), 2, median) else NA,
  mc_sd            = apply(est, 2, sd),
  coverage_largeT  = colMeans(stack("cov_largeT")),
  coverage_fixedT  = colMeans(stack("cov_fixedT")),
  coverage_boot    = if (bootstrap_B > 0) colMeans(stack("cov_boot")) else NA)
cat("\nInference: median standard errors, Monte Carlo s.d., and coverage of 95% intervals\n")
print(format(inf_tab, digits = 3), row.names = FALSE)

## Published values for comparison (i.i.d. normal errors, 1,000 replications)
published <- list(
  `3`  = list(mean = c(0.391, 0.099, 0.163), misclass = 9.50,
              se = rbind(c(0.035, 0.051, 0.068, 0.043), c(0.0093, 0.0132, 0.0156, 0.0137), c(0.013, 0.022, 0.030, 0.021)),
              cov = rbind(c(0.847, 0.965, 0.970), c(0.723, 0.883, 0.914), c(0.693, 0.917, 0.900))),
  `5`  = list(mean = c(0.262, 0.086, 0.117), misclass = 9.68,
              se = rbind(c(0.037, 0.068, 0.097, 0.058), c(0.0088, 0.0135, 0.0160, 0.0112), c(0.011, 0.022, 0.035, 0.022)),
              cov = rbind(c(0.848, 0.961, 0.973), c(0.788, 0.932, 0.943), c(0.710, 0.936, 0.928))),
  `10` = list(mean = c(0.286, 0.078, 0.110), misclass = 44.73,
              se = rbind(c(0.037, 0.048, 0.091, 0.059), c(0.0074, 0.0095, 0.0156, 0.0103), c(0.009, 0.012, 0.026, 0.015)),
              cov = rbind(c(0.798, 0.902, 0.992), c(0.841, 0.912, 0.996), c(0.783, 0.904, 0.986))))
if (dgp == "iid" && as.character(G0) %in% names(published)) {
  pb <- published[[as.character(G0)]]
  cat("\nPublished (Tables S.III and S.IV): GFE mean", pb$mean, "; misclassified", pb$misclass, "%\n")
  cat("  median s.e. (large-T, fixed-T, bootstrap) and MC s.d., by parameter:\n"); print(pb$se)
  cat("  coverage (large-T, fixed-T, bootstrap), by parameter:\n"); print(pb$cov)
}

tag <- sprintf("G%d_%s", G0, dgp)
write.csv(cbind(bias_tab, inf_tab[, -1]), sprintf("output/montecarlo_%s.csv", tag), row.names = FALSE)

## -----------------------------------------------------------------------------
## 5. Figure S.1: true and mean estimated group-specific time effects
## -----------------------------------------------------------------------------
alpha_mean <- Reduce(`+`, lapply(sims, `[[`, "alpha")) / n_rep
pdf(sprintf("output/figureS1_alpha_bias_%s.pdf", tag), width = 5, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))
matplot(1:T, truth$alpha, type = "l", lty = 1, lwd = 2, col = "black",
        ylim = range(truth$alpha, alpha_mean), xlab = "t", ylab = "group-specific time effects")
matlines(1:T, alpha_mean, lty = 2, lwd = 2, col = gray(0.4))
legend("bottomright", c("true", "mean of estimates"), lty = 1:2, lwd = 2,
       col = c("black", gray(0.4)), bty = "n")
invisible(dev.off())

## -----------------------------------------------------------------------------
## 6. Table S.VIII: BIC choice of the number of groups (optional)
## -----------------------------------------------------------------------------
## As in BIC_G3_MC.m, sigma2_hat is averaged over the simulated samples
## (sigma22 = mean(sigma2)) and the same value is used in every sample.
if (run_bic) {
  NT <- N * T
  Qm <- stack("bic_Q")
  s2 <- mean(sapply(sims, `[[`, "bic_sigma2"))
  bic <- sweep(Qm / NT, 2, s2 * log(NT) * (bic_G * T + N + K) / NT, "+")
  choice <- bic_G[apply(bic, 1, which.min)]
  cat("\nBIC choice of G across", n_rep, "samples (true G0 =", G0, "; sigma2_hat =", round(s2, 4), ", source", bic_sigma2, "):\n")
  tab <- table(factor(choice, levels = bic_G))
  print(tab)
  write.csv(data.frame(G = bic_G, samples = as.vector(tab)),
            file.path("output", sprintf("tableS8_bic_choice_G%d_%s.csv", G0, bic_sigma2)), row.names = FALSE)
}
