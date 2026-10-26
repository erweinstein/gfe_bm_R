## =============================================================================
## scripts/03_exact_algorithm.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184,
## Supplementary Appendix S.1.2 and Table S.I.
##
## Purpose. The GFE objective is minimized over all G^N groupings, so the
## heuristics (Algorithms 1 and 2) come with no guarantee of reaching the
## global minimum. Bonhomme and Manresa check them against an EXACT solution
## computed by their covariate extension of Brusco's (2006) repetitive branch
## and bound (Exact algorithm/Brusco_cov.f90). This script repeats that check:
##
##   1. A small simulated problem, where the exact solution can also be found
##      by brute-force enumeration: brute force = R engine = C++ engine.
##   2. The first N = 30 countries of the application, G = 2, with the data
##      rounded to 4 decimals as in the authors' input files inputY.txt and
##      inputX.txt. Target (Exact algorithm/RESULTS2.DAT and outputtheta.txt):
##      objective 6.1591..., theta = (0.6070, 0.0689).
##   3. Optionally (run_N90 = TRUE) the full sample, N = 90, G = 2: target
##      objective 19.846 (Table S.I). This takes a long time.
##
## In each case we also run Algorithm 1 and VNS, and report whether they reach
## the exact optimum: that is the comparison made in Table S.I.
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

run_N90 <- FALSE        # TRUE: N = 90, G = 2 (hours, even with the C++ engine)
engine  <- "auto"       # "cpp" (needs Rcpp + Rtools on Windows) or "R"

## -----------------------------------------------------------------------------
## 1. A brute-force check on a toy problem (N = 10, T = 4, K = 2, G = 3)
## -----------------------------------------------------------------------------
set.seed(2015)
N <- 10; T <- 4; G <- 3
g0 <- sample(1:G, N, replace = TRUE)
a0 <- matrix(rnorm(T * G), T, G)
X  <- matrix(rnorm(N * T * 2), N * T, 2)
y  <- drop(X %*% c(0.5, -0.3)) + a0[cbind(rep(1:T, N), rep(g0, each = T))] + rnorm(N * T, sd = 0.7)
toy <- gfe_panel(y, X, N, T)

bf   <- gfe_bruteforce(toy, G)              # G^N = 59,049 groupings
exR  <- gfe_exact_bb(toy, G, engine = "R")
exC  <- if (engine != "R") try(gfe_exact_bb(toy, G, engine = engine), silent = TRUE) else NULL
cat("\nToy problem (N = 10, G = 3):\n")
cat(sprintf("  brute force          %.10f\n", bf$objective))
cat(sprintf("  branch & bound (R)   %.10f   nodes = %d\n", exR$objective, exR$nodes))
if (inherits(exC, "gfe"))
  cat(sprintf("  branch & bound (C++) %.10f   nodes = %d\n", exC$objective, exC$nodes))

## -----------------------------------------------------------------------------
## 2. N = 30, G = 2 (Table S.I, first row)
## -----------------------------------------------------------------------------
## The authors' input files hold the data with 4 decimals. With unrounded data
## the optimum is 6.15958; with the rounded data it is 6.15913, the published
## value. Either way the optimal grouping is the same.
dat <- read.csv("data/bm_final_data.csv", stringsAsFactors = FALSE)
r4  <- function(v) round(v, 4)
n30 <- 30; rows <- seq_len(n30 * 7)
p30 <- gfe_panel(r4(dat$dem[rows]), r4(as.matrix(dat[rows, c("lag_dem", "lag_income")])),
                 N = n30, T = 7, unit_labels = unique(dat$code)[1:n30])

set.seed(1)
t0 <- Sys.time()
ex30 <- gfe_exact_bb(p30, 2, engine = engine, verbose = FALSE)
t_exact <- Sys.time() - t0
a1_30  <- gfe(p30, 2, algorithm = "alg1", n_starts = 1000)
vns_30 <- gfe(p30, 2, algorithm = "vns", n_starts = 10, neighmax = 10, itermax = 10)

res <- data.frame(
  N = n30, G = 2,
  method    = c("exact (branch & bound)", "Algorithm 1 (1,000 starts)", "Algorithm 2, VNS (10; 10; 10)"),
  objective = c(ex30$objective, a1_30$objective, vns_30$objective),
  theta1    = c(ex30$theta[1], a1_30$theta[1], vns_30$theta[1]),
  theta2    = c(ex30$theta[2], a1_30$theta[2], vns_30$theta[2]),
  misclassified_vs_exact = c(0, match_groups(a1_30$groups, ex30$groups, 2)$misclassification * n30,
                             match_groups(vns_30$groups, ex30$groups, 2)$misclassification * n30))
cat("\nN = 30, G = 2 (exact solution in", format(t_exact, digits = 3), "; nodes =", ex30$nodes, ")\n")
print(format(res, digits = 8), row.names = FALSE)
cat("Published (RESULTS2.DAT / outputtheta.txt): objective 6.1591, theta = (0.6070, 0.0689)\n")

## -----------------------------------------------------------------------------
## 3. N = 90, G = 2 (Table S.I): optional
## -----------------------------------------------------------------------------
if (run_N90) {
  p90 <- load_bm_panel()
  set.seed(1)
  ex90 <- gfe_exact_bb(p90, 2, engine = engine, verbose = TRUE)
  vns90 <- gfe(p90, 2, algorithm = "vns", n_starts = 10, neighmax = 10, itermax = 10)
  res <- rbind(res, data.frame(
    N = 90, G = 2, method = c("exact (branch & bound)", "Algorithm 2, VNS (10; 10; 10)"),
    objective = c(ex90$objective, vns90$objective),
    theta1 = c(ex90$theta[1], vns90$theta[1]), theta2 = c(ex90$theta[2], vns90$theta[2]),
    misclassified_vs_exact = c(0, match_groups(vns90$groups, ex90$groups, 2)$misclassification * 90)))
  cat("Published (Table S.I): 19.846\n")
}

write.csv(res, "output/tableS1_exact_vs_heuristics.csv", row.names = FALSE)
