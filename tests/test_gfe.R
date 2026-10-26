## =============================================================================
## tests/test_gfe.R  --  Unit tests for the R reimplementation of
## Bonhomme and Manresa (2015), "Grouped Patterns of Heterogeneity in Panel
## Data," Econometrica 83(3), 1147-1184.
##
## Each test checks one building block against an independent computation
## (brute force, lm(), a direct matrix formula). No packages needed:
##   Rscript tests/test_gfe.R          (from the project folder)
## Validation against the published numbers is in scripts/01_validate_*.R.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

source("gfe.R")
n_pass <- 0
check <- function(desc, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", desc, call. = FALSE)
  n_pass <<- n_pass + 1
  cat("ok  ", desc, "\n")
}

## A small simulated GFE panel with K = 2 covariates.
sim_panel <- function(N, T, G, K = 2, seed = 1, unbalanced = FALSE) {
  set.seed(seed)
  g0 <- sample(1:G, N, replace = TRUE)
  a0 <- matrix(rnorm(T * G, sd = 1.5), T, G)
  X  <- matrix(rnorm(N * T * K), N * T, K, dimnames = list(NULL, paste0("x", 1:K)))
  y  <- drop(X %*% seq(0.5, -0.3, length.out = K)) + a0[cbind(rep(1:T, N), rep(g0, each = T))] + rnorm(N * T)
  d  <- if (unbalanced) as.numeric(runif(N * T) > 0.15) else NULL
  gfe_panel(y, X, N, T, d = d)
}

## 1. The sufficient-statistics objective equals the OLS objective -------------
p <- sim_panel(40, 6, 3)
for (s in 1:5) {
  g <- sample(1:3, p$N, replace = TRUE)
  check(sprintf("gfe_objective() = SSR of the update step (draw %d)", s),
        abs(gfe_objective(p, g, 3) - gfe_update(p, g, 3)$objective) < 1e-8)
}

## 2. The update step is OLS on group x time dummies (unbalanced panel) --------
pu <- sim_panel(30, 5, 3, unbalanced = TRUE)
g  <- sample(1:3, pu$N, replace = TRUE)
up <- gfe_update(pu, g, 3)
obs <- pu$d == 1
cell <- factor(paste(g[pu$unit], pu$time))[obs]
ref <- lm(pu$y[obs] ~ pu$X[obs, ] + cell - 1)
check("gfe_update() theta = lm() with group x time dummies (unbalanced)",
      max(abs(up$theta - coef(ref)[1:2])) < 1e-8)
check("gfe_update() objective = lm() residual sum of squares (unbalanced)",
      abs(up$objective - sum(resid(ref)^2)) < 1e-8)

## 3. Branch and bound = brute force; VNS finds the same optimum ---------------
for (G in 2:3) {
  pt <- sim_panel(9, 4, G, seed = 10 + G)
  bf <- gfe_bruteforce(pt, G)
  ex <- gfe_exact_bb(pt, G, engine = "R")
  vn <- gfe(pt, G, algorithm = "vns", n_starts = 5, neighmax = 5, itermax = 5)
  check(sprintf("exact branch and bound (R engine) = brute force, G = %d", G),
        abs(ex$objective - bf$objective) < 1e-8)
  check(sprintf("VNS reaches the global minimum on a toy problem, G = %d", G),
        abs(vn$objective - bf$objective) < 1e-8)
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    exc <- try(gfe_exact_bb(pt, G, engine = "cpp"), silent = TRUE)
    if (inherits(exc, "gfe"))
      check(sprintf("C++ engine = R engine (objective and node count), G = %d", G),
            abs(exc$objective - ex$objective) < 1e-10 && exc$nodes == ex$nodes)
  }
}

## 4. Local search never increases the objective, and ends at a local minimum -
p <- sim_panel(40, 6, 3, seed = 3)
g <- sample(1:3, p$N, replace = TRUE)
ls <- gfe_local_search(p, g, 3)
check("local search does not increase the objective", ls$objective <= gfe_objective(p, g, 3) + 1e-10)
no_better <- TRUE
for (i in seq_len(p$N)) for (h in setdiff(1:3, ls$groups[i])) {
  g2 <- ls$groups; g2[i] <- h
  if (all(tabulate(g2, 3) > 0) && gfe_objective(p, g2, 3) < ls$objective - 1e-9) no_better <- FALSE
}
check("no single-unit move improves on the local-search solution", no_better)

## 5. Hungarian algorithm = minimum over all G! relabelings --------------------
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(r) c(v[i], r))))
set.seed(7)
for (s in 1:5) {
  G <- 5
  truth <- sample(1:G, 60, replace = TRUE)
  est <- ifelse(runif(60) < 0.3, sample(1:G, 60, replace = TRUE), sample(G)[truth])
  brute <- min(sapply(perms(1:G), function(pm) mean(pm[est] != truth)))
  check(sprintf("match_groups() = minimum over all 120 permutations (draw %d)", s),
        abs(match_groups(est, truth, G)$misclassification - brute) < 1e-12)
}

## 6. Large-T clustered variance = direct sandwich formula --------------------
p <- sim_panel(50, 5, 2, seed = 4)
fit <- gfe(p, 2, algorithm = "alg1", n_starts = 50)
D <- matrix(0, p$N * p$T, 2 * p$T)
D[cbind(seq_len(p$N * p$T), (fit$groups[p$unit] - 1) * p$T + p$time)] <- 1
Xt <- cbind(D, p$X)
e <- resid(lm(p$y ~ Xt - 1))
B <- solve(crossprod(Xt))
M <- Reduce(`+`, lapply(seq_len(p$N), function(i) {
  r <- (i - 1) * p$T + seq_len(p$T)
  crossprod(Xt[r, ], e[r]) %*% crossprod(e[r], Xt[r, ])
}))
V <- B %*% M %*% B * (p$N * p$T) / (p$N * p$T - 2 * p$T - p$K)
check("gfe_vcov_largeT() = direct clustered sandwich",
      max(abs(gfe_vcov_largeT(fit) - V)) < 1e-10)

## 7. Within transformation: G = 1 GFE in deviations = two-way fixed effects --
p <- sim_panel(30, 6, 2, seed = 5)
fe <- gfe(gfe_within(p), 1)
ref <- lm(p$y ~ p$X + factor(p$unit) + factor(p$time))
check("GFE (G = 1) in deviations from unit means = two-way fixed effects",
      max(abs(fe$theta - coef(ref)[2:3])) < 1e-8)

## 8. Two-layer model: constrained LS = unrestricted LS on its own parametrization
p <- sim_panel(36, 5, 3, seed = 6)
H <- c(2, 2, 1)
sg <- sample(seq_len(sum(H)), p$N, replace = TRUE)
while (length(unique(sg)) < sum(H)) sg <- sample(seq_len(sum(H)), p$N, replace = TRUE)
tl <- gfe_twolayer(p, H, groups = sg)
l1 <- rep(seq_along(H), H)[sg]
ref <- lm(p$y ~ p$X + factor(paste(l1[p$unit], p$time)) + factor(sg[p$unit]))
check("two-layer constrained LS = OLS on (layer-1 x time) and subgroup dummies",
      max(abs(tl$theta - coef(ref)[2:3])) < 1e-8 && abs(tl$objective - sum(resid(ref)^2)) < 1e-8)

cat("\nAll", n_pass, "tests passed.\n")
