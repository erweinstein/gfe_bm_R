## =============================================================================
## R/gfe_utils.R  --  Helpers: group matching, misclassification, figures,
##                    and small linear-algebra / regression utilities
##
## Part of an R reimplementation of Bonhomme and Manresa (2015), "Grouped
## Patterns of Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================


## -----------------------------------------------------------------------------
## Matching estimated groups to true groups (label invariance)
## -----------------------------------------------------------------------------
##
## Group labels are not identified: any permutation of {1,...,G} gives the same
## fit. To measure misclassification, footnote 21 of the Supplementary
## Appendix takes, in each sample, "the labeling that yields the minimum
## amount of misclassification across all G! permutations" (and, for G = 10,
## the minimum over 500,000 random permutations). That minimum is a linear
## assignment problem, solved exactly here by the Hungarian algorithm.

#' Hungarian algorithm (Kuhn-Munkres) for a square cost matrix.
#' @return integer vector `assign` with assign[i] = column matched to row i,
#'   minimizing sum(cost[cbind(i, assign[i])]).
hungarian <- function(cost) {
  n <- nrow(cost); stopifnot(ncol(cost) == n)
  INF <- .Machine$double.xmax / 4
  u <- numeric(n + 1); v <- numeric(n + 1)
  p <- integer(n + 1); way <- integer(n + 1)          # 1-based, index 0 = dummy
  for (i in seq_len(n)) {
    p[1] <- i
    j0 <- 0L
    minv <- rep(INF, n + 1); used <- rep(FALSE, n + 1)
    repeat {
      used[j0 + 1] <- TRUE
      i0 <- p[j0 + 1]; delta <- INF; j1 <- 0L
      for (j in seq_len(n)) if (!used[j + 1]) {
        cur <- cost[i0, j] - u[i0 + 1] - v[j + 1]
        if (cur < minv[j + 1]) { minv[j + 1] <- cur; way[j + 1] <- j0 }
        if (minv[j + 1] < delta) { delta <- minv[j + 1]; j1 <- j }
      }
      for (j in 0:n) {
        if (used[j + 1]) { u[p[j + 1] + 1] <- u[p[j + 1] + 1] + delta; v[j + 1] <- v[j + 1] - delta }
        else minv[j + 1] <- minv[j + 1] - delta
      }
      j0 <- j1
      if (p[j0 + 1] == 0) break
    }
    repeat {
      j1 <- way[j0 + 1]; p[j0 + 1] <- p[j1 + 1]; j0 <- j1
      if (j0 == 0) break
    }
  }
  assign <- integer(n)
  for (j in seq_len(n)) if (p[j + 1] > 0) assign[p[j + 1]] <- j
  assign
}

#' Relabel estimated groups to best match true groups.
#' @return list(groups = relabeled estimate, map = old -> new labels,
#'   misclassification = share of units with relabeled group != true group).
match_groups <- function(est, true, G = max(est, true)) {
  overlap <- table(factor(est, levels = seq_len(G)), factor(true, levels = seq_len(G)))
  a <- hungarian(-unclass(overlap))
  relabeled <- a[est]
  list(groups = relabeled, map = a, misclassification = mean(relabeled != true))
}


## -----------------------------------------------------------------------------
## Group averages of the data (MYbar_gt and MXbar_gt in GFE_estimates.m)
## -----------------------------------------------------------------------------

#' T x G matrix of group-by-period means of a variable stacked like y.
group_means <- function(fit, v) {
  p <- fit$panel
  M <- matrix(v, p$T, p$N)
  D <- p$Dm
  out <- vapply(seq_len(fit$G), function(g) {
    idx <- fit$groups == g
    rowSums((M * D)[, idx, drop = FALSE]) / rowSums(D[, idx, drop = FALSE])
  }, numeric(p$T))
  out <- matrix(out, p$T, fit$G, dimnames = list(p$time_labels, paste0("g", seq_len(fit$G))))
  out
}


## -----------------------------------------------------------------------------
## Figures in the style of the paper (black, white, and grays; base graphics)
## -----------------------------------------------------------------------------

#' Plot T x G profiles (group time effects, or group averages of a variable),
#' in the grayscale line style of GFE_estimates.m (thick black and gray lines,
#' solid / dashed / dash-dot).
plot_profiles <- function(M, years = as.numeric(rownames(M)), ylab = "time effects",
                          xlab = "years", ylim = range(M, na.rm = TRUE),
                          col = NULL, lty = NULL, lwd = 4, main = "", ...) {
  G <- ncol(M)
  if (is.null(col)) col <- grDevices::gray(c(0, 0.75, 0.4, 0.6, 0.2, 0.5, 0.3, 0.65, 0.1, 0.85)[(seq_len(G) - 1) %% 10 + 1])
  if (is.null(lty)) lty <- rep(c(1, 2, 4, 1, 2, 4), length.out = G)
  graphics::matplot(years, M, type = "l", lty = lty, lwd = lwd, col = col,
                    xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
  invisible(NULL)
}

#' Parameter estimates against the number of groups, with +/- 1.96 se bands
#' (the "Figure of parameter estimates" block in GFE_estimates.m).
plot_by_G <- function(G, est, se, ylab, ylim = NULL, ...) {
  lo <- est - 1.96 * se; hi <- est + 1.96 * se
  if (is.null(ylim)) ylim <- range(lo, hi, na.rm = TRUE)
  graphics::plot(G, est, pch = 1, cex = 1.6, ylim = ylim, xaxt = "n",
                 xlab = "number of groups G", ylab = ylab, ...)
  graphics::axis(1, at = G)
  graphics::lines(G, est, lwd = 4)
  graphics::lines(G, lo, lty = 3, lwd = 3)
  graphics::lines(G, hi, lty = 3, lwd = 3)
  invisible(NULL)
}


## -----------------------------------------------------------------------------
## Linear algebra and regression helpers
## -----------------------------------------------------------------------------

#' Moore-Penrose pseudo-inverse (MATLAB pinv), via the SVD.
pinv <- function(A, tol = NULL) {
  s <- svd(A)
  if (is.null(tol)) tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  if (!any(pos)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

#' OLS / 2SLS with standard errors clustered by `cluster`.
#'
#' 2SLS: y on X (exogenous and endogenous regressors), instruments Z (exogenous
#' regressors plus excluded instruments).
#'
#' @param small small-sample factor applied to the cluster-robust variance:
#'   "none" (Stata's ivregress, vce(cluster), without the `small` option;
#'   this reproduces the standard errors in Table S.XV), "clusters" (M/(M-1)),
#'   or "regress" (M/(M-1) * (n-1)/(n-k), as Stata's regress, vce(cluster)).
ols_iv_cluster <- function(y, X, Z = X, cluster, small = c("none", "clusters", "regress")) {
  small <- match.arg(small)
  X <- as.matrix(X); Z <- as.matrix(Z)
  Pz_X <- Z %*% solve(crossprod(Z), crossprod(Z, X))      # first-stage fitted X
  b <- drop(solve(crossprod(Pz_X, X), crossprod(Pz_X, y)))
  e <- drop(y - X %*% b)
  bread <- solve(crossprod(Pz_X))
  score <- rowsum(Pz_X * e, cluster)
  M <- nrow(score); n <- length(y); k <- ncol(X)
  adj <- switch(small, none = 1, clusters = M / (M - 1), regress = M / (M - 1) * (n - 1) / (n - k))
  V <- bread %*% crossprod(score) %*% bread * adj
  dimnames(V) <- list(colnames(X), colnames(X))
  names(b) <- colnames(X)
  list(coef = b, vcov = V, se = sqrt(diag(V)), residuals = e, n_clusters = M)
}
