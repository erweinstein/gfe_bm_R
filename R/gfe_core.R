## =============================================================================
## R/gfe_core.R  --  The grouped fixed-effects (GFE) estimator
##
## An R reimplementation of the estimator and computational algorithms in
##
##   Bonhomme, Stephane and Elena Manresa (2015), "Grouped Patterns of
##   Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184.
##   doi:10.3982/ECTA11319
##
## The method, the algorithms, and the original MATLAB / Fortran / Stata code
## are by Stephane Bonhomme and Elena Manresa. This file translates their
## replication code (Bonhomme_Manresa_codes.zip, Econometric Society website)
## into R. Any errors introduced in translation are ours, not theirs.
##
## Original files this file follows:
##   Application/Data_paths.m          Algorithm 1 (iterative), MATLAB version
##   Application/Fortran codes/*.f90   Algorithms 1 and 2 (VNS), with and
##                                     without covariates, balanced/unbalanced
##   Application/Exact algorithm/heuristics_sub.f90   Algorithms 1-2 with
##                                     two covariates (used by the exact B&B)
##   Supplementary Appendix, Section S.1 (algorithm statements)
##
## -----------------------------------------------------------------------------
## The model (equation (1) in the paper):
##
##     y_it = x_it' theta + alpha_{g_i t} + v_it,   i = 1..N, t = 1..T,
##
## where g_i in {1,...,G} is unit i's (unknown) group and alpha_{gt} is a
## group-specific time effect. The GFE estimator solves
##
##     (theta_hat, alpha_hat, gamma_hat) =
##        argmin_{theta, alpha, gamma}  sum_i sum_t (y_it - x_it' theta - alpha_{g_i t})^2
##
## over all groupings gamma = {g_1, ..., g_N}. For a FIXED grouping the problem
## is ordinary least squares on group x time dummies; for FIXED (theta, alpha)
## each unit is simply assigned to the group whose time profile is closest.
## Algorithm 1 alternates between those two steps (it is k-means when K = 0).
##
## Data layout (identical to the original MATLAB/Fortran code):
##   y : vector of length N*T, stacked by unit, then time: row (i-1)*T + t
##   X : (N*T) x K matrix, same row ordering
##   d : optional 0/1 vector of length N*T marking observed (i,t) cells in an
##       unbalanced panel. As in GFE_code.do, missing cells carry y = x = 0 and
##       d = 0, and drop out of every sum.
## =============================================================================


## -----------------------------------------------------------------------------
## 0. Panel container
## -----------------------------------------------------------------------------

#' Build the panel object used by every GFE routine.
#'
#' @param y numeric vector, length N*T, stacked unit by unit.
#' @param X numeric matrix (N*T x K), or NULL when there are no covariates.
#' @param N,T number of units and of time periods (T = maximal span).
#' @param d optional 0/1 vector of observed cells (unbalanced panels).
#' @param unit_labels optional names for the N units (e.g. country codes).
#' @param time_labels optional labels for the T periods (e.g. years).
gfe_panel <- function(y, X = NULL, N, T, d = NULL,
                      unit_labels = NULL, time_labels = NULL) {
  stopifnot(length(y) == N * T)
  if (is.null(X)) X <- matrix(0, N * T, 0)
  X <- as.matrix(X)
  stopifnot(nrow(X) == N * T)
  if (is.null(d)) d <- rep(1, N * T)
  stopifnot(length(d) == N * T, all(d %in% c(0, 1)))

  ## Convention of GFE_code.do: unobserved cells are set to zero.
  y[d == 0] <- 0
  X[d == 0, ] <- 0
  if (anyNA(y) || anyNA(X)) stop("NA in an observed cell; set d = 0 there.")

  K <- ncol(X)

  ## Unit-by-unit arrays used for the sufficient-statistics representation
  ## of the objective (see gfe_stats below): Z[t, , i] = d_it * (y_it, x_it').
  Z <- array(0, dim = c(T, K + 1, N))
  Z[, 1, ] <- matrix(y, T, N)
  if (K > 0) for (k in seq_len(K)) Z[, k + 1, ] <- matrix(X[, k], T, N)

  structure(list(
    y = y, X = X, d = d, N = N, T = T, K = K,
    Dm = matrix(d, T, N),               # T x N matrix of observation indicators
    Z  = Z,
    W  = crossprod(cbind(y, X)),        # sum_it z_it z_it' (total cross products)
    unit = rep(seq_len(N), each = T),   # unit index of each row
    time = rep(seq_len(T), times = N),  # period index of each row
    unit_labels = if (is.null(unit_labels)) as.character(seq_len(N)) else unit_labels,
    time_labels = if (is.null(time_labels)) seq_len(T) else time_labels,
    balanced = all(d == 1)
  ), class = "gfe_panel")
}

#' Build a gfe_panel from a long data frame.
#'
#' Rows with a missing outcome or covariate become unobserved cells (d = 0),
#' as in the "AUTOMATIC CODE" part of GFE_code.do. Units are sorted by `id`.
gfe_panel_from_df <- function(df, id, time, y, x = character(0)) {
  ids <- sort(unique(df[[id]]))
  times <- sort(unique(df[[time]]))
  N <- length(ids); T <- length(times)
  full <- expand.grid(time = times, id = ids)[, c("id", "time")]
  key_full <- paste(full$id, full$time)
  key_df <- paste(df[[id]], df[[time]])
  if (anyDuplicated(key_df)) stop("duplicate (id, time) rows")
  m <- match(key_full, key_df)
  yy <- df[[y]][m]
  XX <- if (length(x)) as.matrix(df[m, x, drop = FALSE]) else NULL
  obs <- !is.na(m) & !is.na(yy)
  if (length(x)) obs <- obs & stats::complete.cases(XX)
  yy[!obs] <- 0
  if (length(x)) { XX[!obs, ] <- 0; colnames(XX) <- x }
  gfe_panel(yy, XX, N, T, d = as.numeric(obs),
            unit_labels = as.character(ids), time_labels = times)
}


## -----------------------------------------------------------------------------
## 1. Update step: OLS given a grouping
## -----------------------------------------------------------------------------

#' Least-squares estimates of (theta, alpha) for a given grouping.
#'
#' This is the "update" step of Algorithm 1. With group x time dummies, OLS is
#' equivalent to within-(group, time) demeaning, exactly as in the loops that
#' build `Ybar_gt` and `Xbar_gt` in GFE_estimates.m:
#'
#'   theta_hat = [sum_it (x_it - xbar_{g_i t})(x_it - xbar_{g_i t})']^{-1}
#'               sum_it (x_it - xbar_{g_i t})(y_it - ybar_{g_i t})
#'   alpha_hat_gt = ybar_gt - xbar_gt' theta_hat
#'
#' @return list(theta, alpha (T x G, NA where a group-period cell is empty),
#'   objective, residuals, n_gt (T x G cell counts), ok)
gfe_update <- function(panel, groups, G) {
  N <- panel$N; T <- panel$T; K <- panel$K; d <- panel$d
  gt <- (groups[panel$unit] - 1L) * T + panel$time      # cell index 1..G*T
  ## Group-time sums of (d, d*y, d*x); every cell is represented, even if empty.
  sums <- matrix(0, G * T, 2 + K)
  s <- rowsum(cbind(d, panel$y, panel$X), gt)
  sums[as.integer(rownames(s)), ] <- s
  n_gt <- sums[, 1]
  ok <- all(n_gt[gt[d == 1]] > 0)                      # always TRUE in practice
  denom <- ifelse(n_gt > 0, n_gt, NA)
  ybar <- sums[, 2] / denom
  xbar <- if (K > 0) sums[, 2 + seq_len(K), drop = FALSE] / denom else matrix(0, G * T, 0)

  ## Demeaned data (unobserved cells contribute zero).
  yd <- d * (panel$y - ybar[gt])
  yd[d == 0] <- 0
  if (K > 0) {
    Xd <- d * (panel$X - xbar[gt, , drop = FALSE])
    Xd[d == 0, ] <- 0
    XtX <- crossprod(Xd)
    theta <- tryCatch(drop(solve(XtX, crossprod(Xd, yd))), error = function(e) rep(NA_real_, K))
    if (anyNA(theta)) ok <- FALSE
    e <- yd - drop(Xd %*% theta)
    alpha <- ybar - drop(xbar %*% theta)
  } else {
    theta <- numeric(0)
    e <- yd
    alpha <- ybar
  }
  names(theta) <- colnames(panel$X)
  list(theta = theta,
       alpha = matrix(alpha, T, G),      # column g = time profile of group g
       objective = sum(e^2),
       residuals = e,
       n_gt = matrix(n_gt, T, G),
       ok = ok)
}


## -----------------------------------------------------------------------------
## 2. Assignment step: closest group profile
## -----------------------------------------------------------------------------

#' Distance of every unit to every group profile, for given (theta, alpha):
#'   dist[i, g] = sum_t d_it (y_it - x_it' theta - alpha_gt)^2.
gfe_distances <- function(panel, theta, alpha) {
  r <- panel$y
  if (panel$K > 0) r <- r - drop(panel$X %*% theta)
  R <- matrix(r, panel$T, panel$N)                   # column i = unit i's residual path
  G <- ncol(alpha)
  dist <- matrix(0, panel$N, G)
  for (g in seq_len(G)) {
    a <- alpha[, g]
    if (anyNA(a)) {
      ## A group with an empty period cannot host units observed in that period.
      bad <- is.na(a); a[bad] <- 0
      dist[, g] <- colSums(panel$Dm * (R - a)^2)
      dist[colSums(panel$Dm[bad, , drop = FALSE]) > 0, g] <- Inf
    } else {
      dist[, g] <- colSums(panel$Dm * (R - a)^2)
    }
  }
  dist
}

#' Assignment step of Algorithm 1, with the empty-group repair of the
#' Fortran code: if a group receives no unit, the unit farthest from its
#' assigned profile (taken from a group with more than one member) is moved
#' into it.
gfe_assign <- function(panel, theta, alpha) {
  dist <- gfe_distances(panel, theta, alpha)
  G <- ncol(alpha)
  groups <- max.col(-dist, ties.method = "first")    # Fortran minloc: first minimum
  best <- dist[cbind(seq_len(panel$N), groups)]
  for (g in seq_len(G)) {
    if (!any(groups == g)) {
      size <- tabulate(groups, G)
      cand <- which(size[groups] > 1)
      i <- cand[which.max(best[cand])]
      groups[i] <- g
      best[i] <- 0
    }
  }
  groups
}


## -----------------------------------------------------------------------------
## 3. Objective through sufficient statistics (used by the local search)
## -----------------------------------------------------------------------------
##
## For a grouping gamma, let z_it = (y_it, x_it')', S_gt = sum_{i in g} d_it z_it,
## n_gt = sum_{i in g} d_it, and
##
##     M(gamma) = sum_it d_it z_it z_it' - sum_{g,t} S_gt S_gt' / n_gt .
##
## M is the within-(group, time) cross-product matrix of (y, x). Concentrating
## out theta,
##
##     Q(gamma) = M_yy - M_yx M_xx^{-1} M_xy ,
##
## which equals the GFE objective at the best (theta, alpha) for gamma. Moving
## one unit changes only two groups' S_g, so the effect of a move on M can be
## computed in O(T (K+1)^2) operations. This is the same "objective as a
## function of the partition only" idea behind the matrices AMAT, BMAT1-2, and
## CMAT1-3 in Exact algorithm/Brusco_cov.f90.

gfe_stats <- function(panel, groups, G) {
  T <- panel$T; K1 <- panel$K + 1
  S <- array(0, dim = c(T, K1, G))
  n <- matrix(0, T, G)
  for (g in seq_len(G)) {
    idx <- which(groups == g)
    if (length(idx)) {
      S[, , g] <- rowSums(panel$Z[, , idx, drop = FALSE], dims = 2)
      n[, g] <- rowSums(panel$Dm[, idx, drop = FALSE])
    }
  }
  B <- vapply(seq_len(G), function(g) .between(S[, , g], n[, g]), numeric(K1 * K1))
  list(S = S, n = n, B = array(B, dim = c(K1, K1, G)))
}

## Between-cell cross products of one group: sum_t S_gt S_gt' / n_gt.
.between <- function(Sg, ng) {
  Sg <- matrix(Sg, ncol = length(Sg) / length(ng))
  ok <- ng > 0
  if (!any(ok)) return(numeric(ncol(Sg)^2))
  as.vector(crossprod(Sg[ok, , drop = FALSE] / sqrt(ng[ok])))
}

## Concentrated objective from the within cross-product matrix M.
.obj_from_M <- function(M) {
  K <- nrow(M) - 1
  if (K == 0) return(M[1, 1])
  if (K == 1) {                                  # closed forms for speed (K = 1, 2)
    if (M[2, 2] <= 1e-12 * abs(M[1, 1])) return(Inf)
    return(M[1, 1] - M[2, 1]^2 / M[2, 2])
  }
  if (K == 2) {
    a <- M[2, 2]; b <- M[2, 3]; c <- M[3, 3]; det <- a * c - b * b
    if (det <= 1e-12 * a * c) return(Inf)
    u <- M[2, 1]; v <- M[3, 1]
    return(M[1, 1] - (c * u * u - 2 * b * u * v + a * v * v) / det)
  }
  Mxx <- M[-1, -1, drop = FALSE]
  if (rcond(Mxx) < 1e-12) return(Inf)
  M[1, 1] - drop(crossprod(M[-1, 1], solve(Mxx, M[-1, 1])))
}

#' GFE objective Q(gamma) for a given grouping (theta and alpha concentrated out).
gfe_objective <- function(panel, groups, G) {
  st <- gfe_stats(panel, groups, G)
  .obj_from_M(panel$W - rowSums(st$B, dims = 2))
}


## -----------------------------------------------------------------------------
## 4. Algorithm 1: the simple iterative algorithm
## -----------------------------------------------------------------------------

#' Run Algorithm 1 to convergence from one starting value.
#'
#' Steps (paper, Section 2; Supplementary Appendix S.1):
#'   assignment: g_i = argmin_g sum_t (y_it - x_it' theta - alpha_gt)^2
#'   update:     (theta, alpha) = OLS of y on x and group x time dummies
#' until the grouping no longer changes (equivalently, until the parameter
#' change `deltapar` in the original code is exactly zero).
#'
#' @param theta0,alpha0 starting values; alternatively `groups0`, a starting
#'   grouping (used inside VNS, where Algorithm 1 restarts from a partition).
.alg1_run <- function(panel, G, theta0 = NULL, alpha0 = NULL, groups0 = NULL,
                      max_iter = 1000) {
  if (!is.null(groups0)) {
    groups <- groups0
    upd <- gfe_update(panel, groups, G)
    if (!upd$ok) return(list(objective = Inf, groups = groups))
    theta <- upd$theta; alpha <- upd$alpha
  } else {
    theta <- theta0; alpha <- alpha0; groups <- rep(0L, panel$N)
  }
  for (iter in seq_len(max_iter)) {
    new_groups <- gfe_assign(panel, theta, alpha)
    if (identical(new_groups, groups)) break          # deltapar == 0
    groups <- new_groups
    upd <- gfe_update(panel, groups, G)
    if (!upd$ok) return(list(objective = Inf, groups = groups))  # empty cell: restart
    theta <- upd$theta; alpha <- upd$alpha
  }
  upd <- gfe_update(panel, groups, G)
  list(theta = upd$theta, alpha = upd$alpha, groups = groups,
       objective = upd$objective, iterations = iter)
}

#' Random starting value, as in Data_paths.m and heuristics_sub.f90:
#'   1. draw theta0 ~ N(0, I_K);
#'   2. draw G distinct units at random and set alpha0_gt = y_{i_g t} - x_{i_g t}' theta0.
.random_start <- function(panel, G, theta_sd = 1) {
  theta0 <- stats::rnorm(panel$K, sd = theta_sd)
  r <- panel$y
  if (panel$K > 0) r <- r - drop(panel$X %*% theta0)
  R <- matrix(r, panel$T, panel$N)
  centers <- sample.int(panel$N, G)
  list(theta = theta0, alpha = R[, centers, drop = FALSE])
}

#' Algorithm 1 with many random starting values.
#'
#' Returns the solution with the lowest objective, together with the vector of
#' objectives reached from every start. Inspecting that vector is the check the
#' authors recommend in GFE_code.do ("Inspection of this file provides an idea
#' of the numerical stability of the solution"): a reliable solution is reached
#' from many different starting values.
#'
#' @param n_starts number of starting values (Ns; 1,000 in the paper's tables).
gfe_alg1 <- function(panel, G, n_starts = 1000, theta_sd = 1, max_iter = 1000,
                     verbose = FALSE) {
  best <- list(objective = Inf)
  objs <- numeric(n_starts)
  for (s in seq_len(n_starts)) {
    st <- .random_start(panel, G, theta_sd)
    fit <- .alg1_run(panel, G, theta0 = st$theta, alpha0 = st$alpha, max_iter = max_iter)
    objs[s] <- fit$objective
    if (fit$objective < best$objective) best <- fit
    if (verbose && s %% 100 == 0) message("Algorithm 1: start ", s, ", best objective ", signif(best$objective, 8))
  }
  .as_gfe(panel, G, best, objs, algorithm = "Algorithm 1 (iterative)",
          tuning = list(n_starts = n_starts))
}


## -----------------------------------------------------------------------------
## 5. Algorithm 2: Variable Neighborhood Search (VNS)
## -----------------------------------------------------------------------------

#' Local search (Step 5 of Algorithm 2).
#'
#' Systematically checks every reassignment of a unit i to another group g and
#' keeps it whenever the objective decreases, until no single reassignment
#' improves the objective. The resulting grouping is a "local minimum" in the
#' sense of the Supplementary Appendix: no single unit can be moved profitably,
#' which Algorithm 1 alone does not guarantee.
#'
#' For each unit we move it to the best improving group (best-improvement per
#' unit). The original Fortran tries groups in order and accepts any
#' improvement; both satisfy the stopping rule above.
gfe_local_search <- function(panel, groups, G, tol = 1e-10) {
  st <- gfe_stats(panel, groups, G)
  M <- panel$W - rowSums(st$B, dims = 2)
  obj <- .obj_from_M(M)
  size <- tabulate(groups, G)
  K1 <- panel$K + 1
  repeat {
    improved <- FALSE
    for (i in seq_len(panel$N)) {
      g0 <- groups[i]
      if (size[g0] <= 1) next                          # never empty a group
      zi <- panel$Z[, , i]; di <- panel$Dm[, i]
      Sg0 <- st$S[, , g0] - zi; ng0 <- st$n[, g0] - di
      Bg0_new <- .between(Sg0, ng0)
      best_obj <- obj; best_g <- g0
      best_parts <- NULL
      for (g in seq_len(G)[-g0]) {
        Sg <- st$S[, , g] + zi; ng <- st$n[, g] + di
        Bg_new <- .between(Sg, ng)
        Mnew <- M + st$B[, , g0] + st$B[, , g] -
          matrix(Bg0_new, K1, K1) - matrix(Bg_new, K1, K1)
        o <- .obj_from_M(Mnew)
        if (o < best_obj - tol * max(1, abs(best_obj))) {
          best_obj <- o; best_g <- g
          best_parts <- list(Mnew = Mnew, Sg = Sg, ng = ng, Bg = Bg_new)
        }
      }
      if (best_g != g0) {
        g <- best_g
        st$S[, , g0] <- Sg0; st$n[, g0] <- ng0; st$B[, , g0] <- Bg0_new
        st$S[, , g] <- best_parts$Sg; st$n[, g] <- best_parts$ng; st$B[, , g] <- best_parts$Bg
        M <- best_parts$Mnew; obj <- best_obj
        groups[i] <- g
        size[g0] <- size[g0] - 1; size[g] <- size[g] + 1
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  list(groups = groups, objective = obj)
}

#' Neighborhood jump (Step 3 of Algorithm 2): relocate n randomly selected
#' units to n randomly selected (different) groups, never emptying a group.
.vns_jump <- function(groups, G, n) {
  N <- length(groups)
  for (k in seq_len(n)) {
    size <- tabulate(groups, G)
    movable <- which(size[groups] > 1)
    if (!length(movable)) break
    i <- movable[sample.int(length(movable), 1)]
    groups[i] <- sample(setdiff(seq_len(G), groups[i]), 1)
  }
  groups
}

#' Algorithm 2: Variable Neighborhood Search.
#'
#' Supplementary Appendix, Section S.1 (extending Pacheco and Valencia (2003)
#' and Brusco and Steinley (2007) to allow for covariates):
#'
#'   1. From a starting value, obtain an initial grouping gamma*.
#'      (As in the Fortran code, we run Algorithm 1 to convergence from a
#'      random start.)  Set j = 0.
#'   2. Set n = 1.
#'   3. Neighborhood jump: relocate n random units to n random groups (gamma').
#'      One update step gives (theta', alpha').
#'   4. Apply Algorithm 1 from (theta', alpha').
#'   5. Local search: reassign single units while the objective decreases.
#'   6. If the objective improves on gamma*, set gamma* to the new grouping and
#'      go to 2; otherwise set n = n + 1.
#'   7. If n <= neighmax go to 3; otherwise go to 8.
#'   8. Set j = j + 1. If j > itermax stop; otherwise go to 2.
#'
#' Tuning parameters (Ns; neighmax; itermax): the paper uses (10; 10; 10) in
#' the main estimates and (5; 10; 5) for bootstrap and Monte Carlo exercises.
#' Rule of thumb from the authors: check that different starting values yield
#' the exact same solution.
gfe_vns <- function(panel, G, n_starts = 10, neighmax = 10, itermax = 10,
                    theta_sd = 1, max_iter = 1000, verbose = FALSE) {
  tol <- 1e-10
  best <- list(objective = Inf)
  objs <- numeric(n_starts)
  for (s in seq_len(n_starts)) {
    ## Step 1: initial grouping from Algorithm 1 (one random start).
    st <- .random_start(panel, G, theta_sd)
    cur <- .alg1_run(panel, G, theta0 = st$theta, alpha0 = st$alpha, max_iter = max_iter)
    if (!is.finite(cur$objective)) { objs[s] <- Inf; next }
    cur <- gfe_local_search(panel, cur$groups, G)
    j <- 0
    while (j <= itermax) {
      n <- 1                                                     # Step 2
      while (n <= neighmax) {                                    # Step 7
        prop <- .vns_jump(cur$groups, G, n)                      # Step 3
        a1 <- .alg1_run(panel, G, groups0 = prop, max_iter = max_iter)  # Step 4
        if (!is.finite(a1$objective)) { n <- n + 1; next }
        ls <- gfe_local_search(panel, a1$groups, G)              # Step 5
        if (ls$objective < cur$objective - tol * max(1, abs(cur$objective))) {
          cur <- ls                                              # Step 6: improvement
          n <- 1
        } else {
          n <- n + 1
        }
      }
      j <- j + 1                                                 # Step 8
    }
    objs[s] <- cur$objective
    if (cur$objective < best$objective) best <- cur
    if (verbose) message("VNS: start ", s, "/", n_starts, ", objective ", signif(cur$objective, 10),
                         ", best ", signif(best$objective, 10))
  }
  upd <- gfe_update(panel, best$groups, G)
  best <- c(upd[c("theta", "alpha", "objective")], list(groups = best$groups))
  .as_gfe(panel, G, best, objs, algorithm = "Algorithm 2 (VNS)",
          tuning = list(n_starts = n_starts, neighmax = neighmax, itermax = itermax))
}


## -----------------------------------------------------------------------------
## 6. Main entry point and the "gfe" object
## -----------------------------------------------------------------------------

#' Grouped fixed-effects estimation.
#'
#' @param panel a gfe_panel.
#' @param G number of groups.
#' @param algorithm "vns" (Algorithm 2, default) or "alg1" (Algorithm 1).
#' @param groups optional: a grouping to evaluate directly (no search), e.g.
#'   the published groupings in data/bm_published_groupings.csv.
#' @param ... tuning parameters passed to gfe_vns() or gfe_alg1().
gfe <- function(panel, G, algorithm = c("vns", "alg1"), groups = NULL, ...) {
  algorithm <- match.arg(algorithm)
  if (G == 1) groups <- rep(1L, panel$N)
  if (!is.null(groups)) {
    groups <- as.integer(groups)
    upd <- gfe_update(panel, groups, G)
    return(.as_gfe(panel, G, c(upd[c("theta", "alpha", "objective")], list(groups = groups)),
                   upd$objective, algorithm = "given grouping", tuning = list()))
  }
  if (algorithm == "vns") gfe_vns(panel, G, ...) else gfe_alg1(panel, G, ...)
}

.as_gfe <- function(panel, G, best, objs, algorithm, tuning) {
  groups <- as.integer(best$groups)
  upd <- gfe_update(panel, groups, G)
  rownames(upd$alpha) <- panel$time_labels
  colnames(upd$alpha) <- paste0("g", seq_len(G))
  names(groups) <- panel$unit_labels
  structure(list(
    theta = upd$theta, alpha = upd$alpha, groups = groups,
    objective = upd$objective, residuals = upd$residuals,
    G = G, N = panel$N, T = panel$T, K = panel$K,
    group_sizes = tabulate(groups, G),
    objectives_by_start = objs, algorithm = algorithm, tuning = tuning,
    panel = panel
  ), class = "gfe")
}

print.gfe <- function(x, digits = 4, ...) {
  cat("Grouped fixed-effects estimator (Bonhomme and Manresa, 2015)\n")
  cat("  N =", x$N, " T =", x$T, " K =", x$K, " G =", x$G, "\n")
  cat("  computed with:", x$algorithm, "\n")
  cat("  objective:", format(x$objective, digits = 10), "\n")
  cat("  group sizes:", x$group_sizes, "\n")
  if (x$K > 0) { cat("  theta:\n"); print(round(x$theta, digits)) }
  ob <- x$objectives_by_start
  if (length(ob) > 1) {
    share <- mean(abs(ob - min(ob)) < 1e-8 * max(1, abs(min(ob))))
    cat(sprintf("  starts reaching the best objective: %d of %d (%.0f%%)\n",
                round(share * length(ob)), length(ob), 100 * share))
  }
  invisible(x)
}

#' Relabel groups 1..G in order of a chosen criterion (by default, by the
#' first-period value of alpha), so solutions from different runs can be
#' compared. Group labels are arbitrary in the GFE problem.
gfe_relabel <- function(fit, order_by = fit$alpha[1, ]) {
  o <- order(order_by)
  map <- integer(fit$G); map[o] <- seq_len(fit$G)
  fit$groups[] <- map[fit$groups]
  fit$alpha <- fit$alpha[, o, drop = FALSE]
  colnames(fit$alpha) <- paste0("g", seq_len(fit$G))
  fit$group_sizes <- tabulate(fit$groups, fit$G)
  fit
}
