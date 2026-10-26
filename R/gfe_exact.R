## =============================================================================
## R/gfe_exact.R  --  Exact GFE solution by repetitive branch and bound
##
## R translation of Bonhomme and Manresa's extension of Brusco's (2006)
## repetitive branch-and-bound algorithm for minimum within-cluster sums of
## squares partitioning, modified to allow for covariates:
##
##   Original: Bonhomme_Manresa_codes/Application/Exact algorithm/Brusco_cov.f90
##             (with heuristics_sub.f90 for the starting heuristic and the
##              COMPUTE_MATRIX subroutine for the pairwise matrices)
##   Method:   Brusco, M. J. (2006), "A Repetitive Branch-and-Bound Procedure
##             for Minimum Within-Cluster Sums of Squares Partitioning,"
##             Psychometrika 71(2), 347-363.
##   Use:      Supplementary Appendix of Bonhomme and Manresa (2015), Section
##             S.1.2 and Table S.I ("exact" column, single star): the exact
##             solution certifies that the heuristics found the global minimum.
##
## The authors report that their implementation was practical for small
## problems (G = 2 with N <= 90, and G = 3 with N = 30). The number of nodes
## grows very quickly with N, so the pure-R stage loop below (kept as the
## readable reference) is practical only up to N of about 20. The same loop in
## C++ (src/gfe_exact_bb.cpp, engine = "cpp") plays the role of their Fortran.
##
## -----------------------------------------------------------------------------
## Key identity. For a grouping gamma with n_g units in group g, the
## within-group sum of squares of any T-vectors w_i equals
##
##   sum_g sum_{i in g} || w_i - wbar_g ||^2 = sum_g (1/n_g) sum_{i<j in g} || w_i - w_j ||^2 .
##
## Applying it to w_i = y_i - x_i theta and concentrating out theta gives the GFE
## objective as a function of PAIRWISE quantities only (balanced panels):
##
##   A_ij      = sum_t (y_it - y_jt)^2
##   B_ij[k]   = sum_t (y_it - y_jt)(x_itk - x_jtk)
##   C_ij[k,l] = sum_t (x_itk - x_jtk)(x_itl - x_jtl)
##
##   a = sum_g (1/n_g) sum_{i<j in g} A_ij,  b and C defined alike,
##   Q(gamma) = a - b' C^{-1} b .
##
## These are the matrices AMAT, BMAT1-2, CMAT1-3 read by Brusco_cov.f90.
##
## Why the bound is valid with covariates. Restricting a grouping to a subset
## of units can only lower the within sum of squares for every theta, and
## min_theta [f1(theta) + f2(theta)] >= min f1 + min f2. Hence, if the first P
## units are assigned, the objective of ANY completion is at least
##   Q(partial assignment of the first P units) + Q*(remaining units),
## where Q*(remaining) is the optimum for the remaining units on their own.
## The "repetitive" part of Brusco's method computes these Q*(.) by solving
## the problem for the last G+1, G+2, ..., N units in turn.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

#' Exact GFE solution by repetitive branch and bound (balanced panels).
#'
#' @param panel a balanced gfe_panel.
#' @param G number of groups (the authors used G = 2, and G = 3 for N = 30).
#' @param heuristic optional "gfe" fit used for the object ordering and the
#'   upper bounds. By default a VNS fit is computed.
#' @param stage_starts Algorithm 1 starts used to get an upper bound at each
#'   stage (Brusco_cov.f90 uses heuristics_simple with 1,000 starts).
#' @param engine "auto" (C++ if Rcpp is available, else R), "cpp", or "R".
#' @param verbose print the optimum of each stage, as the Fortran code does.
#' @return a "gfe" object (algorithm = "exact branch and bound"), plus
#'   `stage_optima` and `nodes` (number of branch-and-bound nodes visited).
gfe_exact_bb <- function(panel, G, heuristic = NULL, stage_starts = 50,
                         engine = c("auto", "cpp", "R"), verbose = FALSE) {
  if (!panel$balanced) stop("gfe_exact_bb() requires a balanced panel.")
  engine <- match.arg(engine)
  if (engine != "R") {
    ok <- .gfe_load_cpp()
    if (!ok && engine == "cpp") stop("C++ engine unavailable (needs Rcpp and a compiler).")
    engine <- if (ok) "cpp" else "R"
  }
  N <- panel$N; T <- panel$T; K <- panel$K
  t_start <- Sys.time()

  ## ---------------------------------------------------------------------------
  ## STEP A. A good heuristic solution (Brusco_cov.f90: CALL heuristics, i.e.
  ## VNS with (Nsim; neighmax; itermax) = (10; 10; 2)).
  ## ---------------------------------------------------------------------------
  if (is.null(heuristic)) heuristic <- gfe_vns(panel, G, n_starts = 10, neighmax = 10, itermax = 2)
  hg <- as.integer(heuristic$groups)

  ## ---------------------------------------------------------------------------
  ## STEP B. Object ordering (Brusco_cov.f90, "observation sequencing").
  ## Put clusters in increasing order of size, take their members in turns,
  ## then reverse the sequence.
  ## ---------------------------------------------------------------------------
  members <- split(seq_len(N), factor(hg, levels = seq_len(G)))
  members <- members[order(lengths(members))]
  maxlen <- max(lengths(members))
  S <- unlist(lapply(seq_len(maxlen), function(k)
    unlist(lapply(members, function(m) if (k <= length(m)) m[k] else NULL))))
  S <- rev(S)

  ## ---------------------------------------------------------------------------
  ## STEP C. Pairwise matrices for the reordered units (COMPUTE_MATRIX),
  ## stacked in one array PW[i, j, q], q = 1..Qd:
  ##   q = 1                  : A_ij
  ##   q = 1 + k              : B_ij[k]       (k = 1..K)
  ##   q = 1 + K + (l-1)K + k : C_ij[k, l]    (column-major K x K)
  ## ---------------------------------------------------------------------------
  Ym <- matrix(panel$y, T, N)[, S, drop = FALSE]                  # T x N
  Xa <- array(panel$Z[, -1, S, drop = FALSE], dim = c(T, K, N))   # T x K x N
  pair_cross <- function(U, V) {
    ## sum_t (u_it - u_jt)(v_it - v_jt) for all pairs (i, j)
    su <- colSums(U * V)
    outer(su, rep(1, N)) + outer(rep(1, N), su) - crossprod(U, V) - crossprod(V, U)
  }
  Qd <- 1 + K + K * K
  PW <- array(0, dim = c(N, N, Qd))
  PW[, , 1] <- pair_cross(Ym, Ym)
  for (k in seq_len(K)) {
    PW[, , 1 + k] <- pair_cross(Ym, Xa[, k, ])
    for (l in seq_len(K)) PW[, , 1 + K + (l - 1) * K + k] <- pair_cross(Xa[, k, ], Xa[, l, ])
  }

  ## opt[r + 1] = optimal objective for the last r objects of the sequence
  ## (BUP2 in Brusco_cov.f90); with r <= G objects each can sit alone, so the
  ## optimum is zero.
  opt <- numeric(N + 1)
  nodes <- 0
  XB_final <- NULL

  ## ---------------------------------------------------------------------------
  ## STEP D. Repetitive branch and bound: E = G+1, ..., N objects.
  ## ---------------------------------------------------------------------------
  for (E in (G + 1):N) {
    sub <- (N - E + 1):N                  # positions of the last E objects

    ## Upper bound for this stage (ZXB): the best of the heuristic grouping
    ## restricted to these units and Algorithm 1 on the subset. The optimum
    ## of a subset never exceeds the objective of any feasible grouping of it,
    ## so these are valid incumbents. A small slack lets the exact search find
    ## a solution equal to the bound (Fortran: ZXB + 0.00001).
    sub_units <- S[sub]
    sub_panel <- .subset_panel(panel, sub_units)
    restricted <- hg[sub_units]
    restricted <- match(restricted, unique(restricted))     # relabel 1..G'
    ub <- if (max(restricted) == G) gfe_objective(sub_panel, restricted, G) else Inf
    if (stage_starts > 0) ub <- min(ub, gfe_alg1(sub_panel, G, n_starts = stage_starts)$objective)
    Z <- ub + 1e-6 * max(1, abs(ub))

    PWsub <- PW[sub, sub, , drop = FALSE]
    res <- if (engine == "cpp") bb_stage_cpp(as.vector(PWsub), E, G, K, Z, opt)
           else .bb_stage_R(PWsub, E, G, K, Z, opt)
    nodes <- nodes + res$nodes
    if (!isTRUE(res$found))
      stop("Branch and bound found no solution below the upper bound at stage E = ", E, ".")
    opt[E + 1] <- res$Z
    if (verbose) message(sprintf("NUMBER OF OBJECTS %3d  Z = %.6f  (nodes so far: %.0f)", E, res$Z, nodes))
    if (E == N) XB_final <- res$XB
  }

  ## Map back to the original unit ordering (XREOR in Brusco_cov.f90).
  groups <- integer(N)
  groups[S] <- XB_final
  upd <- gfe_update(panel, groups, G)
  fit <- .as_gfe(panel, G, c(upd[c("theta", "alpha", "objective")], list(groups = groups)),
                 opt[N + 1], algorithm = paste0("exact branch and bound (Brusco, 2006; with covariates; ",
                                                engine, " engine)"),
                 tuning = list(stage_starts = stage_starts, engine = engine))
  fit$stage_optima <- opt[-1]
  fit$nodes <- nodes
  fit$elapsed <- difftime(Sys.time(), t_start, units = "secs")
  fit$heuristic_objective <- heuristic$objective
  fit
}

#' One stage of the branch and bound, in R (reference implementation).
#'
#' The steps and their numbers follow Brusco_cov.f90 (GO TO 100, 200, ...).
#' Objects 1..E of the subproblem are assigned in order; clusters are opened
#' in order 1, 2, ..., G, which removes label symmetry.
.bb_stage_R <- function(PWsub, E, G, K, Z, opt) {
  Qd <- 1 + K + K * K
  rows <- lapply(seq_len(E), function(P) matrix(PWsub[P, , ], E, Qd))

  ## Concentrated partial objective a - b' C^{-1} b from the stacked cluster
  ## averages v = (a, b, vec(C)); NA when C is singular (the Fortran code
  ## skips the bound test in that case). Fortran tests ABS(DISC) > 1*10**(-6);
  ## with integer arithmetic 10**(-6) evaluates to 0 there, so the effective
  ## test is DISC != 0. We use a small relative tolerance instead.
  partial_obj <- function(v) {
    if (K == 0) return(v[1])
    if (K == 1) return(if (v[3] > 1e-12) v[1] - v[2]^2 / v[3] else NA_real_)
    if (K == 2) {
      c11 <- v[4]; c21 <- v[5]; c22 <- v[7]; dt <- c11 * c22 - c21 * c21
      if (!(dt > 1e-12 * c11 * c22)) return(NA_real_)
      return(v[1] - (c22 * v[2]^2 - 2 * c21 * v[2] * v[3] + c11 * v[3]^2) / dt)
    }
    Cm <- matrix(v[(K + 2):Qd], K, K); b <- v[2:(K + 1)]
    if (rcond(Cm) < 1e-12) return(NA_real_)
    v[1] - drop(crossprod(b, solve(Cm, b)))
  }

  ## ---- STEP 0: initialize ----------------------------------------------------
  P <- 0; Q <- G                          # Q = number of still-empty clusters
  n <- integer(G)
  x <- integer(E)
  XB <- NULL
  sums <- matrix(0, G, Qd)                # cluster sums (SUMA, SUMB1-2, SUMC1-3)
  nodes <- 0

  add_obj <- function(P, M, sign) {
    ## add (sign = +1) or remove (sign = -1) the pairs (P, j), j < P, j in M
    if (P <= 1) return(invisible())
    js <- which(x[seq_len(P - 1)] == M)
    if (!length(js)) return(invisible())
    sums[M, ] <<- sums[M, ] + sign * colSums(rows[[P]][js, , drop = FALSE])
  }

  step <- 1
  M <- 1L
  repeat {
    if (step == 1) {
      ## ---- STEP 1: increment search depth -----------------------------------
      P <- P + 1; M <- 1L
      n[M] <- n[M] + 1L; if (n[M] == 1L) Q <- Q - 1
      x[P] <- M
      add_obj(P, M, +1)
      step <- 2
    } else if (step == 2) {
      ## ---- STEP 2: feasibility (enough objects left to fill empty groups) --
      step <- if (E - P < Q) 5 else 3
    } else if (step == 3) {
      ## ---- STEP 3: sub-optimality (bound test) -------------------------------
      nodes <- nodes + 1
      used <- seq_len(G - Q)              # clusters are opened in order 1, 2, ...
      ZB <- partial_obj(colSums(sums[used, , drop = FALSE] / n[used]))
      if (!is.na(ZB) && Z <= ZB + opt[E - P + 1]) { step <- 5; next }
      ## ---- STEP 4: update incumbent at a complete assignment -----------------
      if (P != E) { step <- 1; next }
      if (!is.na(ZB)) { Z <- ZB; XB <- x }
      step <- 5
    } else if (step == 5) {
      ## ---- STEP 5: action after fathoming -------------------------------------
      ## Branching right is pointless if M is the last cluster, or if P is
      ## alone in M (moving it would only relabel a new empty cluster).
      step <- if (M == G || n[M] == 1L) 7 else 6
    } else if (step == 6) {
      ## ---- STEP 6: branch right (move object P to the next cluster) ---------
      add_obj(P, M, -1)
      n[M] <- n[M] - 1L
      M <- M + 1L
      n[M] <- n[M] + 1L; if (n[M] == 1L) Q <- Q - 1
      x[P] <- M
      add_obj(P, M, +1)
      step <- 2
    } else if (step == 7) {
      ## ---- STEP 7: retract depth ----------------------------------------------
      x[P] <- 0L
      add_obj(P, M, -1)
      n[M] <- n[M] - 1L
      P <- P - 1
      if (n[M] == 0L) Q <- Q + 1
      if (P > 0) { M <- x[P]; step <- 5 } else break
    }
  }
  list(Z = Z, found = !is.null(XB), XB = XB, nodes = nodes)
}

## Compile src/gfe_exact_bb.cpp once per session (needs Rcpp and a compiler;
## on Windows, Rtools).
.gfe_load_cpp <- function() {
  if (exists("bb_stage_cpp", mode = "function")) return(TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
  src <- file.path(getOption("gfe.root", "."), "src", "gfe_exact_bb.cpp")
  if (!file.exists(src)) return(FALSE)
  ok <- tryCatch({ Rcpp::sourceCpp(src, env = globalenv()); TRUE },
                 error = function(e) { message("Could not compile ", src, ": ", conditionMessage(e)); FALSE })
  ok
}

## Panel restricted to a subset of units (balanced panels).
.subset_panel <- function(panel, units) {
  rows <- as.vector(outer(seq_len(panel$T), (units - 1) * panel$T, "+"))
  gfe_panel(panel$y[rows], panel$X[rows, , drop = FALSE], N = length(units), T = panel$T,
            unit_labels = panel$unit_labels[units], time_labels = panel$time_labels)
}

#' Brute-force GFE solution by enumerating all groupings (tiny N only).
#' Used in tests/ to check gfe_exact_bb() and the heuristics.
gfe_bruteforce <- function(panel, G) {
  N <- panel$N
  if (G^(N - 1) > 2e6) stop("too many groupings to enumerate")
  best <- list(objective = Inf)
  ## Unit 1 is fixed in group 1 to remove part of the label symmetry.
  grid <- as.matrix(expand.grid(rep(list(seq_len(G)), N - 1)))
  for (r in seq_len(nrow(grid))) {
    g <- c(1L, as.integer(grid[r, ]))
    if (length(unique(g)) < G) next
    o <- gfe_objective(panel, g, G)
    if (o < best$objective) best <- list(objective = o, groups = g)
  }
  best
}
