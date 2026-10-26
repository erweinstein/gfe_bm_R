## =============================================================================
## R/gfe_inference.R  --  Standard errors for the GFE estimator
##
## R translation of the inference procedures in Bonhomme and Manresa (2015):
##
##   (1) Large-N,T variance, clustered at the unit level (Arellano, 1987):
##       GFE_estimates.m ("Large - T variance") and Supplementary Appendix
##       equations (S.2)-(S.3). Group membership is treated as known, which
##       Corollary 1 of the paper justifies as N and T grow.
##
##   (2) Large-N, fixed-T variance: fixedT_function.m, "an extension of
##       Pollard's (1982) fixed-T variance formula" to models with covariates
##       (Supplementary Appendix S.2.2, equations (S.7)-(S.16)). It adds the
##       contribution of units near the frontier between two groups, which the
##       large-T formula ignores.
##
##   (3) Bootstrap: resampling unit-specific blocks (y_i, x_i) with replacement
##       and re-estimating the GFE model in each bootstrap sample (Supplementary
##       Appendix S.2.2.2; Bootstrap_version.exe / MC_bootstrap.exe).
##
## The original code is by Stephane Bonhomme and Elena Manresa; this is a
## translation into R.
## =============================================================================


## -----------------------------------------------------------------------------
## Regressor matrix of the "known groups" regression: [group x time dummies, X]
## (Xtot = [kron(gi, eye(T)) X] in GFE_estimates.m). Columns: alpha_{g,t} in
## group-major order (g1t1, ..., g1tT, g2t1, ...), then theta.
## -----------------------------------------------------------------------------
.gfe_design <- function(fit) {
  p <- fit$panel; T <- p$T; G <- fit$G
  D <- matrix(0, p$N * T, G * T)
  col <- (fit$groups[p$unit] - 1L) * T + p$time
  D[cbind(seq_len(p$N * T), col)] <- 1
  D <- D * p$d
  colnames(D) <- paste0("alpha_g", rep(seq_len(G), each = T), "_t", rep(p$time_labels, G))
  Xtot <- cbind(D, p$X)
  Xtot
}

#' Large-T variance clustered at the unit level (GFE_estimates.m).
#'
#'   V = (Xtot'Xtot)^{-1} [ sum_i Xtot_i' e_i e_i' Xtot_i ] (Xtot'Xtot)^{-1}
#'         * NT / (NT - GT - K),
#'
#' with e the GFE residuals. The degrees-of-freedom factor is the one used in
#' the original code. Returns the full (GT + K) x (GT + K) matrix; the theta
#' block is the large-N,T variance of theta_hat in (S.3).
gfe_vcov_largeT <- function(fit) {
  p <- fit$panel
  Xtot <- .gfe_design(fit)
  e <- p$d * (p$y - drop(Xtot %*% c(as.vector(fit$alpha), fit$theta)))
  e[is.na(e)] <- 0
  keep <- colSums(Xtot != 0) > 0               # drop empty group-period cells
  Xk <- Xtot[, keep, drop = FALSE]
  bread <- solve(crossprod(Xk))
  ## Meat: sum over units of (X_i' e_i)(X_i' e_i)'.
  score <- rowsum(Xk * e, p$unit)              # N x ncol: X_i' e_i for each unit
  meat <- crossprod(score)
  nobs <- sum(p$d)
  V <- bread %*% meat %*% bread * nobs / (nobs - fit$G * p$T - p$K)
  out <- matrix(NA_real_, ncol(Xtot), ncol(Xtot), dimnames = list(colnames(Xtot), colnames(Xtot)))
  out[keep, keep] <- V
  out
}

#' Convenience: standard errors of theta from a variance matrix.
gfe_se_theta <- function(fit, V) {
  k <- seq_len(fit$K) + fit$G * fit$T
  setNames(sqrt(diag(V)[k]), names(fit$theta))
}

#' Delta-method variance of a smooth function h(theta) (e.g. the cumulative
#' income effect theta_2 / (1 - theta_1) in the application).
#' @param grad gradient of h at theta_hat (length K).
gfe_delta <- function(fit, V, grad) {
  k <- seq_len(fit$K) + fit$G * fit$T
  drop(crossprod(grad, V[k, k] %*% grad))
}


## -----------------------------------------------------------------------------
## Fixed-T variance (fixedT_function.m)
## -----------------------------------------------------------------------------

#' Large-N, fixed-T variance of (alpha, theta): Gamma^{-1} V Gamma^{-1}.
#'
#' Faithful translation of fixedT_function.m. With a_g the T-vector of group
#' g's time effects and r_i = y_i - x_i theta (a T-vector):
#'
#'   * frontier weight of unit i between groups g and h (g != h):
#'       w_igh = phi( u_igh / eps ) / eps
#'               * 1{ g and h are both at least as close to r_i as every other group l },
#'       u_igh = (a_h - a_g)' (r_i - (a_g + a_h)/2) / ||a_h - a_g||,
#'     i.e. a Gaussian kernel in the signed distance from r_i to the
#'     hyperplane equidistant from a_g and a_h  (I_igh(eps_N) in (S.13)-(S.16)).
#'   * Gamma_xx = X'X/N - (1/2N) sum_{g != h} sum_i w_igh x_i' (a_h-a_g)(a_h-a_g)' x_i / ||a_h-a_g||
#'   * Gamma_xg = (1/N) sum_{i in g} x_i' + (1/N) sum_{h != g} sum_i w_igh x_i'(a_g-a_h)(r_i-a_g)' / ||a_h-a_g||
#'   * Gamma_gg = (n_g/N) I_T - (1/N) sum_{h != g} sum_i w_igh (r_i-a_g)(r_i-a_g)' / ||a_h-a_g||
#'   * Gamma_gh = (1/N) sum_i w_igh (r_i-a_g)(r_i-a_h)' / ||a_h-a_g||       (g != h)
#'   * V = sum_i W_i' e_i e_i' W_i with W_i = [g-dummies (x) I_T, x_i] (clustered).
#'
#' As in the original, Gamma is multiplied by N and V is not divided by N, so
#' Gamma^{-1} V Gamma^{-1} is directly the variance of the estimator.
#'
#' Notes on the original code, kept as is:
#'   - the "1/2" in Gamma_xx follows equation (S.13);
#'   - the blocks Gamma_gh are filled for g < h and then symmetrized;
#'   - alpha is recomputed from the given grouping and theta (delta_aux).
#'
#' @param eps bandwidth epsilon_N (kernel = standard normal density).
gfe_vcov_fixedT <- function(fit, eps) {
  p <- fit$panel
  if (!p$balanced) stop("gfe_vcov_fixedT() is written for balanced panels, as in the original.")
  N <- p$N; T <- p$T; K <- p$K; G <- fit$G
  theta <- fit$theta
  groups <- fit$groups
  r <- p$y - if (K > 0) drop(p$X %*% theta) else 0
  R <- matrix(r, T, N)                                   # column i: r_i
  a <- vapply(seq_len(G), function(g) rowMeans(R[, groups == g, drop = FALSE]), numeric(T))
  a <- matrix(a, T, G)
  Xi <- function(i) p$X[(i - 1) * T + seq_len(T), , drop = FALSE]   # T x K

  ## Squared distances of every unit to every group profile.
  dist <- vapply(seq_len(G), function(g) colSums((R - a[, g])^2), numeric(N))
  dist <- matrix(dist, N, G)

  ## Frontier weights w[i, g, h].
  w <- array(0, dim = c(N, G, G))
  for (g in seq_len(G)) for (h in seq_len(G)) if (h != g) {
    ind <- rep(1, N)
    for (l in seq_len(G)) if (l != g && l != h)
      ind <- ind * (dist[, g] <= dist[, l]) * (dist[, h] <= dist[, l])
    dvec <- a[, h] - a[, g]
    nrm <- sqrt(sum(dvec^2))
    u <- colSums((dvec / nrm) * (R - (a[, g] + a[, h]) / 2))
    w[, g, h] <- ind * stats::dnorm(u / eps) / eps
  }

  ## Gamma blocks.
  G_xx <- if (K > 0) crossprod(p$X) / N else matrix(0, 0, 0)
  G_xg <- array(0, dim = c(K, T, G))
  G_gg <- array(0, dim = c(T, T, G))
  G_gh <- array(0, dim = c(T, T, G, G))
  for (g in seq_len(G)) {
    G_gg[, , g] <- sum(groups == g) * diag(T) / N
    if (K > 0) for (i in which(groups == g)) G_xg[, , g] <- G_xg[, , g] + t(Xi(i)) / N
  }
  for (g in seq_len(G)) for (h in seq_len(G)) if (h != g) {
    ag <- a[, g]; ah <- a[, h]
    nrm <- sqrt(sum((ah - ag)^2))
    act <- which(w[, g, h] != 0)
    for (i in act) {
      wi <- w[i, g, h]
      ri <- R[, i]
      if (K > 0) {
        x <- Xi(i)
        G_xx <- G_xx + wi * t(x) %*% (-0.5 * tcrossprod(ah - ag)) %*% x / nrm / N
        G_xg[, , g] <- G_xg[, , g] + wi * t(x) %*% tcrossprod(ag - ah, ri - ag) / nrm / N
      }
      G_gg[, , g] <- G_gg[, , g] - wi * tcrossprod(ri - ag) / nrm / N
      G_gh[, , g, h] <- G_gh[, , g, h] + wi * tcrossprod(ri - ag, ri - ah) / nrm / N
    }
  }

  ## Assemble BigGamma ((GT + K) x (GT + K)): alpha blocks first, then theta.
  dimG <- G * T + K
  BG <- matrix(0, dimG, dimG)
  idx <- function(g) (g - 1) * T + seq_len(T)
  kk <- G * T + seq_len(K)
  if (K > 0) {
    BG[kk, kk] <- G_xx
    for (g in seq_len(G)) BG[idx(g), kk] <- t(G_xg[, , g])
  }
  for (g in seq_len(G)) for (h in g:G) {
    BG[idx(g), idx(h)] <- if (g == h) G_gg[, , g] else G_gh[, , g, h]
  }
  BG[lower.tri(BG)] <- t(BG)[lower.tri(BG)]              # symmetrize from the upper triangle
  BG <- N * BG

  ## Clustered "meat" V = Xtot' Omega Xtot (Omega block-diagonal by unit).
  Xtot <- .gfe_design(fit)
  e <- p$y - drop(Xtot %*% c(as.vector(a), theta))
  score <- rowsum(Xtot * e, p$unit)
  Vm <- crossprod(score)

  BGinv <- solve(BG)
  out <- BGinv %*% Vm %*% BGinv
  dimnames(out) <- list(colnames(Xtot), colnames(Xtot))
  out
}

#' Bandwidth rule for the fixed-T variance (GFE_estimates.m; footnote 10 of
#' the Supplementary Appendix): a Silverman-type rule of thumb,
#'
#'   eps_N = 1.06 * min( sigma_Y, IQR(y - x theta) / 1.34 ) * N^{-1/5},
#'
#' where sigma_Y is the smallest, over pairs of groups (g, h), of the standard
#' deviation of the residual paths projected on the unit normal vector
#' (a_g - a_h) / ||a_g - a_h||.
#'
#' GFE_estimates.m, as distributed, then picks the value closest to eps_N on
#' the grid linspace(0.05, 0.15, 10) (`snap_to_grid = TRUE`). The fixed-T
#' standard errors published in Table S.XI are reproduced by the unsnapped
#' rule (the default here): with it, they agree to within 0.001 for all G,
#' whereas snapping moves several of them by more (the formula is sensitive to
#' eps_N; see scripts/01_validate_against_published.R).
gfe_bandwidth <- function(fit, snap_to_grid = FALSE, grid = seq(0.05, 0.15, length.out = 10)) {
  p <- fit$panel; T <- p$T; N <- p$N; G <- fit$G
  r <- p$y - if (p$K > 0) drop(p$X %*% fit$theta) else 0
  R <- matrix(r, T, N)
  sig <- c()
  for (g in seq_len(G - 1)) for (h in (g + 1):G) {
    nv <- fit$alpha[, g] - fit$alpha[, h]
    nv <- nv / sqrt(sum(nv^2))
    sig <- c(sig, stats::sd(colSums(nv * R)))
  }
  sig_Y <- if (length(sig)) min(sig) else stats::sd(r)
  iqr <- diff(stats::quantile(r, c(0.25, 0.75), type = 5))   # MATLAB quantile convention
  bw <- 1.06 * min(sig_Y, iqr / 1.34) * N^(-1 / 5)
  if (snap_to_grid) grid[which.min(abs(bw - grid))] else unname(bw)
}


## -----------------------------------------------------------------------------
## Bootstrap
## -----------------------------------------------------------------------------

#' Unit-level (block) bootstrap of the GFE estimator.
#'
#' Resamples units (y_i, x_i) with replacement, re-estimates the GFE model in
#' each bootstrap sample, and returns the draws of theta (and of any function
#' of theta supplied in `fun`). Standard errors are the standard deviations
#' of the draws (as std_cluster_b in GFE_estimates.m); percentile intervals use
#' the 2.5% and 97.5% quantiles (as in Processing_bootstrap_G10_dgp1.m).
#'
#' Bootstrap consistency for minimum sum-of-squares partitioning is due to
#' Arcones and Gine (1992), building on Pollard (1982); see Section S.2.2.2.
#'
#' The paper's bootstrap used Algorithm 2 with (5; 10; 5) in the application,
#' and Algorithm 1 with 1,000 starts in the Monte Carlo exercises. Both are
#' expensive in R; reduce B or the tuning for a quick look.
#'
#' @param B number of bootstrap replications (100 in the paper).
#' @param algorithm,... passed to gfe() for each bootstrap sample.
#' @param fun optional function of theta returning extra quantities.
#' @param cores number of parallel workers (uses parallel::parLapply; works on
#'   Windows too). 1 = sequential.
gfe_bootstrap <- function(fit, B = 100, algorithm = c("vns", "alg1"), fun = NULL,
                          seed = NULL, cores = 1, ...) {
  algorithm <- match.arg(algorithm)
  p <- fit$panel; T <- p$T; N <- p$N; G <- fit$G
  if (!is.null(seed)) set.seed(seed)
  draws_units <- replicate(B, sample.int(N, N, replace = TRUE), simplify = FALSE)
  one <- function(units) {
    rows <- as.vector(outer(seq_len(T), (units - 1) * T, "+"))
    pb <- gfe_panel(p$y[rows], p$X[rows, , drop = FALSE], N = N, T = T, d = p$d[rows])
    fb <- gfe(pb, G, algorithm = algorithm, ...)
    c(fb$theta, if (!is.null(fun)) fun(fb$theta))
  }
  if (cores > 1) {
    ## Each worker sources the GFE code; `one` travels with its environment
    ## (the data and tuning parameters). Seeds: one stream per worker.
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterCall(cl, function(files) { for (f in files) source(f); NULL }, .gfe_source_files())
    parallel::clusterSetRNGStream(cl, if (is.null(seed)) NULL else seed)
    res <- parallel::parLapply(cl, draws_units, one)
  } else {
    res <- lapply(draws_units, one)
  }
  draws <- do.call(rbind, res)
  list(draws = draws,
       se = apply(draws, 2, stats::sd),
       percentile = apply(draws, 2, stats::quantile, probs = c(0.025, 0.975)),
       B = B, algorithm = algorithm)
}

## Files to source on parallel workers (set by gfe.R when it loads the code).
.gfe_source_files <- function() getOption("gfe.source_files", character(0))
