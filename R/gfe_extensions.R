## =============================================================================
## R/gfe_extensions.R  --  Extensions of the baseline GFE model
##
## R translation of the extensions in Bonhomme and Manresa (2015),
## Supplementary Appendix Section S.4 and the corresponding replication files
## (original code by Stephane Bonhomme and Elena Manresa):
##
##   (1) Unit-specific heterogeneity: GFE in deviations from unit means,
##       y_it - ybar_i = (x_it - xbar_i)' theta + mu_{g_i t} + (v_it - vbar_i)
##       [GFE_deviations_estimates.m; GFE_replication_code.do with Replic = 5;
##        equations (S.21)-(S.22)], and the first-differenced IV estimator that
##       corrects the small-T bias with a lagged outcome
##       [GFE_replication_code_IV.do; equation (S.30); Table S.XV].
##
##   (2) Heterogeneous (group-specific) coefficients:
##       y_it = x_it' theta_{g_i} + alpha_{g_i t} + v_it
##       [Heterogeneous_coefficientsA.m (all slopes group-specific) and
##        Heterogeneous_coefficientsB.m (common lagged-outcome coefficient,
##        group-specific income coefficient); equation (S.32)].
##
##   (3) Two layers of groups: G "large" groups, each split into H_g subgroups
##       whose time profiles differ from the group profile only by a constant:
##       alpha_{(g,h) t} = a_{g t} + b_{g h}
##       [GFE_twolayers_computation.m and GFE_twolayers_estimates.m; Figure S.6].
## =============================================================================


## -----------------------------------------------------------------------------
## (1a) Deviations from unit means
## -----------------------------------------------------------------------------

#' Transform a panel into deviations from unit means (within transformation),
#' so that GFE removes time-invariant unit effects (model (5) in the paper).
gfe_within <- function(panel) {
  T <- panel$T; N <- panel$N
  dm <- function(v) {
    M <- matrix(v, T, N)
    mu <- colSums(M * panel$Dm) / colSums(panel$Dm)
    as.vector((M - rep(mu, each = T)) * panel$Dm)
  }
  X <- if (panel$K > 0) apply(panel$X, 2, dm) else NULL
  if (!is.null(X)) X <- matrix(X, ncol = panel$K, dimnames = list(NULL, colnames(panel$X)))
  gfe_panel(dm(panel$y), X, N, T, d = panel$d,
            unit_labels = panel$unit_labels, time_labels = panel$time_labels)
}


## -----------------------------------------------------------------------------
## (1b) First-differenced IV with estimated groups (Table S.XV)
## -----------------------------------------------------------------------------

#' First-differenced 2SLS with group x time effects, as in
#' GFE_replication_code_IV.do (Replic = 5):
#'
#'   Delta y_it = theta_1 Delta y_{i,t-1} + theta_2 Delta x_{i,t-1} + mu_{g_i t} + Delta v_it,
#'
#' instrumenting Delta y_{i,t-1} with y_{i,t-2}, and including (estimated)
#' group x period dummies as exogenous regressors. Standard errors are
#' clustered by unit (vce(cluster code_numeric)), without a small-sample
#' factor, as in Stata's ivregress; this reproduces Table S.XV.
#'
#' As in the .do file, lags are taken within the estimation sample (the 1970-
#' 2000 balanced panel), so the first two periods are lost.
#'
#' @param df long data frame with columns id, year, y (outcome level), x
#'   (covariate level, e.g. log income), and group (estimated GFE group).
gfe_fd_iv <- function(df, id = "id", year = "year", y = "y", x = "x", group = "group") {
  df <- df[order(df[[id]], df[[year]]), ]
  lagw <- function(v, k = 1) {
    out <- v; n <- length(v)
    same <- c(rep(FALSE, k), df[[id]][-(1:k)] == df[[id]][1:(n - k)])
    out[] <- NA; out[same] <- v[which(same) - k]
    out
  }
  L1y <- lagw(df[[y]], 1); L2y <- lagw(df[[y]], 2)
  L1x <- lagw(df[[x]], 1); L2x <- lagw(df[[x]], 2)
  dd <- data.frame(id = df[[id]], year = df[[year]], group = df[[group]],
                   delta_y = df[[y]] - L1y,
                   lag_delta_y = L1y - L2y,
                   lag_delta_x = L1x - L2x,
                   lag2_y = L2y)
  dd <- dd[stats::complete.cases(dd), ]
  cell <- factor(paste(dd$group, dd$year))
  Dgt <- stats::model.matrix(~ cell - 1)
  X <- cbind(lag_delta_y = dd$lag_delta_y, lag_delta_x = dd$lag_delta_x, Dgt)
  Z <- cbind(lag2_y = dd$lag2_y, lag_delta_x = dd$lag_delta_x, Dgt)
  fit <- ols_iv_cluster(dd$delta_y, X, Z, cluster = dd$id)
  k <- c("lag_delta_y", "lag_delta_x")
  list(coef = fit$coef[k], se = fit$se[k], vcov = fit$vcov[k, k],
       n = nrow(dd), n_clusters = fit$n_clusters)
}


## -----------------------------------------------------------------------------
## (2) Heterogeneous coefficients (Heterogeneous_coefficientsA/B.m)
## -----------------------------------------------------------------------------

#' GFE with group-specific slope coefficients, computed with "a counterpart to
#' Algorithm 1" (footnote 17 of the Supplementary Appendix), as in
#' Heterogeneous_coefficientsA.m / B.m:
#'
#'   starting values: each unit assigned to a random group;
#'   update:     OLS of y on [group x time dummies, common covariates,
#'               heterogeneous covariates x group dummies]   (pinv, as in MATLAB)
#'   assignment: each unit to the group with the smallest sum of squared
#'               residuals, using that group's slopes and time effects;
#'   repeat until the parameters no longer change; keep the best of n_starts.
#'
#' @param hetero indices of the covariates with group-specific coefficients
#'   (A.m: all of them; B.m: only income, i.e. hetero = 2).
gfe_hetero <- function(panel, G, hetero = seq_len(panel$K), n_starts = 1000,
                       max_iter = 1000, polish_top = 0, start_groups = NULL) {
  N <- panel$N; T <- panel$T; K <- panel$K
  common <- setdiff(seq_len(K), hetero)
  Xc <- panel$X[, common, drop = FALSE]
  Xh <- panel$X[, hetero, drop = FALSE]
  kh <- length(hetero); kc <- length(common)

  design <- function(groups) {
    Dg <- matrix(0, N * T, G * T)
    Dg[cbind(seq_len(N * T), (groups[panel$unit] - 1L) * T + panel$time)] <- 1
    ## heterogeneous slopes: Xh interacted with group dummies (kron(X_i, gi) in
    ## the MATLAB code), ordered covariate-major: (k1:g1..gG, k2:g1..gG, ...)
    Xhg <- matrix(0, N * T, kh * G)
    for (k in seq_len(kh)) Xhg[cbind(seq_len(N * T), (k - 1) * G + groups[panel$unit])] <- Xh[, k]
    cbind(Dg, Xc, Xhg)
  }
  unpack <- function(par) {
    alpha <- matrix(par[seq_len(G * T)], T, G)
    theta_c <- par[G * T + seq_len(kc)]
    beta <- matrix(par[G * T + kc + seq_len(kh * G)], G, kh)   # row g: slopes of group g
    list(alpha = alpha, theta_c = theta_c, beta = beta)
  }
  ssr_by_group <- function(pr) {
    vapply(seq_len(G), function(g) {
      r <- panel$y - pr$alpha[panel$time, g]
      if (kc) r <- r - drop(Xc %*% pr$theta_c)
      r <- r - drop(Xh %*% pr$beta[g, ])
      colSums(matrix(r^2, T, N))
    }, numeric(N))
  }

  ssr_of <- function(groups) {
    if (length(unique(groups)) < G) return(Inf)
    sum(stats::.lm.fit(design(groups), panel$y)$residuals^2)
  }
  run_alg1 <- function(groups) {
    for (iter in seq_len(max_iter)) {
      if (length(unique(groups)) < G) return(list(objective = Inf, groups = groups))
      Xtot <- design(groups)
      par <- drop(pinv(Xtot) %*% panel$y)
      new_groups <- max.col(-ssr_by_group(unpack(par)), ties.method = "first")
      if (identical(new_groups, groups)) break
      groups <- new_groups
    }
    list(objective = ssr_of(groups), groups = groups)
  }

  ## Random starts (or user-supplied starting groupings).
  starts <- if (is.null(start_groups)) n_starts else length(start_groups)
  sols <- vector("list", starts)
  objs <- numeric(starts)
  for (s in seq_len(starts)) {
    g0 <- if (is.null(start_groups)) sample.int(G, N, replace = TRUE) else start_groups[[s]]
    sols[[s]] <- run_alg1(g0)                          # random initial assignment
    objs[s] <- sols[[s]]$objective
  }

  ## Optional refinement (not in the original code): single-unit local
  ## search, as in Step 5 of Algorithm 2, applied to the best `polish_top`
  ## distinct solutions. Useful when the number of starts is far below the
  ## 1,000,000 used in the paper.
  if (polish_top > 0) {
    keys <- vapply(sols, function(s) paste(s$groups, collapse = ""), "")
    cand <- which(!duplicated(keys) & is.finite(objs))
    cand <- cand[order(objs[cand])][seq_len(min(polish_top, length(cand)))]
    for (s in cand) {
      g <- sols[[s]]$groups; obj <- sols[[s]]$objective
      repeat {
        improved <- FALSE
        for (i in seq_len(N)) {
          if (sum(g == g[i]) <= 1) next
          for (h in seq_len(G)[-g[i]]) {
            g2 <- g; g2[i] <- h
            o <- ssr_of(g2)
            if (o < obj - 1e-10 * obj) { g <- g2; obj <- o; improved <- TRUE }
          }
        }
        if (!improved) break
        a1 <- run_alg1(g)                              # re-run the iterations
        if (a1$objective < obj) { g <- a1$groups; obj <- a1$objective }
      }
      sols[[s]] <- list(objective = obj, groups = g)
      objs[s] <- obj
    }
  }
  b <- which.min(objs)
  Xtot <- design(sols[[b]]$groups)
  best <- list(objective = objs[b], groups = sols[[b]]$groups,
               par = drop(pinv(Xtot) %*% panel$y))

  ## Large-T clustered standard errors as in the MATLAB code
  ## (including its degrees-of-freedom factor NT/(NT - T - K)).
  Xtot <- design(best$groups)
  e <- drop(panel$y - Xtot %*% best$par)
  bread <- pinv(crossprod(Xtot))
  score <- rowsum(Xtot * e, panel$unit)
  V <- bread %*% crossprod(score) %*% bread * N * T / (N * T - T - K)
  pr <- unpack(best$par)
  se_all <- sqrt(pmax(diag(V), 0))
  se <- unpack(se_all)
  dimnames(pr$beta) <- dimnames(se$beta) <- list(paste0("g", seq_len(G)), colnames(panel$X)[hetero])
  structure(list(groups = best$groups, objective = best$objective,
                 alpha = pr$alpha, theta_common = setNames(pr$theta_c, colnames(panel$X)[common]),
                 beta = pr$beta, se_theta_common = se$theta_c, se_beta = se$beta,
                 vcov = V, objectives_by_start = objs, G = G, panel = panel),
            class = "gfe_hetero")
}


## -----------------------------------------------------------------------------
## (3) Two layers of groups (GFE_twolayers_computation.m / _estimates.m)
## -----------------------------------------------------------------------------

#' Linear restrictions R c = 0 on the stacked subgroup time effects
#' c = (c_{s,t}), s = 1..sum(H) (subgroup-major), imposing
#'   c_{(g,h),t} - c_{(g,1),t} - c_{(g,h),1} + c_{(g,1),1} = 0,  h >= 2, t >= 2,
#' i.e. alpha_{(g,h)t} = a_{gt} + b_{gh} (Mat_cons in the MATLAB code).
.twolayer_constraints <- function(H, T, K) {
  G1 <- length(H); S <- sum(H)
  start <- c(0, cumsum(H))[seq_len(G1)]            # subgroups before group g
  R <- matrix(0, sum(H - 1) * (T - 1), S * T + K)
  row <- 0
  for (g in seq_len(G1)) for (h in seq_len(H[g])[-1]) for (t in 2:T) {
    row <- row + 1
    s1 <- start[g] + 1; sh <- start[g] + h
    R[row, (sh - 1) * T + t] <- 1
    R[row, (s1 - 1) * T + t] <- -1
    R[row, (sh - 1) * T + 1] <- -1
    R[row, (s1 - 1) * T + 1] <- 1
  }
  R
}

#' Two-layer GFE (Figure S.6, lower panel; Table S.XIII rightmost column).
#'
#' Units are assigned to one of sum(H) subgroups; the subgroup time profiles
#' are constrained so that, within layer-1 group g, they differ only by
#' constants. Computation follows GFE_twolayers_computation.m: Algorithm 1
#' with the update step replaced by constrained least squares,
#'   par = P pinv(Xtot P) y,  P = I - pinv(R) R.
#'
#' @param H vector of subgroup counts, e.g. c(5, 2, 2).
#' @param groups optional subgroup labels (1..sum(H)) to evaluate directly.
gfe_twolayer <- function(panel, H, n_starts = 1000, groups = NULL, max_iter = 1000) {
  N <- panel$N; T <- panel$T; K <- panel$K; S <- sum(H)
  R <- .twolayer_constraints(H, T, K)
  P <- diag(S * T + K) - pinv(R) %*% R
  design <- function(groups) {
    Dg <- matrix(0, N * T, S * T)
    Dg[cbind(seq_len(N * T), (groups[panel$unit] - 1L) * T + panel$time)] <- 1
    cbind(Dg, panel$X)
  }
  fit_given <- function(groups) {
    Xtot <- design(groups)
    par <- drop(P %*% pinv(Xtot %*% P) %*% panel$y)
    list(par = par, c = matrix(par[seq_len(S * T)], T, S), theta = par[S * T + seq_len(K)],
         objective = sum((panel$y - Xtot %*% par)^2))
  }
  if (is.null(groups)) {
    best <- list(objective = Inf)
    for (s in seq_len(n_starts)) {
      st <- .random_start(panel, S)                    # S random "centers"
      theta <- st$theta; cc <- st$alpha; g <- rep(0L, N)
      for (iter in seq_len(max_iter)) {
        new_g <- max.col(-gfe_distances(panel, theta, cc), ties.method = "first")
        if (identical(new_g, g)) break
        g <- new_g
        if (length(unique(g)) < S) break               # an empty subgroup: stop this start
        f <- fit_given(g); theta <- f$theta; cc <- f$c
      }
      if (length(unique(g)) < S) next
      f <- fit_given(g)
      if (f$objective < best$objective) best <- c(f, list(groups = g))
    }
  } else {
    best <- c(fit_given(as.integer(groups)), list(groups = as.integer(groups)))
  }
  layer1 <- rep(seq_along(H), H)[best$groups]
  within <- unlist(lapply(H, seq_len))[best$groups]
  names(best$theta) <- colnames(panel$X)
  c(best, list(H = H, layer1 = layer1, layer2 = within, panel = panel))
}
