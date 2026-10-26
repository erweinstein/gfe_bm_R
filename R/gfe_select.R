## =============================================================================
## R/gfe_select.R  --  Choosing the number of groups
##
## R translation of the information criteria in Bonhomme and Manresa (2015),
## Supplementary Appendix Section S.3.2, as computed in GFE_estimates.m,
## GFE_deviations_estimates.m, and MonteCarlo/BIC_G3_MC.m (original code by
## Stephane Bonhomme and Elena Manresa).
##
##   BIC(G) = Q(G)/(NT) + sigma2_hat * (G T + N + K)/(NT) * ln(NT)          (S.20)
##
## where Q(G) is the minimized GFE objective with G groups and sigma2_hat is a
## consistent estimate of Var(v_it). Following footnote 12 of the Supplementary
## Appendix (and the code), sigma2_hat is computed once, from the fit with the
## largest number of groups considered, Gmax:
##
##   sigma2_hat = Q(Gmax) / (NT - Gmax T - N - K).
##
## BIC_G3_MC.m also computes an alternative with a steeper penalty,
##
##   BIC2(G) = Q(G)/(NT) + sigma2_hat * G (N + T + K - G)/(NT) * ln(NT),
##
## reported in the Supplementary Appendix as "an alternative choice with a
## steeper penalty". Both are returned.
##
## Consistency: G_hat is consistent for G0 if N and T grow at the same rate;
## if T/N -> 0 the BIC need not be consistent (Section S.3.2).
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

#' Information criteria for a list of GFE fits with different G.
#'
#' @param fits a list of "gfe" objects for G = 1, ..., Gmax (any order). The
#'   G = 1 fit is the model with common time effects only (two-way FE without
#'   unit effects), as the OLS benchmark in GFE_estimates.m.
#' @return data frame with G, objective, BIC, BIC2, and the selected G.
gfe_bic <- function(fits) {
  Gs <- vapply(fits, `[[`, numeric(1), "G")
  fits <- fits[order(Gs)]; Gs <- sort(Gs)
  N <- fits[[1]]$N; T <- fits[[1]]$T; K <- fits[[1]]$K
  NT <- sum(fits[[1]]$panel$d)
  Gmax <- max(Gs)
  Qmax <- fits[[which.max(Gs)]]$objective
  sigma2 <- Qmax / (NT - Gmax * T - N - K)
  Q <- vapply(fits, `[[`, numeric(1), "objective")
  out <- data.frame(
    G = Gs,
    objective = Q,
    BIC = Q / NT + sigma2 * log(NT) * (Gs * T + N + K) / NT,
    BIC2 = Q / NT + sigma2 * log(NT) * Gs * (N + T + K - Gs) / NT
  )
  attr(out, "sigma2") <- sigma2
  attr(out, "G_BIC") <- out$G[which.min(out$BIC)]
  attr(out, "G_BIC2") <- out$G[which.min(out$BIC2)]
  out
}
