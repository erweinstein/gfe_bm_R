## =============================================================================
## scripts/06_explaining_the_waves.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184,
## Supplementary Appendix S.7.4 ("Explaining the estimated grouped patterns"):
## Table S.XVII (descriptive statistics by group) and Table S.XVIII
## (multinomial logits of group membership).
##
## Follows Explaining_the_waves.do (original code by Stephane Bonhomme and
## Elena Manresa), which runs Stata's mlogit on explaining_the_waves_G4.dta,
## shipped here as data/explaining_the_waves_G4.csv. The G = 4 groups are those
## of Figure 2: 1 high democracy, 2 low democracy (the base category),
## 3 early transition, 4 late transition.
##
## Treating estimated groups as data. Section S.5.3 of the Supplementary
## Appendix justifies using the GFE group estimates as a dependent variable,
## and computing standard errors as if they were known, when N and T are both
## large: misclassification then vanishes at an exponential rate.
##
## Variables (Acemoglu, Johnson, Robinson, and Yared, 2008):
##   GDP_1500_cons  log GDP per capita in 1500 (footnote 34 of the Appendix)
##   indyear        year of independence (displayed per 100 years)
##   consfirstaug   constraints on the executive at independence
##   lag_dem1970    democracy in 1965 (the lag of the 1970 observation)
##   lag_income1970 log GDP per capita in 1965
##   education1970  average years of schooling, 1970
##   rel_catho80, rel_protmg80  shares of Catholic and Protestant, 1980
## =============================================================================

source("gfe.R")
dir.create("output", showWarnings = FALSE)

w <- read.csv("data/explaining_the_waves_G4.csv", stringsAsFactors = FALSE)
w <- w[!is.na(w$assignment), ]                       # the 90 countries of the panel
w$indyear100 <- w$indyear / 100                      # "Independence Year/100"
group_names <- c("1 (High Dem.)", "2 (Low Dem.)", "3 (Early Trans.)", "4 (Late Trans.)")

## -----------------------------------------------------------------------------
## 1. Table S.XVII: means and standard deviations by group
## -----------------------------------------------------------------------------
vars <- c(`log GDP p.c. (1500)` = "GDP_1500_cons", `Independence Year` = "indyear",
          `Constraints` = "consfirstaug", `Democracy (1965)` = "lag_dem1970",
          `log GDP p.c. (1965)` = "lag_income1970", `Education (1970)` = "education1970",
          `Share Catholic (1980)` = "rel_catho80", `Share Protestant (1980)` = "rel_protmg80")
desc <- do.call(rbind, lapply(names(vars), function(lab) {
  v <- w[[vars[lab]]]
  m <- tapply(v, w$assignment, mean, na.rm = TRUE)
  s <- tapply(v, w$assignment, sd, na.rm = TRUE)
  data.frame(variable = lab, statistic = c("mean", "sd"), rbind(m, s), check.names = FALSE)
}))
desc <- rbind(desc, data.frame(variable = "Number of observations", statistic = "n",
                               rbind(table(w$assignment)), check.names = FALSE))
names(desc)[3:6] <- group_names
cat("Table S.XVII: descriptive statistics by group\n")
print(cbind(desc[, 1:2], round(desc[, 3:6], 3)), row.names = FALSE)
write.csv(desc, "output/tableS17_descriptives.csv", row.names = FALSE)
## Every entry agrees with Table S.XVII to the published precision, except two
## for log GDP p.c. in 1500: the Group 1 s.d. (0.290 here, 0.300 published) and
## the Group 3 mean (6.478 here, 6.49 published).

## -----------------------------------------------------------------------------
## 2. Multinomial logit by Newton-Raphson
## -----------------------------------------------------------------------------
## P(g_i = j | z_i) = exp(z_i' b_j) / sum_k exp(z_i' b_k),  b_base = 0.
## Written out (instead of calling nnet::multinom) because the log-likelihood
## is globally concave, so Newton's method converges in a few steps from zero,
## and the score and Hessian give both the usual and the robust ("clustered
## by country", one observation per country) variance directly.
mlogit <- function(y, Z, base, tol = 1e-12, max_iter = 100) {
  levs <- sort(unique(y)); alt <- setdiff(levs, base)
  n <- nrow(Z); k <- ncol(Z); J <- length(alt)
  Y <- sapply(alt, function(j) as.numeric(y == j))           # n x J
  B <- matrix(0, k, J)
  for (iter in seq_len(max_iter)) {
    E <- exp(Z %*% B); P <- E / (1 + rowSums(E))               # n x J
    score_i <- do.call(cbind, lapply(seq_len(J), function(j) Z * (Y[, j] - P[, j])))   # n x kJ
    H <- matrix(0, k * J, k * J)                               # minus the Hessian
    for (a in seq_len(J)) for (b in seq_len(J)) {
      wgt <- P[, a] * ((a == b) - P[, b])
      H[(a - 1) * k + 1:k, (b - 1) * k + 1:k] <- crossprod(Z * wgt, Z)
    }
    step <- solve(H, colSums(score_i))
    B <- B + matrix(step, k, J)
    if (max(abs(step)) < tol) break
  }
  E <- exp(Z %*% B); P <- E / (1 + rowSums(E))
  score_i <- do.call(cbind, lapply(seq_len(J), function(j) Z * (Y[, j] - P[, j])))
  V_oim <- solve(H)
  V_rob <- V_oim %*% crossprod(score_i) %*% V_oim * n / (n - 1)   # Stata's vce(cluster) factor
  nm <- as.vector(outer(colnames(Z), alt, function(v, j) paste0("g", j, ":", v)))
  dimnames(B) <- list(colnames(Z), paste0("g", alt))
  dimnames(V_oim) <- dimnames(V_rob) <- list(nm, nm)
  list(coef = B, se = matrix(sqrt(diag(V_oim)), k, J, dimnames = dimnames(B)),
       se_robust = matrix(sqrt(diag(V_rob)), k, J, dimnames = dimnames(B)),
       loglik = sum(log(cbind(1, E) / (1 + rowSums(E)))[cbind(seq_len(n), match(y, c(base, alt)))]),
       n = n, iterations = iter)
}

## -----------------------------------------------------------------------------
## 3. Table S.XVIII: columns (1)-(10), base category = Group 2
## -----------------------------------------------------------------------------
specs <- list(
  `(1)`  = c("GDP_1500_cons"),
  `(2)`  = c("GDP_1500_cons", "indyear100", "consfirstaug"),
  `(3)`  = c("GDP_1500_cons", "indyear100", "consfirstaug", "rel_catho80", "rel_protmg80"),
  `(4)`  = c("lag_dem1970", "lag_income1970"),
  `(5)`  = c("lag_dem1970", "education1970"),
  `(6)`  = c("lag_dem1970", "lag_income1970", "education1970"),
  `(7)`  = c("GDP_1500_cons", "indyear100", "consfirstaug", "education1970"),
  `(8)`  = c("GDP_1500_cons", "indyear100", "consfirstaug", "lag_dem1970", "education1970"),
  `(9)`  = c("GDP_1500_cons", "indyear100", "consfirstaug", "lag_dem1970", "lag_income1970", "education1970"),
  `(10)` = c("GDP_1500_cons", "indyear100", "consfirstaug", "lag_dem1970", "lag_income1970",
             "education1970", "rel_catho80", "rel_protmg80"))

fits <- lapply(specs, function(v) {
  d <- w[stats::complete.cases(w[, v]), ]               # listwise deletion, as in Stata
  Z <- cbind(`(Intercept)` = 1, as.matrix(d[, v, drop = FALSE]))
  mlogit(d$assignment, Z, base = 2)
})

rows <- c(names(vars)[names(vars) != "Independence Year"], "Independence Year/100")
row_var <- c(setdiff(vars, "indyear"), "indyear100")
row_order <- c(1, 8, 2:7)                               # the paper's row order
long <- do.call(rbind, lapply(names(fits), function(col) {
  f <- fits[[col]]
  do.call(rbind, lapply(c(1, 3, 4), function(g) {
    gname <- paste0("g", g)
    data.frame(column = col, group = g, variable = rows[row_order], var = row_var[row_order],
               coef = f$coef[match(row_var[row_order], rownames(f$coef)), gname],
               se = f$se[match(row_var[row_order], rownames(f$coef)), gname],
               se_robust = f$se_robust[match(row_var[row_order], rownames(f$coef)), gname],
               n = f$n)
  }))
}))
write.csv(long, "output/tableS18_mlogit.csv", row.names = FALSE)

## Display in the layout of Table S.XVIII: coefficient (standard error).
show <- function(g) {
  sub <- long[long$group == g, ]
  cell <- ifelse(is.na(sub$coef), "-", sprintf("%.3f (%.3f)", sub$coef, sub$se))
  out <- matrix(cell, ncol = length(specs), dimnames = list(rows[row_order], names(specs)))
  out
}
for (g in c(1, 3, 4)) {
  cat("\nTable S.XVIII: Group", g, "vs. Group 2 (low democracy)\n")
  print(noquote(show(g)))
}
cat("\nSample size in column (10):", fits[["(10)"]]$n, "(paper: N = 68)\n")
cat("Published, column (1), Group 1: 1.39 (0.971); Group 3: 0.959 (1.19); Group 4: -1.06 (1.14).\n")
cat("Explaining_the_waves.do calls mlogit without vce(cluster ...), i.e. Stata's default\n",
    "(observed information) standard errors, reported in column `se`; `se_robust` in the CSV\n",
    "file gives the clustered alternative mentioned in the notes to the table.\n")
