## =============================================================================
## scripts/00_build_data.R
##
## R reimplementation of Bonhomme and Manresa (2015), "Grouped Patterns of
## Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184.
##
## Builds the estimation sample of the empirical application from the raw
## Acemoglu, Johnson, Robinson, and Yared (2008, AER) 5-year panel, following
## the "AUTOMATIC CODE" part of GFE_replication_code.do (Replic = 1), and
## checks that it reproduces the authors' final_data.mat, which we ship as
## data/bm_final_data.csv.
##
##   outcome     y_it : Freedom House democracy index (fhpolrigaug)
##   covariates  x_it : lagged democracy, lagged log GDP per capita (lrgdpch),
##                      lags of one 5-year period, taken within country
##   sample          : balanced panel (samplebalancefe == 1), 90 countries,
##                      7 periods (1970, 1975, ..., 2000)
##
## Input: data/ajry_5yearpanel_extract.csv (variables extracted from the
## authors' copy of 5yearpanel.dta, restricted to the 90 sample countries and
## 1965-2000, the years needed for the lags), or the full .dta via `dta_path`.
## Original data: http://economics.mit.edu/faculty/acemoglu/data/ajry2008
## =============================================================================

source("gfe.R")

dta_path <- NULL   # e.g. "C:/.../Bonhomme_Manresa_codes/Application/5yearpanel.dta"

raw <- if (!is.null(dta_path)) {
  as.data.frame(haven::zap_labels(haven::read_dta(dta_path)))
} else {
  read.csv("data/ajry_5yearpanel_extract.csv", stringsAsFactors = FALSE)
}

## tsset ID TIME; sort ID TIME
raw <- raw[order(raw$code_numeric, raw$year), ]

## gen lag_dem = fhpolrigaug[_n-1] if ID[_n] == ID[_n-1]   (and likewise for income)
same_unit <- c(FALSE, raw$code_numeric[-1] == raw$code_numeric[-nrow(raw)])
lag1 <- function(v) { out <- c(NA, v[-length(v)]); out[!same_unit] <- NA; out }
raw$lag_dem <- lag1(raw$fhpolrigaug)
raw$lag_income <- lag1(raw$lrgdpch)

## drop if year < 1960; keep if samplebalancefe == 1
est <- raw[raw$year >= 1960 & !is.na(raw$samplebalancefe) & raw$samplebalancefe == 1, ]
est <- est[order(est$code_numeric, est$year), ]
stopifnot(nrow(est) == 630, all(table(est$year) == 90))

## Compare with the authors' final_data.mat (shipped as a CSV).
fin <- read.csv("data/bm_final_data.csv", stringsAsFactors = FALSE)
stopifnot(identical(est$code, fin$code), identical(as.integer(est$year), as.integer(fin$year)))
cat("Max absolute difference with final_data.mat:\n")
print(c(dem        = max(abs(est$fhpolrigaug - fin$dem)),
        lag_dem    = max(abs(est$lag_dem - fin$lag_dem)),
        lag_income = max(abs(est$lag_income - fin$lag_income))))
cat("Differences are at the 1e-7 level (single-precision storage in final_data.mat).\n",
    "All later scripts use data/bm_final_data.csv, i.e. the authors' exact numbers.\n")
