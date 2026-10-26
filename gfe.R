## =============================================================================
## gfe.R  --  Load the R implementation of the grouped fixed-effects estimator
##
##   Bonhomme, Stephane and Elena Manresa (2015), "Grouped Patterns of
##   Heterogeneity in Panel Data," Econometrica 83(3), 1147-1184.
##   doi:10.3982/ECTA11319
##
## The method and the original MATLAB, Fortran, and Stata replication code are
## by Stephane Bonhomme and Elena Manresa. This folder is an R translation for
## study and reuse; please cite their paper when you use it (see CITATION.md).
##
## Usage (from the project root, e.g. in RStudio with the .Rproj open):
##   source("gfe.R")
##
## Dependencies: base R only. Optional: Rcpp (compiled exact algorithm; on
## Windows it needs Rtools), haven (reading the original .dta files). The
## parallel package used for the bootstrap and Monte Carlo ships with R.
##
## R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
## GPL, version 2 or (at your option) any later version; see LICENSE. The
## license covers the translation only, not Bonhomme and Manresa's method or
## original code.
## =============================================================================

local({
  root <- getOption("gfe.root", ".")
  files <- file.path(root, "R", c("gfe_core.R", "gfe_exact.R", "gfe_inference.R",
                                  "gfe_select.R", "gfe_extensions.R", "gfe_utils.R"))
  missing <- files[!file.exists(files)]
  if (length(missing))
    stop("Cannot find ", paste(missing, collapse = ", "),
         ". Run from the project root or set options(gfe.root = '<path>').")
  for (f in files) sys.source(f, envir = globalenv(), keep.source = TRUE)
  options(gfe.root = normalizePath(root), gfe.source_files = normalizePath(files))
})

## Analysis data used in the paper: the 5-year balanced panel of Acemoglu,
## Johnson, Robinson, and Yared (2008), N = 90 countries, T = 7 (1970-2000),
## exactly as stored in the authors' final_data.mat.
load_bm_panel <- function(path = file.path(getOption("gfe.root", "."), "data", "bm_final_data.csv")) {
  dat <- read.csv(path, stringsAsFactors = FALSE)
  dat <- dat[order(dat$unit, dat$year), ]
  N <- length(unique(dat$unit)); T <- length(unique(dat$year))
  p <- gfe_panel(dat$dem, as.matrix(dat[, c("lag_dem", "lag_income")]), N = N, T = T,
                 unit_labels = dat$code[dat$year == min(dat$year)],
                 time_labels = sort(unique(dat$year)))
  attr(p, "data") <- dat
  p
}

## Published groupings and bootstrap standard errors (transcribed from the
## authors' MATLAB files by tools/make_published_groupings.R).
load_bm_published <- function(root = getOption("gfe.root", ".")) {
  list(groupings = read.csv(file.path(root, "data", "bm_published_groupings.csv"), stringsAsFactors = FALSE),
       bootstrap_se = read.csv(file.path(root, "data", "bm_published_bootstrap_se.csv"), stringsAsFactors = FALSE))
}
