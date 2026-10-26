# Grouped fixed-effects (GFE) estimation in R

### An R translation of the replication code for Bonhomme and Manresa (2015)

> Bonhomme, Stéphane and Elena Manresa (2015), "Grouped Patterns of Heterogeneity
> in Panel Data," *Econometrica* 83(3), 1147–1184. doi:10.3982/ECTA11319

**R translation by**  
Eliot R. Weinstein  
eliotw@uchicago.edu  
[erweinstein.com](https://erweinstein.com)

**Credit.** The grouped fixed-effects estimator, its asymptotic theory, the
computational algorithms (Algorithm 1, variable neighborhood search, the exact
branch-and-bound extension), the fixed-T variance, the Monte Carlo design, the
empirical application, and the original MATLAB, Fortran, and Stata code are the
work of **Stéphane Bonhomme and Elena Manresa**. This folder translates their
replication package (`Bonhomme_Manresa_codes`, distributed with the paper by the
Econometric Society) into R, for study and reuse. Every source file names the
original file it follows. Please cite their paper when you use this code
(see `CITATION.md`). The R translation is by Eliot R. Weinstein; it was not
reviewed by the original authors, and errors introduced in translation are his.

The code has two goals:

1. **Preserve the method.** The estimator, the algorithms, and the inference
   procedures follow the original code step by step, and the scripts check the
   results against the numbers published in the paper and its Supplementary
   Appendix (see "Validation" below).
2. **Serve as a model of the research workflow.** The scripts follow the
   structure of the paper: the estimator as an optimization problem, the
   algorithms and a check of their accuracy, the main estimates, inference
   with a calibrated Monte Carlo, robustness organized around the main
   finding, and an attempt to explain it. `GFE_walkthrough.qmd` reads the
   paper alongside the code.

---

## Quick start (RStudio on Windows, macOS, or Linux)

1. Open `gfe_bm_R.Rproj` in RStudio (this sets the working directory).
2. In the console:

```r
source("gfe.R")                       # loads the functions in R/
p   <- load_bm_panel()                # 90 countries x 7 periods (1970-2000)
fit <- gfe(p, G = 4, algorithm = "vns", n_starts = 10, neighmax = 10, itermax = 10)
fit                                   # objective 14.31867, theta = (0.302, 0.082)
```

3. Run the scripts in `scripts/` in order, or render `GFE_walkthrough.qmd`
   (Render button in RStudio; needs Quarto, which ships with recent RStudio).

**Dependencies.** Base R only. Optional: `Rcpp` for the compiled exact
algorithm (on Windows, also install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R
version); `haven` to read the original `.dta` files; `knitr`/`rmarkdown` for
the walkthrough. Without Rcpp, `gfe_exact_bb()` falls back to a pure-R engine
that is practical for about N = 20.

**Parallel computing.** `gfe_bootstrap(..., cores = 4)` and
`scripts/04_montecarlo.R` (`cores <- 4`) use `parallel::parLapply`, which works
on Windows.

---

## What is here

```
LICENSE                   GNU GPL (>= 2) for the R translation, with scope note
CITATION.md               how to cite the paper, the data, and this translation
EXAMPLE.md                a first exercise: change G and the starting values
gfe.R                     loads the code; load_bm_panel(), load_bm_published()
R/gfe_core.R              the estimator: update step, assignment step,
                          Algorithm 1, local search, VNS (Algorithm 2)
R/gfe_exact.R             exact solution: Brusco's repetitive branch and bound
                          with covariates; brute force for tiny problems
src/gfe_exact_bb.cpp      the branch-and-bound stage loop in C++ (via Rcpp)
R/gfe_inference.R         large-T clustered variance, fixed-T (Pollard)
                          variance and bandwidth rule, unit-level bootstrap
R/gfe_select.R            BIC for the number of groups
R/gfe_extensions.R        deviations from unit means, first-differenced IV,
                          heterogeneous coefficients, two layers of groups
R/gfe_utils.R             Hungarian matching of group labels, group means,
                          plotting in the paper's style, pinv, cluster-robust IV
scripts/00_build_data.R   builds the panel from the AJRY (2008) data
scripts/01_validate_against_published.R   checks against Tables S.I, S.III,
                          S.XI, S.XIV, S.XV
scripts/02_application_income_democracy.R Table S.XI, Figures 1, 2, S.4
scripts/03_exact_algorithm.R              Table S.I (exact vs heuristics)
scripts/04_montecarlo.R                   Tables S.III-S.VIII, Figure S.1
scripts/05_extensions.R                   Tables S.XIV-S.XVI, Figures S.5-S.7
scripts/06_explaining_the_waves.R         Tables S.XVII-S.XVIII
tests/test_gfe.R          unit tests against brute force, lm(), direct formulas
GFE_walkthrough.qmd       a reading of the paper with the code (Quarto)
GFE_walkthrough.html      the same, rendered (opens in any browser)
tools/make_published_groupings.R  transcribes the authors' groupings
data/                     see below
output/                   tables (CSV) and figures (PDF) written by the scripts
```

### Data

| file | content | source |
|---|---|---|
| `bm_final_data.csv` | estimation sample: `dem`, `lag_dem`, `lag_income`, N = 90, T = 7 | the authors' `final_data.mat` |
| `ajry_5yearpanel_extract.csv` | raw variables used to build the sample (90 sample countries, 1965–2000) | the authors' copy of AJRY's `5yearpanel.dta` |
| `explaining_the_waves_G4.csv` | covariates for Tables S.XVII–S.XVIII | the authors' `explaining_the_waves_G4.dta` |
| `bm_published_groupings.csv` | optimal groupings G = 2..15, continents, deviations G = 2..5, two-layer (5,2,2) | hard-coded in `GFE_estimates.m`, `GFE_deviations_estimates.m`, `GFE_twolayers_estimates.m` |
| `bm_published_bootstrap_se.csv` | bootstrap standard errors, G = 2..15 | `std_cluster_b` in `GFE_estimates.m` |
| `rcode_hetero_groupings.csv` | groupings for the heterogeneous-coefficient models | found by this R code (not published by the authors) |

The underlying data are from Acemoglu, Johnson, Robinson, and Yared (2008),
"Income and Democracy," *AER* 98(3), 808–842,
[jstor.org/stable/29730096](https://www.jstor.org/stable/29730096). They remain
subject to their original terms (see "License" below).

---

## Map from the original files to the R code

| original (Bonhomme and Manresa) | R |
|---|---|
| `Application/Fortran codes/*.f90` (Algorithms 1 and 2) | `gfe_alg1()`, `gfe_vns()`, `gfe_local_search()` in `R/gfe_core.R` |
| `Application/Data_paths.m` (Algorithm 1 in MATLAB) | `gfe_alg1()` |
| `Application/Exact algorithm/Brusco_cov.f90`, `heuristics_sub.f90` | `gfe_exact_bb()`, `src/gfe_exact_bb.cpp` |
| `Application/GFE_estimates.m` | `gfe_update()`, `gfe_vcov_largeT()`, `gfe_bic()`; `scripts/02` |
| `Application/fixedT_function.m` | `gfe_vcov_fixedT()`, `gfe_bandwidth()` |
| `Application/GFE_code.do`, `GFE_replication_code.do` | `gfe_panel_from_df()`; `scripts/00`; map in `scripts/02` |
| `Application/GFE_deviations_estimates.m` | `gfe_within()`; `scripts/05` |
| `Application/GFE_replication_code_IV.do` | `gfe_fd_iv()`, `ols_iv_cluster()`; `scripts/05` |
| `Application/Heterogeneous_coefficientsA.m`, `B.m` | `gfe_hetero()`; `scripts/05` |
| `Application/GFE_twolayers_computation.m`, `_estimates.m` | `gfe_twolayer()`; `scripts/05` |
| `Application/Explaining_the_waves.do` | `scripts/06` |
| `MonteCarlo/DGP_final.m`, `MC1.exe`, `MC_bootstrap.exe` | `scripts/04` |
| `MonteCarlo/Processing_bootstrap_G10_dgp1.m`, `BIC_G3_MC.m` | `scripts/04` |

---

## Validation

All numbers below are produced by the scripts in this folder.

**Main estimates (Table S.XI).** Given the authors' groupings, the objective,
BIC, coefficients, cumulative income effect, and large-T standard errors match
the published table for G = 1..15, and the BIC selects G = 10 as in the paper.
Of the 168 comparisons with Table S.XI (including the fixed-T standard errors
and the objectives found by search), 162 agree to the third decimal and the
other 6 differ by at most 0.001 (see `output/validation_report.csv`).

**Searching from scratch.** VNS with (10; 10; 10) recovers the published
optimal groupings for G = 2, 3, 4, 5 with no country classified differently;
Algorithm 1 with 1,000 starts reaches the same optima for G = 2 and 3.

**Exact algorithm (Table S.I).** For the first 30 countries and G = 2, the
branch and bound (160 million nodes, under a minute in C++) gives objective
6.1591339 and theta = (0.6070, 0.0689), as in the authors' `RESULTS2.DAT` and
`outputtheta.txt`; Algorithm 1 and VNS reach the same solution. The R and C++
engines and brute-force enumeration agree on small problems, node for node.

**Monte Carlo (Tables S.III–S.IV).** The "true" values of the DGP match
Table S.III for G0 = 3, 5, 10. With the script's defaults (G0 = 3, i.i.d.
errors, 50 replications, about 5 minutes), the results line up with the
paper's 1,000 replications. Average misclassification is 9.2% (paper: 9.5%).
Median large-T standard errors are 0.035, 0.0094, 0.013 for theta1, theta2,
and the cumulative effect (paper: 0.035, 0.0093, 0.013), and median fixed-T
standard errors are 0.052, 0.014, 0.024 (paper: 0.051, 0.013, 0.022). Coverage
of the 95% large-T intervals is 0.78–0.88; for fixed-T intervals it is
0.92–0.98. With 50 replications the means and coverage rates are noisy (see
`output/montecarlo_G3_iid.csv`).

**Extensions.** Table S.XIV (deviations from country means): 20 of 20 entries.
Table S.XV (first-differenced IV): 20 of 20 coefficients and standard errors.
Table S.XVI version B (heterogeneous income coefficient): the R search
(20,000 random starts plus local search; 11 minutes) reproduces the published estimates and
standard errors, 0.288 (0.054) for lagged democracy and 0.103 (0.019),
0.047 (0.014), 0.087 (0.018), 0.082 (0.013) for income. Version A (both
coefficients group-specific) is reproduced too: 0.644, 0.319, 0.016, 0.248 for
lagged democracy and 0.070, 0.041, 0.122, 0.090 for income, with the published
standard errors. Table S.XIV's G = 1 (two-way fixed effects) row has slightly
larger published standard errors than the MATLAB formula gives (0.058 vs 0.053
for theta1), plausibly from Stata. The paper's statements
that the deviation and baseline classifications agree for 80 of 90 countries,
and the two-layer and deviation classifications "except in five cases", are
also reproduced.

**Explaining the waves (Tables S.XVII–S.XVIII).** All coefficients and
standard errors of the ten multinomial logits match, with N = 68 in column (10).
The descriptive statistics match except two entries for log GDP per capita in
1500 (Group 1 s.d. 0.290 vs 0.300; Group 3 mean 6.478 vs 6.49).

---

## Differences from the original code, and things we noticed

These are documented in the code where they occur.

- **Fixed-T bandwidth.** `GFE_estimates.m`, as distributed, snaps the
  rule-of-thumb bandwidth to the grid `linspace(0.05, 0.15, 10)`. The fixed-T
  standard errors published in Table S.XI are reproduced by the unsnapped rule,
  which is our default (`gfe_bandwidth(fit, snap_to_grid = TRUE)` gives the
  snapped version). The fixed-T variance is sensitive to the bandwidth.
- **Bootstrap standard errors for G = 7.** In `GFE_estimates.m` the G = 7
  entry of the plotted bootstrap bands reads `std_cluster2_b` instead of
  `std_cluster_b`, so Figure 1 reuses the G = 8 values at G = 7. The
  `note` column of `data/bm_published_bootstrap_se.csv` flags this.
- **G = 1, large-T standard error of theta1.** The formula in the MATLAB code
  gives 0.0483; the table reports 0.049, which is what Stata's `regress,
  vce(cluster)` small-sample factors would give.
- **First-differenced IV.** Stata's `ivregress` applies no small-sample factor
  to clustered standard errors; `ols_iv_cluster(small = "none")` matches
  Table S.XV. The standard errors of the cumulative effect in Table S.XV are not
  reproduced: the code that computed them is not in the replication files.
- **Exact algorithm tolerance.** `Brusco_cov.f90` tests `ABS(DISC) >
  1*10**(-6)`; with integer arithmetic, `10**(-6)` evaluates to 0 in Fortran.
  We use a small relative tolerance. The authors' input files store the data
  with 4 decimals; the published objective 6.1591 (N = 30) is for the rounded
  data (unrounded: 6.1596; the optimal grouping is the same).
- **Monte Carlo errors, alternative DGP.** The text of Section S.6.2 describes
  resampling unit-specific residual vectors with replacement; `DGP_final.m`
  (DGP = 2) multiplies each unit's residual vector by a standard normal draw
  (a wild bootstrap). `scripts/04_montecarlo.R` offers both (`dgp = "wild"`,
  `dgp = "resample"`).
- **BIC in the Monte Carlo.** `BIC_G3_MC.m` averages the variance estimate
  sigma2_hat across simulated samples and uses the average in every sample; the
  script does the same.
- **Misclassification.** The paper minimizes over all G! relabelings (500,000
  random ones for G = 10). `match_groups()` uses the Hungarian algorithm, which
  finds the exact minimum for any G.
- **Heterogeneous coefficients.** The paper uses 1,000,000 random starts.
  `gfe_hetero(..., polish_top = k)` adds a single-unit local search (Step 5 of
  Algorithm 2) to the best k solutions, which finds the published optimum for
  version B with 20,000 starts.
- **Rounding.** Published values have 3 decimals and appear to be rounded from
  4-decimal MATLAB output (e.g. 0.406464 → 0.4065 → 0.407); the validation
  script accepts either rounding path.

---

## Coding conventions (following the original)

- Data are stacked unit by unit: row `(i - 1) * T + t` holds unit `i`, period
  `t`, exactly as in the MATLAB and Fortran code, so that the original formulas
  translate line by line.
- Functions are named after the step of the paper they implement, and their
  comments give the equation or algorithm number and the original file.
- Figures use the paper's grayscale line styles and base graphics.
- Each script states what it reproduces, sets its expensive options at the top
  (with the paper's values in comments), and defaults to settings that run in
  seconds or minutes.

---

## License

The R translation is Copyright (C) 2026 Eliot R. Weinstein and is licensed
under the GNU General Public License, version 2 or (at your option) any later
version; see `LICENSE`.

The license covers only the translation: the R and C++ code, the scripts, the
tests, and the documentation. It grants no rights in the grouped fixed-effects
method and the original MATLAB, Fortran, and Stata code, which are Stéphane
Bonhomme and Elena Manresa's, or in the data, which come from Acemoglu,
Johnson, Robinson, and Yared (2008) via Bonhomme and Manresa's replication
files. Both remain subject to their own terms. If you use this code, please
cite Bonhomme and Manresa (2015); see `CITATION.md`.
