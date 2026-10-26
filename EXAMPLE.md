# Example: watching the estimator work

R translation by  
Eliot R. Weinstein  
eliotw@uchicago.edu  
[erweinstein.com](https://erweinstein.com)

A good first exercise with this code is to change two things, the number of
groups G and the starting values of the search, and watch how the estimated
grouping and the estimates of theta move. It takes a few minutes and shows
quickly why Bonhomme and Manresa (2015) put so much care into computation.

## Why this is worth doing

The GFE objective is minimized over all G^N ways of assigning N units to G
groups, so no search is guaranteed to find the global minimum. Bonhomme and
Manresa's answer is a model of careful applied econometrics. They define the
estimator as a single optimization problem, then check the computation instead
of assuming it works:

- they run the heuristics from many starting values, with an explicit rule that
  different starts should reach the exact same solution;
- they certify the heuristics with an exact branch-and-bound algorithm on
  subsamples;
- they pair the asymptotic theory with Monte Carlo simulations calibrated to
  their own data;
- they organize the robustness checks around the main finding, asking in each
  one whether the grouped patterns survive a specific objection;
- and they release replication files that let a reader rerun every step.

The exercises below let you see, on their data, why each of those steps
matters.

## Setup

Open `gfe_bm_R.Rproj` in RStudio, then:

```r
source("gfe.R")
p   <- load_bm_panel()                    # 90 countries, 1970-2000
pub <- load_bm_published()$groupings      # the authors' optimal groupings
cum <- function(th) th[2] / (1 - th[1])   # cumulative income effect
```

## 1. Same G, different starting values

Run Algorithm 1 from 200 random starting values with G = 3, one start at a
time, and keep the objective and theta from each run.

```r
set.seed(2015)
runs <- t(replicate(200, {
  f <- gfe(p, G = 3, algorithm = "alg1", n_starts = 1)
  c(objective = f$objective, theta1 = f$theta[[1]], theta2 = f$theta[[2]],
    cumulative = cum(f$theta))
}))
runs <- as.data.frame(runs)

summary(runs$objective)
sum(runs$objective < 16.5988)   # starts that reach the global minimum, 16.59873

plot(runs$objective, runs$theta1, pch = 16, col = gray(0, 0.4),
     xlab = "objective at convergence (G = 3)",
     ylab = "estimate of theta1 (lagged democracy)")
abline(v = 16.59873, lty = 2)   # global minimum (published solution)
```

Each point is a local minimum. With this seed, none of the 200 single starts
reaches the global minimum (16.5987, the published solution; with other seeds a
few do), and the local minima give quite different estimates of theta. Taking the first answer the algorithm returns would be a
mistake. Taking the best of many starts, and checking that the best is reached
repeatedly, is what the paper recommends.

Variable neighborhood search (Algorithm 2) adds a local search and random
jumps to escape these local minima. Compare a handful of single-start VNS runs:

```r
set.seed(2015)
vns_runs <- sapply(1:5, function(s) gfe(p, 3, algorithm = "vns", n_starts = 1,
                                        neighmax = 10, itermax = 10)$objective)
round(vns_runs, 5)   # most or all reach 16.59873
```

## 2. Same data, different G

Now hold the search fixed and change G. With the authors' groupings, this is
instant:

```r
fits <- lapply(2:6, function(G) gfe(p, G, groups = pub[[paste0("G", G)]]))
data.frame(G          = 2:6,
           objective  = sapply(fits, `[[`, "objective"),
           theta1     = sapply(fits, function(f) f$theta[[1]]),
           theta2     = sapply(fits, function(f) f$theta[[2]]),
           cumulative = sapply(fits, function(f) cum(f$theta)))
```

The objective always falls as G grows; that alone says nothing about the right
G (see `gfe_bic()`). The coefficient on lagged democracy moves a lot between
G = 2 and G = 4, because the group-specific time effects absorb persistence
that a model with fewer groups attributes to lagged democracy. The cumulative
income effect is more stable.

To see how the groups themselves change, cross-tabulate the groupings for
neighboring values of G:

```r
table(G3 = pub$G3, G4 = pub$G4)
table(G4 = pub$G4, G5 = pub$G5)
```

Most countries stay together: each G = 3 group maps mainly onto one or two of
the G = 4 groups, and only a handful of countries move elsewhere. Group labels
are arbitrary, so focus on the large cells, not on the numbers on the
margins.

## 3. Put the two together

Finally, search from scratch for several G, and check each search against the
published solution:

```r
set.seed(1)
for (G in 2:4) {
  v <- gfe(p, G, algorithm = "vns", n_starts = 5, neighmax = 10, itermax = 5)
  m <- match_groups(v$groups, pub[[paste0("G", G)]], G)
  cat(sprintf("G = %d: objective %.5f, %d of 5 starts reach it, %d countries classified differently\n",
              G, v$objective, sum(abs(v$objectives_by_start - v$objective) < 1e-8),
              round(m$misclassification * p$N)))
}
```

As G grows, each search takes longer and fewer starts tend to reach the best
solution. With 10 starts each, `scripts/01_validate_against_published.R` finds
9, 8, 4, and 3 of 10 reaching the optimum for G = 2, 3, 4, and 5.
This is why the paper uses (10; 10; 10) for its main estimates and an exact
algorithm to check the heuristics where that is feasible
(`scripts/03_exact_algorithm.R`).

## Going further

- `GFE_walkthrough.qmd` goes through the paper section by section.
- `scripts/04_montecarlo.R` repeats the exercise above many times on simulated
  data, where the true grouping is known, and shows what misclassification
  does to inference.
- To try the estimator on your own panel, build it with
  `gfe_panel_from_df(df, id, time, y, x)` and call `gfe()` on the result.

## Reference

Bonhomme, Stéphane and Elena Manresa (2015), "Grouped Patterns of Heterogeneity
in Panel Data," *Econometrica* 83(3), 1147–1184. doi:10.3982/ECTA11319
