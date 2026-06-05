# NMABayesianalltypes-

A single-file R **Shiny** application for **Bayesian network meta-analysis** (and
pairwise meta-analysis) across multiple effect-measure types.

## What it does

The app (`app.R`) fits Bayesian meta-analysis models via JAGS (through `rjags`)
and supports both:

- **Pairwise** meta-analysis (random-effects hierarchical normal model), and
- **Network meta-analysis** (contrast-based model with one mean parameter per
  treatment, optional **meta-regression**).

Supported summary measures: **SMD, MD, log Odds Ratio, log Risk Ratio,
log Hazard Ratio** (effects/SEs supplied on the analysis scale).

It also constructs a frequentist `netmeta` object for the network branch and a
`meta` (`metagen`) object for the pairwise branch, and provides:

- Sample CSVs (pairwise and network, one per summary measure) as downloads.
- Plots: traceplot, density, autocorrelation, Gelman, forest, funnel,
  caterpillar, residual, bubble, pairs, posterior predictive check, network
  plot, ECDF, netrank table, netheat, netsplit forest.
- Convergence diagnostics (Gelman-Rubin, Geweke, effective sample size),
  study-level and heterogeneity summaries, and ECDF threshold probabilities.

UI is built with `bs4Dash`; MCMC handling uses `coda`.

## Prerequisites

### 1) JAGS (external system binary — required for the Bayesian engine)

The Bayesian models run through `rjags`, which is an R interface to **JAGS**
(Just Another Gibbs Sampler). JAGS is a separate program that must be installed
on your operating system **before** installing the `rjags` R package:

- Download/install JAGS: https://mcmc-jags.sourceforge.io/

Without the JAGS binary, `rjags` will not load and the Bayesian analyses
(pairwise and network) cannot run. The frequentist `netmeta` path does not
require JAGS.

### 2) R packages

```r
install.packages(c(
  "shiny", "bs4Dash", "rjags", "coda", "ggplot2",
  "igraph", "netmeta", "dmetar", "meta"
))
```

(`grid` is part of base R and does not need installing. `dmetar` may need to be
installed from its GitHub source on some systems.)

## Run

From the repository root:

```r
shiny::runApp(".")
```

## Files

- `app.R` — the Shiny application (UI + server + `shinyApp()`).
- `test_nma.R` — a stat-core smoke test that runs a small known network through
  `netmeta::netmeta()` with `stopifnot()` property checks (frequentist path
  only; does not exercise JAGS/rjags). Run with:
  `Rscript test_nma.R`
- `CC ZERO LICENSE` — license text.

## Notes

- The main file was previously named `NMA Bayseian all types` (with spaces and a
  typo); it has been renamed to `app.R` to follow the Shiny single-file
  convention and remove the access problem the spaced name caused.
- This environment cannot run the full app (JAGS/`rjags`, `shiny`, and `bs4Dash`
  are not installed here), so only `app.R` parsing and the frequentist `netmeta`
  smoke test were verified locally.
