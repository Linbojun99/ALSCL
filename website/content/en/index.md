> Age- and length-structured catch-at-length models for fisheries stock assessment

**ALSCL** is an R package for fitting population models to fishery-independent survey length data using [Template Model Builder (TMB)](https://github.com/kaskr/adcomp). It provides an age-structured model (**ACL**) and a joint age–length model (**ALSCL**), together with simulation, diagnostics, model comparison and plotting tools.

Start with a [worked example](articles/first-model.html), browse the [basic functions](#contents) and [worked case studies](articles/index.html), or look up a function in the [reference](reference/index.html). The model framework is based on [Zhang & Cadigan (2022)](https://doi.org/10.1111/faf.12673).

<h2 id="contents">Contents</h2>

- [Installation](#installation)
- [Two model structures](#models)
- [Basic use](#basic-use)
- [Check and interpret the fit](#check-fit)
<!-- ARTICLE_CONTENTS -->
- [Getting help](#help)
- [Citation](#citation)

<h2 id="installation">Installation</h2>

Install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
```

Fitting compiles the bundled TMB models, so an R-compatible C++ toolchain is required. See [installation and compilation](articles/getting-started.html) for platform-specific guidance.

<h2 id="models">Two model structures</h2>

| | ACL | ALSCL |
|---|---|---|
| Population state | Numbers by age and time | Numbers by length, age and time |
| Growth representation | Age–length probability matrix | Growth transition matrix |
| Fishing mortality | Estimated by age | Estimated by length; age-specific F is derived |
| Fitting function | `run_acl()` | `run_alscl()` |

Both models use three aligned input tables: survey numbers-at-length, mean individual weight-at-length, and maturity-at-length. Natural mortality and survey catchability settings are supplied. Their assumptions affect the scale and interpretation of the estimates; see the [model framework](articles/model-theory.html).

<h2 id="basic-use">Basic use</h2>

The bundled `YTF_example` contains synthetic yellowtail-flounder-style data, biological settings and known simulation truth. It is a teaching dataset, not the empirical dataset from the paper.

```r
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]

# Fit both models to the same observations and biological settings
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))

plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![ACL and ALSCL biomass, spawning biomass and recruitment estimates from the synthetic YTF example](assets/manual/figures/ytf/compare_ts.png)

The example fixes selected biological and process parameters. Its uncertainty intervals are conditional on those settings. The [worked example](articles/first-model.html) explains the inputs, fitted objects and data provenance.

<h2 id="check-fit">Check and interpret the fit</h2>

```r
diagnose_model(x$data.CatL, b)
plot_residuals(b, type = "year")
comparison <- compare_models(a, b, x$data.CatL)
comparison$summary
```

Check the optimizer status, maximum absolute gradient, Hessian and parameter boundaries separately. Then inspect residual patterns and sensitivity to fixed assumptions. A successful optimizer message alone does not establish reliable inference. See [diagnostics](articles/diagnostics.html) and [model comparison](articles/model-comparison.html).

<h2 id="learn-more">Basic functions</h2>

<!-- ARTICLE_GUIDE -->

<h2 id="help">Getting help</h2>

Use the [GitHub issue tracker](https://github.com/Linbojun99/ALSCL/issues) for bugs and feature requests. Include a small reproducible example, the package version, input dimensions and the output of `sessionInfo()`. Check [troubleshooting](articles/reproducibility.html) for common data and compilation problems.

<h2 id="citation">Citation</h2>

Zhang, F. & Cadigan, N. G. (2022). An age- and length-structured statistical catch-at-length model for hard-to-age fisheries stocks. *Fish and Fisheries*, **23**(5), 1121–1135. [doi:10.1111/faf.12673](https://doi.org/10.1111/faf.12673).

<h3 id="krill-paper">Related application: Antarctic krill</h3>

Dong, S., Zhang, F. & Zhu, G. (2025). Length-dependent growth and mortality within each cohort cannot be ignored in stock assessment: a case study of Antarctic krill *Euphausia superba*. *Marine Ecology Progress Series*, **769**, 1–22. [doi:10.3354/meps14923](https://doi.org/10.3354/meps14923).

This study applies ACL and ALSCL to Antarctic krill to examine how length-dependent growth and mortality within cohorts affect assessment results.

When reporting analyses, also record the ALSCL version, code revision, data source and fitted assumptions.
