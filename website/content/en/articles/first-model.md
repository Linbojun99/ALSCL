# Built-in YTF

## Fit the example

Run this code in one R session after [installation](getting-started.html). Both calls use the same observations but different population structures. They create `a` (ACL), `b` (ALSCL), `x` and `inputs`, reused throughout the tutorials.

```r
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
# Inspect dimensions and provenance
str(inputs)
x$provenance
# Same observations and conditional settings
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![Model comparison](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)

|Object|Contents and use|
|---|---|
| `YTF` |Original 26-element biological parameter list, no survey table|
| `YTF_example` |23 bins × 20 periods, three tables, parameters, truth, fit settings and provenance|
| `example_data` |Legacy ALSCL tables, 9 bins × 21 years, provenance unverified|


`YTF_example` reuses shared YTF parameters, including survey log SD 0.1, with remaining flatfish defaults. Seed 4, 80 burn-in years and 20 retained years are used. Dates are illustrative; cm and kg are teaching conventions. These are not empirical paper data. The data-raw script rebuilds the object and restores the caller's RNG state.


`x$fit_config` fixes known growth and selected process/correlation parameters. Intervals are conditional on these assumptions and omit their uncertainty. The example does not establish identifiability when all parameters are estimated together.

## Inspect the result

```r
names(b)
names(b$report)
diagnose_model(x$data.CatL, b)
plot_residuals(b, type = "year")
```

`report` contains derived population quantities and survey predictions. `opt` contains optimization results. `est_std` and `vcov` describe approximate uncertainty; `obj` is the live TMB objective. Preserve the original settings when comparing models or running retrospectives. See [diagnostics](diagnostics.html) before interpreting the estimates.
