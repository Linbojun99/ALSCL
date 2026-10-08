# Do ACL and ALSCL give different answers?

## The problem

You want to determine whether changing population structure changes the assessment. Fit ACL and ALSCL to the same synthetic YTF observations while preserving the supplied biological settings. This is a conditional comparison, not a test that one model is universally better.

## Fit the same observations

The bundled data are simulated. Known growth and selected process parameters are fixed by the example configuration; their uncertainty is omitted from the fitted intervals.

```r
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))
```

## Check whether the fits are interpretable

```r
diagnose_model(x$data.CatL, a)
diagnose_model(x$data.CatL, b)
comparison <- compare_models(a, b, x$data.CatL)
comparison$summary
plot_compare_residuals(a, b, x$data.CatL)
```

Inspect numerical status and residual patterns before comparing fitted trajectories. Compare likelihood-based criteria only when the observation data and likelihood definitions are comparable. Do not rank a numerically unreliable fit by its objective value alone.

## Compare population estimates

```r
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![Conditional population estimates from the synthetic YTF example](../assets/manual/figures/ytf/compare_ts.png)

Look for differences in timing, trend and scale. The displayed figure is a bundled illustration; your installed version and numerical settings determine the new fit. Agreement does not establish identifiability, and this single example cannot establish general superiority.

## Follow up the differences

Review fixed biological assumptions, residual structure and sensitivity to starting values. Use the [comparison plotting guide](model-comparison.html) for survey fits, mortality and growth; use [retrospective analysis](case-retrospective-bias.html) to examine revisions when terminal observations are removed.
