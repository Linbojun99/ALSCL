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

![Synthetic YTF: residual comparison for ACL and ALSCL](../assets/manual/figures/ytf/compare_residuals.png)

The panels compare residual distributions and patterns by time and length. Annual error bars describe dispersion, not confidence intervals. A model can follow biomass trends while retaining systematic observation residuals.

## Compare the fitted survey observations

Select the same years in both models and inspect the observations that constrain each fit.

```r
plot_compare_CatL(a, b, x$data.CatL,
                  years = c(2000, 2005, 2010, 2015), ncol = 2)
```

![Synthetic YTF: observed survey indices and fitted medians in four years](../assets/manual/figures/ytf/compare_CatL.png)

Points are observations and curves are fitted medians. Compare the fitted shape, especially sparse larger-length bins, before discussing differences in derived abundance. The figure uses the bundled synthetic fit.

## Compare population estimates

```r
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![Conditional population estimates from the synthetic YTF example](../assets/manual/figures/ytf/compare_ts.png)

Look for differences in timing, trend and scale. The displayed figure is a bundled illustration; your installed version and numerical settings determine the new fit. Agreement does not establish identifiability, and this single example cannot establish general superiority.

## Check the scale of fishing mortality

ACL represents fishing mortality by age; ALSCL represents it by length. The coordinates therefore differ, and individual cells in these panels do not correspond one to one.

```r
plot_compare_F(a, b)
```

![Synthetic YTF: age-based and length-based fishing mortality](../assets/manual/figures/ytf/compare_F.png)

Compare time patterns while accounting for the different state variables. If using `plot_compare_annual_F()`, its name does not mean that quarterly values are automatically annualized.

## Save a comparison that can be reviewed

```r
dir.create("my_model_comparison", showWarnings = FALSE)
write.csv(comparison$summary, "my_model_comparison/summary.csv", row.names = FALSE)
write.csv(comparison$fit_metrics, "my_model_comparison/fit_metrics.csv", row.names = FALSE)
p_population <- plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
ggplot2::ggsave("my_model_comparison/population.png", p_population,
               width = 10, height = 7, dpi = 300)
saveRDS(list(inputs = inputs, fit_args = x$fit_args, fit_config = x$fit_config,
             comparison = comparison, ACL_report = a$report, ALSCL_report = b$report),
        "my_model_comparison/summary.rds")
writeLines(capture.output(sessionInfo()), "my_model_comparison/session.txt")
```

|Observation|What it supports|What it does not establish|
|---|---|---|
|Both numerical fits are acceptable|A comparison of fitted outputs is worth examining|That all assumptions are correct|
|Residual patterns differ|The models explain observed patterns differently|That one is universally better|
|Biomass or recruitment differs|Assessment conclusions depend on population structure or associated assumptions|Which estimate is true for an empirical stock|
|Intervals overlap|The conditional estimates are not clearly separated by those intervals|Equivalence or identical management consequences|

For your own observations, replace both models’ inputs together and review model-specific parameter names. Do not compare an ACL fit to one dataset with an ALSCL fit to another while attributing all differences to model structure. Report time window, fixed biology, zero handling, numerical diagnostics and sensitivity results.

## Follow up the differences

Review fixed biological assumptions, residual structure and sensitivity to starting values. Use the [comparison plotting guide](model-comparison.html) for survey fits, mortality and growth; use [retrospective analysis](case-retrospective-bias.html) to examine revisions when terminal observations are removed.
