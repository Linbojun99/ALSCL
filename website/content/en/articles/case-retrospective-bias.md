# Do estimates change when recent data are removed?

## The problem

An assessment may revise historical biomass or recruitment as new observations arrive. A retrospective analysis asks how estimates differ when the terminal observations are removed and the same model is refitted. It does not measure forecast accuracy.

## Refit shorter series

This complete setup uses the synthetic YTF example and its conditional biological settings. It requires a working compiler and performs multiple fits.

```r
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
rb <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                            list(nyear = 3, train_times = 2, silent = TRUE)))
plot_retro(rb, facet_col = 2, rho_digits = 3)
rb$rho_text
```

## Inspect trajectories and Mohn’s rho

![Bundled ALSCL retrospective example](../assets/manual/figures/ytf/retro_ALSCL.png)

Compare each peel with the full fit at the peel’s final period. Positive rho means shortened fits tend to be higher at those terminal periods; negative rho indicates the opposite. Opposing deviations can cancel, so inspect individual curves and the numerical status of each fit. A zero reference estimate makes its relative difference undefined.

## Compare the ACL retrospective

Use the same observations and peel count, but keep ACL’s own parameter names and maps.

```r
ra <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl,
                          list(nyear = 3, train_times = 2, silent = TRUE)))
plot_retro(ra, facet_col = 2, rho_digits = 3)
ra$rho_text
```

![Synthetic YTF: ACL retrospective trajectories](../assets/manual/figures/ytf/retro_ACL.png)

The two model figures are bundled examples. Compare revisions within each model before comparing rho between models; inspect the same quantity and terminal periods. Different y scales can exaggerate or hide apparent differences.

## Inspect the returned tables

```r
head(rb$results)
rb$last_points
rb$rho_text
names(rb)
```

`results` stores trajectories, `last_points` identifies terminal estimates, and `rho_text` stores the quantity-specific average relative revision. The wrapper does **not** return the underlying fitted model for every peel or a complete table of numerical diagnostics.

### Optional: audit each ALSCL peel explicitly

This additional loop refits the full series and three peels so that numerical status and failures can be recorded. It duplicates the fitting work above; use it when a diagnostic audit is needed. Keep the label column and remove the same terminal time columns from all three tables.

```r
peel_checks <- lapply(0:3, function(k) {
  keep <- seq_len(ncol(inputs$data.CatL) - k)
  peeled <- lapply(inputs, function(d) d[, keep, drop = FALSE])
  tryCatch({
    fit <- do.call(run_alscl, c(peeled, x$fit_args, x$fit_config$alscl,
                               list(train_times = 2, silent = TRUE)))
    data.frame(Peel = k, Terminal = max(fit$year),
      Code = fit$convergence_code, Gradient = fit$max_abs_gradient,
      pdHess = fit$pdHess, Boundary = fit$bound_hit, Error = "")
  }, error = function(e) data.frame(Peel = k,
    Terminal = as.numeric(tail(names(peeled$data.CatL), 1)),
    Code = NA_integer_, Gradient = NA_real_, pdHess = NA,
    Boundary = NA, Error = conditionMessage(e)))
})
peel_checks <- do.call(rbind, peel_checks)
peel_checks
```

Do not average an unreliable peel into a scientific interpretation merely because a curve was returned. Report failures and reconsider the fitting assumptions. There is no single gradient or rho threshold that guarantees a valid assessment.

## Adapt the analysis to your assessment

Replace the inputs and biology with the same configuration used for your full-data fit. Keep fixed parameters, observation handling and time units consistent. `nyear` counts removed terminal steps: three quarterly steps are not three years. Retain sufficient observations in each peel; see [the function guide](retrospectives.html).

## Export the retrospective record

```r
dir.create("my_retrospective", showWarnings = FALSE)
write.csv(rb$results, "my_retrospective/alscl_trajectories.csv", row.names = FALSE)
write.csv(rb$rho_text, "my_retrospective/alscl_rho.csv", row.names = FALSE)
write.csv(ra$rho_text, "my_retrospective/acl_rho.csv", row.names = FALSE)
p_retro <- plot_retro(rb, facet_col = 2, rho_digits = 3)
ggplot2::ggsave("my_retrospective/alscl.png", p_retro,
               width = 10, height = 7, dpi = 300)
saveRDS(list(ACL = ra, ALSCL = rb, fit_args = x$fit_args,
             fit_config = x$fit_config), "my_retrospective/results.rds")
writeLines(capture.output(sessionInfo()), "my_retrospective/session.txt")
```

|Symptom|Interpretation or check|
|---|---|
|Too few time steps after peeling|Reduce peel count; every fit must retain at least two time steps|
|Opposite-signed terminal deviations|Inspect individual peels; average rho can hide cancellation|
|Quarterly analysis interpreted as years|Convert removed steps to duration using `growth_step`|
|A custom `zero_action` is required|The retrospective wrapper does not expose it; use an explicit loop with the fit functions and retain the same observation treatment|
|A short fit is unstable|Inspect its diagnostics before interpreting its terminal revision|

## Decide what to investigate next

A repeated directional revision motivates checks of input changes, model assumptions and parameter sensitivity. Rho alone does not identify a cause or supply a correction factor. Record the full fit, each peel’s diagnostics and its terminal period before drawing conclusions.
