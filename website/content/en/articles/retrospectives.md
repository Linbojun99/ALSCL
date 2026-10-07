# Retrospectives

Complete the [worked YTF example](first-model.html) first. The examples below use its `x`, `inputs`, `a` (ACL) and `b` (ALSCL) objects.

```r
# Refit after peeling 1, 2 and 3 steps
ra <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl,
                          list(nyear = 3, train_times = 2, silent = TRUE)))
rb <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                            list(nyear = 3, train_times = 2, silent = TRUE)))
plot_retro(rb, facet_col = 2, rho_digits = 3)
rb$rho_text
```


Here `nyear` counts terminal **steps**: four quarterly steps equal one year. Retain at least two observations and keep assumptions consistent across peels. Mohn's rho measures retrospective revision, not forecast accuracy; inspect each fit's numerical status.

![ALSCL retrospective](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/retro_ALSCL.png?raw=true)

```r
# Rebuild figures, CSV and session record
source("scripts/ytf_workflow.R", encoding = "UTF-8")
result <- run_ytf_guide(out = "my_ytf_guide", npeels = 3)
# Keep live fits for TMB-dependent plots
plot_ridges(result$alscl)
```


The default output is `docs` and overwrites matching outputs; choose another folder to preserve bundled files. The workflow requires ggplot2, patchwork and jsonlite. Initial compilation needs a compiler, and retrospectives involve multiple fits. A writable TMB cache can be configured. Saved R objects do not preserve usable TMB external pointers across sessions; refit for operations requiring a live `obj`.
## Meaning of Mohn's rho

Compare each peeled estimate with the full fit at that peel's terminal period, then average relative differences:

```math
\rho=\frac{1}{K}\sum_{k=1}^{K}\frac{\hat\theta^{(-k)}_{T-k}-\hat\theta^{(0)}_{T-k}}{\hat\theta^{(0)}_{T-k}}.
```

Positive rho means shorter fits tend to be higher at their terminal periods; signs can cancel. Inspect individual trajectories and fits. A zero full-fit denominator makes the relative difference undefined.


## ACL retrospective analysis

```r
# Draw this figure
plot_retro(ra, facet_col = 2, rho_digits = 3, rho_size = 3, point_size = 1, line_size = 0.6)
```

[![ACL retrospective analysis](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/retro_ACL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/retro_ACL.png?raw=true)



Successively peel terminal observations and refit; colors indicate different terminal periods.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_retro)



## ALSCL retrospective analysis

```r
# Draw this figure
plot_retro(rb, facet_col = 2, rho_digits = 3, rho_size = 3, point_size = 1, line_size = 0.6)
```

[![ALSCL retrospective analysis](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/retro_ALSCL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/retro_ALSCL.png?raw=true)



Successively peel terminal observations and refit; colors indicate different terminal periods.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_retro)
