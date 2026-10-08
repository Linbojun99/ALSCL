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

## Adapt the analysis to your assessment

Replace the inputs and biology with the same configuration used for your full-data fit. Keep fixed parameters, observation handling and time units consistent. `nyear` counts removed terminal steps: three quarterly steps are not three years. Retain sufficient observations in each peel; see [the function guide](retrospectives.html).

## Decide what to investigate next

A repeated directional revision motivates checks of input changes, model assumptions and parameter sensitivity. Rho alone does not identify a cause or supply a correction factor. Record the full fit, each peel’s diagnostics and its terminal period before drawing conclusions.
