# How do I fit my own survey data?

## The problem

You have observed survey numbers by length and time, with mean weight and maturity information, and want to fit ALSCL. The goal is a traceable assessment workflow. No empirical dataset or biological assumptions for your stock are supplied by this case.

## Prepare the three tables

Copy the [Excel template](../assets/manual/data/ALSCL_Data_Entry_Example.xlsx) and replace its synthetic observations. Keep exactly three sheets: `data.CatL`, `data.wgt`, and `data.mat`. Each uses `LengthBin` in its first column and increasing numeric times in the remaining headers. Bins, times and ordering must match.

Record species, survey design, index standardization, length and weight units, bin boundaries and the meaning of zeros. Weight is mean individual weight; maturity is a proportion in [0,1]. See [data preparation](data-preparation.html) for missing periods and unequal bins.

## What the workbook should look like

The following screenshots show the supplied **synthetic YTF workbook**, illustrating layout only. Replace its values with your observations and retain a separate copy of the original template.

![Survey-number input sheet](../assets/manual/figures/excel_CatL.png)

Enter a consistently standardized number index. Blank/NA means missing; do not recode missing observations as zero.

![Mean individual weight input sheet](../assets/manual/figures/excel_wgt.png)

Enter mean individual weight for each cell, not total catch weight. Document the weight unit and any length–weight relationship used.

![Maturity-proportion input sheet](../assets/manual/figures/excel_mat.png)

Enter proportions between zero and one. All three sheets must use exactly the same bins, periods and order.

## Import and validate

Install `readxl` once, and run from the downloaded repository root so the guide helpers can be sourced. Replace the workbook path with your own file.

```r
library(ALSCL)
source("scripts/real_data_workflow.R", encoding = "UTF-8")
my_inputs <- read_survey_excel("my_survey.xlsx")
step_years <- 1
validate_survey_tables(my_inputs, growth_step = step_years, zero_action = "error")
```

Use `step_years <- 0.25` only for a regular quarterly series and consistent biological units. Validation stops on zeros by default. The package default of excluding zero observations is not automatically appropriate for genuine zero catches; a zero-inflated observation model is not implemented.

## First run the supplied workbook end to end

This optional practice run is fully specified. It uses the template's synthetic YTF data and **its matching biology**, so it is useful for checking installation, import and compilation before replacing any observations. It does not estimate your stock.

```r
library(ALSCL)
source("scripts/real_data_workflow.R", encoding = "UTF-8")
data("YTF_example")
demo_inputs <- read_survey_excel("docs/data/ALSCL_Data_Entry_Example.xlsx")
demo_config <- c(YTF_example$fit_args, YTF_example$fit_config$alscl,
                 list(zero_action = "error", train_times = 2, silent = TRUE))
demo_fit <- fit_survey_tables(demo_inputs, demo_config, model_type = "alscl")
diagnose_model(demo_inputs$data.CatL, demo_fit)
plot_CatL(demo_fit, type = "year", exp_transform = FALSE, facet_ncol = 4)
```

![Synthetic YTF illustration: observed and fitted survey index by period](../assets/manual/figures/ytf/CatL_year_FALSE.png)

Compare fitted patterns with observations on the displayed log scale. Systematic under- or overprediction within a range of lengths suggests a structural issue worth investigating. This bundled figure illustrates the template run; it is not a result for your workbook.

## Specify the assessment assumptions

The following is an editable configuration template. Replace every `NA` using evidence for your stock before running the fit. The guard deliberately stops incomplete configurations.

```r
config <- list(
  rec.age = NA_real_, nage = NA_integer_, M = NA_real_,
  sel_L50 = NA_real_, sel_L95 = NA_real_,
  growth_step = step_years, zero_action = "error",
  train_times = 2, silent = TRUE)
required_biology <- c("rec.age", "nage", "M", "sel_L50", "sel_L95")
stopifnot(all(is.finite(unlist(config[required_biology]))))
```

This list is a starting template, not a complete biological specification. Also review actual length boundaries, growth assumptions, initial values, bounds and fixed parameters. Add the corresponding named fit arguments to `config`; see [fitting controls](fitting-controls.html), [estimation parameters](estimation-parameters.html), and [run_alscl()](../reference/run_alscl.html). Do not copy YTF growth or catchability settings solely because your tables have similar dimensions. Mortality must match the model time step.

### Review these inputs before fitting observations

|Input|What must be justified|
|---|---|
|`rec.age`, `nage`|Recruitment age and number of modeled age classes|
|`M`|Natural mortality per model step, with its source and uncertainty|
|`sel_L50`, `sel_L95`|Fixed survey catchability locations; these are not estimated fishery selectivity|
|`growth_step`|Years per step: e.g. 1 or 0.25, matching the time headers|
|`len_mid`, `len_border`|Actual midpoints and internal boundaries; unequal bins need explicit interpretation|
|`parameters`, `parameters.L`, `parameters.U`|Starting values and bounds on the optimizer's scale|
|`map`|Which parameters are fixed, and why; intervals are conditional on fixed settings|

For ALSCL, growth starts use `log_Linf`, `log_vbk` and `log_t0`. This parameterization restricts `t0` to positive values; do not take a logarithm of a negative external estimate. Review model suitability rather than silently changing its sign. An initial value does not fix a parameter: fixing requires a corresponding map entry such as `factor(NA)`.

`train_times` counts successive optimization passes, not independent random starts. Different starting values should be assessed separately; see the [fitting guide](fitting-controls.html).

## Fit and inspect

After completing and reviewing `config`, fit the model and check numerical diagnostics before interpreting population estimates.

```r
my_fit <- fit_survey_tables(my_inputs, config, model_type = "alscl")
checks <- diagnose_model(my_inputs$data.CatL, my_fit)
checks
plot_residuals(my_fit, type = "year")
plot_biomass(my_fit)
```

Inspect optimizer status, gradient, Hessian, parameter boundaries and residual patterns separately. A completed fit is not evidence that all parameters are identifiable. Reconsider assumptions when diagnostics or sensitivity checks are poor; [diagnostics](diagnostics.html) explains the tools.

## Read the diagnostic and population plots

These two figures use the bundled YTF fit to explain the display. Run the commands on `my_fit` to obtain the corresponding figures for your own observations.

![Synthetic YTF illustration: residual patterns by year](../assets/manual/figures/ytf/residuals_year.png)

Look for sustained positive or negative residuals and changes in spread. These are raw log residuals; visible structure should lead to checks of input standardization, catchability assumptions and process structure, not automatic deletion of observations.

![Synthetic YTF illustration: biomass trajectory and conditional uncertainty](../assets/manual/figures/ytf/plot_biomass_B.png)

Interpret the trajectory with the supplied weight units and fixed assumptions. The plotted interval does not include uncertainty in every external biological input. Apparent precision cannot replace sensitivity analysis.

### Export your own results

```r
dir.create("my_assessment", showWarnings = FALSE)
write.csv(checks, "my_assessment/diagnostics.csv", row.names = FALSE)
p_index <- plot_CatL(my_fit, type = "year", facet_ncol = 4)
p_residual <- plot_residuals(my_fit, type = "year")
p_biomass <- plot_biomass(my_fit)
ggplot2::ggsave("my_assessment/survey_fit.png", p_index, width = 10, height = 7, dpi = 300)
ggplot2::ggsave("my_assessment/residuals.png", p_residual, width = 10, height = 6, dpi = 300)
ggplot2::ggsave("my_assessment/biomass.png", p_biomass, width = 8, height = 5, dpi = 300)
```

### If something fails

|Symptom|Next action|
|---|---|
|Workbook sheet or header error|Restore exact sheet names, `LengthBin`, and numeric time headers|
|Tables differ|Align bins and time columns in all three sheets|
|Zero observations rejected|Establish whether zeros are genuine; do not replace them with an arbitrary small constant|
|Compiler error|Follow the platform-specific installation guide and check the writable TMB cache|
|Non-positive-definite Hessian or large gradient|Inspect starts, bounds, fixed assumptions and identifiability before interpreting uncertainty|
|Large changes under different M or q|Report sensitivity to supported external values rather than selecting the most convenient result|

## Preserve the assessment record

```r
saveRDS(list(inputs = my_inputs, config = config,
             diagnostics = checks, report = my_fit$report),
        "my_assessment_summary.rds")
writeLines(capture.output(sessionInfo()), "my_assessment_session.txt")
```

Retain input provenance and the package revision as well. Saved summaries are suitable for reporting; operations requiring a live TMB objective need a fit in the current R session. For a fully specified practice run before using observations, use the [synthetic YTF example](first-model.html).
