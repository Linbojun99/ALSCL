# How do I fit my own survey data?

## The problem

You have observed survey numbers by length and time, with mean weight and maturity information, and want to fit ALSCL. The goal is a traceable assessment workflow. No empirical dataset or biological assumptions for your stock are supplied by this case.

## Prepare the three tables

Copy the [Excel template](../assets/manual/data/ALSCL_Data_Entry_Example.xlsx) and replace its synthetic observations. Keep exactly three sheets: `data.CatL`, `data.wgt`, and `data.mat`. Each uses `LengthBin` in its first column and increasing numeric times in the remaining headers. Bins, times and ordering must match.

Record species, survey design, index standardization, length and weight units, bin boundaries and the meaning of zeros. Weight is mean individual weight; maturity is a proportion in [0,1]. See [data preparation](data-preparation.html) for missing periods and unequal bins.

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

## Preserve the assessment record

```r
saveRDS(list(inputs = my_inputs, config = config,
             diagnostics = checks, report = my_fit$report),
        "my_assessment_summary.rds")
writeLines(capture.output(sessionInfo()), "my_assessment_session.txt")
```

Retain input provenance and the package revision as well. Saved summaries are suitable for reporting; operations requiring a live TMB objective need a fit in the current R session. For a fully specified practice run before using observations, use the [synthetic YTF example](first-model.html).
