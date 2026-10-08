# How do I simulate survey data?

## The problem

You need controlled survey observations to test a workflow or investigate model performance. This case separates a quick input-table check from a simulation experiment that retains generating truth. Neither pathway provides empirical fishery data.

## Prepare the session

Install ALSCL, download or clone the repository, and run R from its root directory. The helper script below is distributed in the repository, not exported by the installed package.

```r
library(ALSCL)
source("scripts/simulation_workflow.R", encoding = "UTF-8")
```

## Path 1: create annual input tables

Use this path to test import, validation and plotting. Choose a new output directory if you want to preserve previous files.

```r
quick <- simulate_example_data(
  years = 2000:2019, seed = 42, bin_breaks = seq(5, 51, 2),
  Linf = 60, vbk = 0.2, t0 = 1/60, M = 0.2, nage = 15,
  L50_sel = 15, L95_sel = 20, L50_mat = 35, L95_mat = 40,
  wgt_a = exp(-12), wgt_b = 3, mean_F = 0.3,
  cv_catch = 0.2, rec_sigma = 0.3,
  save_csv = TRUE, output_dir = "my_annual_simulation")
str(quick)
```

Inspect the returned tables and exported CSV files. This simplified annual simulator differs from the full operating model. Relabeling years as quarters does not create quarterly dynamics.

## Path 2: retain truth and simulation replicates

Use this path when the question concerns bias, uncertainty or differences between an operating model and an estimator.

```r
pa <- initialize_params(species = "flatfish", observation_error = "independent")
bio <- sim_cal(pa)
sm <- sim_data(bio, pa, iter_range = 4:5, return_iter = 4,
               output_dir = "my_full_simulation")
inputs <- simulation_to_tables(sm, pa)
saveRDS(list(parameters = pa, truth = sm, inputs = inputs),
        "my_full_simulation/truth_and_inputs.rds")
str(inputs)
```

`iter_range` controls replicate IDs and seeds; `return_iter` selects the in-memory replicate. `simulation_to_tables()` transposes survey numbers and constructs time headers from `growth_step`. The saved `sim_rep*` files use `load()`, not `readRDS()`.

## What to check and report

Confirm table alignment, time units, retained duration, generating biology and observation error before fitting. Two replicates demonstrate file handling; they do not establish performance. For a simulation study, vary the target assumption deliberately and retain failed fits alongside bias, RMSE and interval coverage.

Continue with the [simulation function guide](simulation.html), [annual and quarterly settings](case-studies.html), or [fitting tables and specifying biology](case-fit-survey.html).
