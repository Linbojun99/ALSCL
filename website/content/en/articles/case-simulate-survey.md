# How do I simulate survey data?

## The problem

You need controlled survey observations to test a workflow or investigate model performance. This case separates a quick input-table check from a simulation experiment that retains generating truth. Neither pathway provides empirical fishery data.

## Choose the experiment before generating data

|Question|Pathway|What to retain|
|---|---|---|
|Does import or plotting work?|Quick annual tables|Seed, settings and three input tables|
|How do estimation and generating truth differ?|Full operating model|Parameters, biological arrays, truth, observations and all fit outcomes|
|Do annual and quarterly settings behave differently?|Full presets with an explicit time step|Time units and all changed biological assumptions|

The figures in this case come from the repository’s full flatfish preset, replicate 4. They illustrate the full pathway below, not the quick-table simulator. The flatfish case has 20 retained annual observations after an 80-year burn-in; the tuna preset has 20 retained quarterly observations after a 95-year burn-in.

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

## Validate and visualize the simulated observations

The converter returns three data frames with length-bin labels in the first column. Validate them before passing them to an estimator. Use the full preset’s actual length midpoints when drawing the heatmap.

```r
source("scripts/real_data_workflow.R", encoding = "UTF-8")
validate_survey_tables(inputs, growth_step = pa$growth_step, zero_action = "error")
library(ggplot2)
years <- as.numeric(names(inputs$data.CatL)[-1])
observations <- data.frame(
  Year = rep(years, each = length(pa$len_mid)),
  Length = rep(pa$len_mid, length(years)),
  Index = as.vector(as.matrix(inputs$data.CatL[-1])))
p_data <- ggplot(observations, aes(Year, Length, fill = log10(Index))) +
  geom_tile() + scale_fill_gradient(low = "white", high = "steelblue") +
  labs(fill = "log10 index") + theme_bw()
p_data
ggsave("my_full_simulation/survey.png", p_data, width = 10, height = 4, dpi = 160)
```

![Full flatfish simulation: survey numbers by length and year](../assets/manual/figures/cases/flatfish_data.png)

**Read the figure:** color is the log10 survey index. Bands moving through length classes reflect cohort progression modified by catchability and noise. This is an observation surface, not total abundance; a change in survey catchability can change it without the same change in population size.

![Generating growth curve, maturity and survey catchability](../assets/manual/figures/cases/flatfish_biology.png)

The upper panel links age to mean length. The lower curves represent maturity and survey catchability, which have different roles. Keep their L50/L95 settings separate. The figure’s biological values are the flatfish preset, not assumptions for an unknown stock.

## Fit both structures to one replicate

For a conditional demonstration, fix the known growth and selected process parameters as in the repository’s case-study script. ACL and ALSCL use different parameter names and transformations. The following block uses the `pa` and `inputs` just generated.

```r
pA <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk), t0 = pa$t0,
           log_cv_len = log(pa$cv_L), log_std_log_F = log(pa$F_sd),
           log_std_log_N0 = log(pa$std_logN0), logit_log_R = qlogis(pa$R_ar))
pB <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk), log_t0 = log(pa$t0),
           log_cv_len = log(pa$cv_L), log_cv_grow = log(pa$cv_inc),
           log_sigma_log_F = log(pa$F_sd), log_sigma_log_N0 = log(pa$std_logN0),
           logit_log_R = qlogis(pa$R_ar))
mA <- lapply(pA, function(value) factor(NA))
mB <- lapply(pB, function(value) factor(NA))
mB$logit_log_F_l <- mB$logit_log_F_y <- factor(NA)
common <- c(inputs, list(
  rec.age = pa$rec.age, nage = pa$nage, M = pa$M,
  sel_L50 = pa$q_surv_L50, sel_L95 = pa$q_surv_L95,
  len_mid = pa$len_mid, len_border = pa$len_border[-c(1, length(pa$len_border))],
  growth_step = pa$growth_step, train_times = 2, ncores = 1, silent = TRUE))
a <- do.call(run_acl, c(common, list(parameters = pA, map = mA)))
b <- do.call(run_alscl, c(common, list(parameters = pB, map = mB)))
checks <- list(ACL = diagnose_model(inputs$data.CatL, a),
               ALSCL = diagnose_model(inputs$data.CatL, b))
```

## Compare estimates with generating truth

Here `TB` is the simulator’s total biomass and `B` is the fitted report’s biomass. Align the same periods and quantities before subtracting or dividing values.

```r
truth <- data.frame(Year = rep(years, 3),
  Quantity = rep(c("B", "SSB", "Rec"), each = length(years)),
  Truth = c(sm$TB, sm$SSB, sm$Rec))
estimates <- rbind(
  data.frame(Year = rep(years, 3), Quantity = truth$Quantity,
    Value = c(a$report$B, a$report$SSB, a$report$Rec), Model = "ACL"),
  data.frame(Year = rep(years, 3), Quantity = truth$Quantity,
    Value = c(b$report$B, b$report$SSB, b$report$Rec), Model = "ALSCL"))
p_truth <- ggplot(estimates, aes(Year, Value, color = Model)) + geom_line() +
  geom_line(data = truth, aes(Year, Truth), inherit.aes = FALSE, linetype = 2) +
  facet_wrap(~ Quantity, scales = "free_y", ncol = 1) + theme_bw()
p_truth
ggsave("my_full_simulation/truth_comparison.png", p_truth,
       width = 10, height = 6, dpi = 160)
saveRDS(list(parameters = pa, truth = sm, inputs = inputs, checks = checks,
             ACL_report = a$report, ALSCL_report = b$report),
        "my_full_simulation/assessment_summary.rds")
```

![Flatfish replicate 4: generating truth and conditional ACL/ALSCL estimates](../assets/manual/figures/cases/flatfish_truth.png)

Dashed lines are generating truth; colored lines are conditional estimates from the bundled case. Compare direction and magnitude within each panel; free y scales prevent comparisons of apparent amplitude across quantities. Because the generator here is age-based and only one replicate is shown, the plot does not rank the models in general.

## Extend to quarters or more replicates

For quarters, start with `initialize_params(species = "tuna")` and regenerate biology, observations and matching fit settings. Do not only relabel the annual columns. For repeated simulations, keep seeds, truth and failures; `sim_acl()` automates ACL fits only. ALSCL replicates need a loop around `run_alscl()` with replicate-specific inputs.

|Symptom|Check first|
|---|---|
|Time headers and `growth_step` disagree|Use the converter with the same parameter list used to generate the replicate|
|Unreasonable trajectories|Check units, burn-in, growth, recruitment age and survey catchability|
|A fit fails|Keep its error and seed; inspect numerical diagnostics before summarizing successful fits|
|Files cannot be read|Use `load()` for `sim_rep*` and `readRDS()` for files saved with `saveRDS()`|

## What to check and report

Confirm table alignment, time units, retained duration, generating biology and observation error before fitting. Two replicates demonstrate file handling; they do not establish performance. For a simulation study, vary the target assumption deliberately and retain failed fits alongside bias, RMSE and interval coverage.

Continue with the [simulation function guide](simulation.html), [annual and quarterly settings](case-studies.html), or [fitting tables and specifying biology](case-fit-survey.html).
