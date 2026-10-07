# Annual and quarterly cases

Both cases are refitted with the current package to teach annual and quarterly workflows. They do not reproduce the paper's data, and one replicate cannot rank the models. General plotting lessons continue to use `YTF_example`.

|Setting| Flatfish | Tuna |
|---|---|---|
|Generator| age_based | length_based |
|Step, years| 1 | 0.25 |
| Total years / burn-in years | 100 / 80 | 100 / 95 |
|Retained steps| 20 | 20 |
| rec.age / nage | 1 / 15 | 0.25 / 20 |
| Linf / annual k | 60 / 0.20 | 152 / 0.38 |
| t0, years | 1/60 | 1/152 |
| M / mean F, per step | 0.20 / 0.30 | 0.20 / 0.20 |
| Maturity L50 / L95 | 35 / 40 | 100 / 120 |
| Survey catchability L50 / L95 | 15 / 20 | 30 / 50 |
|Midpoints| 6, 8, …, 50 | 12.5, 17.5, …, 122.5 |
|Survey log SD| 0.20 | 0.10 |


The flatfish preset uses survey SD 0.20, while `YTF_example` retains YTF's 0.10. Tuna M is per quarter. Units must match the length-weight coefficients; dates starting in 2000 are illustrative.

## Read and rebuild

```r
# Read the supplied cases
source("scripts/real_data_workflow.R", encoding="UTF-8")
flatfish_inputs <- read_survey_csv("docs/data/cases/flatfish")
tuna_inputs <- read_survey_csv("docs/data/cases/tuna")
validate_survey_tables(flatfish_inputs, growth_step=1)
validate_survey_tables(tuna_inputs, growth_step=.25)
pa <- readRDS("docs/data/cases/tuna/truth.rds")$parameters
str(pa)
# Regenerate inputs
source("scripts/simulation_workflow.R", encoding="UTF-8")
simulations <- run_simulation_examples(out="my_simulations", fit_batch=FALSE)
# Refit both cases and regenerate seven figures
source("scripts/case_studies.R", encoding="UTF-8")
case_status <- run_case_studies(out="docs")
```


The final call reads saved truth from `out/data/cases` and writes results under `out`. Copy these inputs first when choosing another directory. The linked R script defines every starting value, map and biological argument.

## Annual flatfish

![Synthetic survey](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_data.png?raw=true)


Colors encode log10 survey numbers. Cohort patterns are modified by catchability and observation noise.

![Biology](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_biology.png?raw=true)


The growth curve links age to mean length. Maturity and survey q use different L50/L95 values and are not interchangeable.

![Truth and estimates](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_truth.png?raw=true)


Black dashed curves are generating truth; colored curves are conditional estimates. Compare trajectories and offsets without treating agreement as proof of identifiability.

## Quarterly tuna

![Synthetic survey](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_data.png?raw=true)


Twenty retained quarters span five years, not twenty. The preset changes dynamics, bins, ages and biological parameters together.

![Biology](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_biology.png?raw=true)

![Truth and estimates](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_truth.png?raw=true)


Tuna is generated with length-based dynamics and fitted with an explicit quarterly step. Known growth, CV and selected process/correlation parameters are fixed; conditional uncertainty omits errors in those assumptions.

## Extend into a simulation study


Use multiple seeds and vary observation noise, M, q, binning or series length systematically. Retain failed fits before summarizing bias, RMSE, coverage and retrospectives; averaging only successful fits introduces selection bias.

[Four-fit convergence record](https://github.com/Linbojun99/ALSCL/blob/main/docs/results/case_convergence.csv)。
