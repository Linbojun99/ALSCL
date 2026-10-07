# Simulation

Load the package before running these examples.

```r
library(ALSCL)
```

### Quick input-table simulator

```r
# Simple annual simulation with CSV export
quick <- simulate_example_data(
  years = 2000:2019, seed = 42, bin_breaks = seq(5, 51, 2),
  Linf = 60, vbk = 0.2, t0 = 1/60, M = 0.2, nage = 15,
  L50_sel = 15, L95_sel = 20, L50_mat = 35, L95_mat = 40,
  wgt_a = exp(-12), wgt_b = 3, mean_F = 0.3,
  cv_catch = 0.2, rec_sigma = 0.3,
  save_csv = TRUE, output_dir = "my_simulation")
str(quick)
```


Use consecutive annual `years`. Changing labels does not create quarterly dynamics. `cv_catch` controls survey observation CV and `rec_sigma` log recruitment variation. This simplified simulator uses common F, independent recruitment perturbations and fixed length CV 0.1; it differs from both the full operating model and the bundled YTF generation recipe.

### Full simulator and two cases

The seed-4 data, input CSV files and complete truth objects are included in the repository. See [download and reading instructions](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/cases/README.md).

```r
# Parameters → biological arrays → simulation
pa <- initialize_params(species = "flatfish", observation_error = "independent")
bio <- sim_cal(pa)
sm <- sim_data(bio, pa, iter_range = 4:5, return_iter = 4,
               output_dir = "simulation_examples/flatfish")
# Length-based quarterly case
tuna_pa <- initialize_params(species = "tuna")
tuna <- sim_data(sim_cal(tuna_pa), tuna_pa, iter_range = 4:5,
                 return_iter = 4, output_dir = "simulation_examples/tuna")
# Convert simulated matrices to input tables
source("scripts/simulation_workflow.R", encoding = "UTF-8")
flatfish_inputs <- simulation_to_tables(sm, pa)
tuna_inputs <- simulation_to_tables(tuna, tuna_pa)
```

|Setting| flatfish / YTF-style | tuna |
|---|---|---|
|Operating model|Age-based|Joint age-length|
| `growth_step` | 1 year | 0.25 year |
| `nyear`, `burn_in` | 100, 80 years | 100, 95 years |
| Retained observations | 20 annual steps | 20 quarterly steps |
| `nage`, `rec.age` | 15 classes, 1 year | 20 classes, 0.25 year |
| `Linf`, `vbk` | 60, 0.2 per year | 152, 0.38 per year |
| `M`, `F_mean` | 0.2, 0.3 per annual step | 0.2, 0.2 per quarterly step |


These are the two package simulation presets, not recovered empirical datasets from the paper. Duration arguments are in **years**, whereas returned `sm$nyear` counts retained **steps**. Transpose `SN_at_len` (time × length); weight and maturity are already length × time.


`iter_range` provides replicate IDs and seeds; `return_iter` selects the in-memory replicate. `sim_rep*` files use `save()/load()`, not RDS. Independent error varies across cells, while shared-time error applies the same perturbation to all length bins in a period.

```r
# Generate both cases, CSV and truth files
cases <- run_simulation_examples(out = "simulation_examples", fit_batch = FALSE)
# Add two ACL fits and retain failures
# run_simulation_examples(out = "simulation_examples", fit_batch = TRUE)
e <- new.env()
load("simulation_examples/flatfish/sim_rep4", envir = e)
str(e$sim.data)
```


`sim_acl()` fits ACL only. For ALSCL, convert each replicate then call `run_alscl()` in a loop. Generator and estimator parameter names differ. Two replicates only verify the workflow; performance studies need more replicates, retained truth, bias, RMSE, coverage and failure rates.
