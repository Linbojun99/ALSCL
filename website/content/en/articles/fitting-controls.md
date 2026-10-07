# Fitting controls

Complete the [worked YTF example](first-model.html) first. The examples below use its `x`, `inputs`, `a` (ACL) and `b` (ALSCL) objects.

|Setting|Meaning and effect|
|---|---|
| `rec.age`, `nage` |Recruitment age and number of age classes|
| `growth_step` |Years per step. ACL NULL uses rec.age if below 1, otherwise 1; ALSCL defaults to 1. Set 0.25 explicitly for quarterly data.|
| `M` |Instantaneous natural mortality per model step|
| `sel_L50`, `sel_L95` |Fixed survey catchability locations|
| `len_mid`, `len_border` |Midpoints and internal boundaries; ALSCL also offers outer bounds|
| `parameters`, `.L`, `.U` |Named starting values and bounds on the optimizer scale|
| `map` |Fix or release parameters at their supplied values|
| `train_times` |Successive optimization passes, not independent random starts|
| `control` |Optimizer controls|
| `ncores` |Independent fit starts versus outer replicate/peel workers|
| `output` |Export diagnostics and plots under output when TRUE|


Simulation presets and estimation starts are separate. Choosing a species in `create_parameters()` does not automatically change the biological arguments passed to `run_*()`.

```r
p <- create_parameters(model_type = "acl", species = "flatfish")
names(p) # Inspect start and bound components
# Explicit fixed parameters, model-specific names
fixed <- list(log_Linf = factor(NA), log_vbk = factor(NA))
# NULL releases a default fixed item
released <- generate_map(list(logit_log_F_y = NULL))
VB_func(Linf = 60, k = 0.2, t0 = 1/60, age = 1:15)
mat_func(L50 = 35, L95 = 40, length = seq(6, 50, 2))
```

| Meaning | ACL parameter | ALSCL parameter |
|---|---|---|
|Growth| `log_Linf`, `log_vbk`, `t0` | `log_Linf`, `log_vbk`, `log_t0` |
|Length dispersion| `log_cv_len` | `log_cv_len`, `log_cv_grow` |
|Process SD| `log_std_log_R`, `log_std_log_F`, `log_std_log_N0` | `log_sigma_log_R`, `log_sigma_log_F`, `log_sigma_log_N0` |
|F correlation| `logit_log_F_a`, `logit_log_F_y` | `logit_log_F_l`, `logit_log_F_y` |
|Recruitment correlation| `logit_log_R` | `logit_log_R` |


Most `log_` parameters use log-transformed values; correlation transforms must match the template. ACL `t0` is untransformed, while ALSCL exponentiates `log_t0`, preventing negative t0 in this parameterization. Do not interchange their parameter lists.
