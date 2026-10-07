# Estimation parameter lists

The keys below apply to starting values and bounds. Transformations map optimizer values to the natural scale; enter starts and bounds on the optimizer scale. The constructor's three named components match the corresponding fit arguments.

| ACL key | ALSCL key | Transform to natural scale | Meaning |
|---|---|---|---|
| `log_init_Z` | `log_init_Z` | exp |Initial equilibrium total mortality|
| `log_std_log_N0` | `log_sigma_log_N0` | exp |Initial log-number deviation SD|
| `mean_log_R` | `mean_log_R` | exp |Location of the log recruitment process|
| `log_std_log_R` | `log_sigma_log_R` | exp |Log recruitment process SD scale|
| `logit_log_R` | `logit_log_R` | plogis |Recruitment correlation restricted to (0,1)|
| `mean_log_F` | `mean_log_F` | exp |Location of the log F process|
| `log_std_log_F` | `log_sigma_log_F` | exp |F deviation SD scale|
| `logit_log_F_y` | `logit_log_F_y` | plogis |Temporal F correlation|
| `logit_log_F_a` | `logit_log_F_l` | plogis |Age/length F correlation|
| `log_vbk` | `log_vbk` | exp |Annual VB k|
| `log_Linf` | `log_Linf` | exp |Asymptotic length|
| `t0` | `log_t0` |Identity vs exp|ALSCL restricts t0 to positive values|
| `log_cv_len` | `log_cv_len` | exp |Length-at-age CV|
| — | `log_cv_grow` | exp |Growth-increment CV|
| `log_std_index` | `log_sigma_index` | exp |Log survey observation SD|


Exponentiating a log location does not give the unconditional arithmetic mean when random deviations are present. Both templates currently restrict AR(1) correlations to positive values. See the theory section for SD conventions in correlated processes.

```r
# Inspect all starts and bounds for a preset
p <- create_parameters(model_type = "alscl", species = "flatfish")
keys <- names(p$parameters)
parameter_table <- data.frame(
  Parameter = keys,
  Start = unlist(p$parameters[keys]),
  Lower = unlist(p$parameters.L[keys]),
  Upper = unlist(p$parameters.U[keys]))
parameter_table
# Change a start and upper bound
p$parameters$log_sigma_index <- log(0.15)
p$parameters.U$log_sigma_index <- log(0.5)
# Explicit map fixes the value
fixed_index <- list(log_sigma_index = factor(NA))
```


Fits construct random effects `dev_log_R`, `dev_log_F` and `dev_log_N0` from data dimensions. These are not ordinary scalar entries above. Generator names such as `std_logR` and `F_mean` are not estimator keys.
