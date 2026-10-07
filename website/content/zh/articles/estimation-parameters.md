# 估计参数列表

`parameters`、`parameters.L`、`parameters.U` 的键名如下。表中变换指从优化参数到自然尺度；初值和界限都要在优化尺度填写。`create_parameters()` 返回的三个命名列表就是传给拟合入口的同名参数。


|ACL 键名|ALSCL 键名|变换|含义|
|---|---|---|---|
| `log_init_Z` | `log_init_Z` | exp |初始平衡总死亡率|
| `log_std_log_N0` | `log_sigma_log_N0` | exp |初始年龄尾数偏差 SD|
| `mean_log_R` | `mean_log_R` | exp |对数补充过程的位置参数|
| `log_std_log_R` | `log_sigma_log_R` | exp |对数补充过程 SD 缩放|
| `logit_log_R` | `logit_log_R` | plogis |补充 AR(1) 相关，限定 (0,1)|
| `mean_log_F` | `mean_log_F` | exp |对数 F 过程的位置参数|
| `log_std_log_F` | `log_sigma_log_F` | exp |F 偏差 SD 缩放|
| `logit_log_F_y` | `logit_log_F_y` | plogis |F 的时间相关|
| `logit_log_F_a` | `logit_log_F_l` | plogis |F 的年龄/体长相关|
| `log_vbk` | `log_vbk` | exp |VB k，单位为每年|
| `log_Linf` | `log_Linf` | exp |渐近体长|
| `t0` | `log_t0` |ACL 原值；ALSCL exp|VB t0；ALSCL 本版限定正值|
| `log_cv_len` | `log_cv_len` | exp |年龄体长 CV|
| — | `log_cv_grow` | exp |增长增量 CV|
| `log_std_index` | `log_sigma_index` | exp |对数调查观测 SD|

`exp(mean_log_R)` 和 `exp(mean_log_F)` 是对数位置的指数转换；有随机偏差时不能直接称作无条件算术平均数。两个模板当前都不支持负的 AR(1) 相关。SD 在相关过程中的具体尺度见指南原理部分。


```r
# 查看某物种的全部初值和上下界
p <- create_parameters(model_type = "alscl", species = "flatfish")
keys <- names(p$parameters)
parameter_table <- data.frame(
  Parameter = keys,
  Start = unlist(p$parameters[keys]),
  Lower = unlist(p$parameters.L[keys]),
  Upper = unlist(p$parameters.U[keys]))
parameter_table
# 改一个初值及上界：此处并未自动固定该参数
p$parameters$log_sigma_index <- log(0.15)
p$parameters.U$log_sigma_index <- log(0.5)
# 要固定该 SD，还需显式 map；固定值取自 parameters
fixed_index <- list(log_sigma_index = factor(NA))
```

拟合会按数据维度建立 `dev_log_R`、`dev_log_F`、`dev_log_N0` 随机效应。它们不是本表的普通标量初值；不要把生成器的 `std_logR`、`F_mean` 等直接作为估计列表键名。
