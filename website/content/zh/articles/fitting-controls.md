# 初值、边界与拟合

先完成 [YTF 拟合示例](first-model.html)。以下代码使用其中的 `x`、`inputs`、`a`（ACL）和 `b`（ALSCL）对象。

|设置|怎么调整及影响|
|---|---|
| `rec.age`, `nage` |起始年龄（年）和年龄组数；年龄网格为 rec.age + (0:nage-1) × growth_step|
| `growth_step` |每步多少年；ACL 的 NULL 按 rec.age<1 时取 rec.age、否则取 1，ALSCL 默认 1，季度应显式设 0.25|
| `M` |每模型步长自然死亡瞬时率；年率转季度需相应换算|
| `sel_L50`, `sel_L95` |固定调查 q 的位置和斜率，L95 必须大于 L50|
| `len_mid`, `len_border` |组中心及内部边界；ALSCL 还提供 `len_lower/len_upper`|
| `parameters`, `.L`, `.U` |命名初值、下界、上界；在程序参数化尺度上提供|
| `map` |固定或释放参数；固定值来自 `parameters`|
| `train_times` |从优化结果继续优化的次数，并非随机多起点|
| `control` |传给 `nlminb` 的控制项，如 `eval.max`、`iter.max`|
| `ncores` |单次拟合为独立初始点数，使用 socket 并行并扰动自由参数；批量/回溯为外层并行数|
| `output` |FALSE 不写文件，TRUE 在 output 目录导出诊断与图|

`initialize_params(species=...)` 建立模拟参数；`create_parameters(model_type=..., species=...)` 建立估计初值及边界。后者的物种预设不会自动替你更改 `run_*` 的 M、年龄、调查 q 或数据。


```r
p <- create_parameters(model_type = "acl", species = "flatfish")
names(p) # 查看初值及上下界列表
# 明确固定参数；不同模型参数名不同
fixed <- list(log_Linf = factor(NA), log_vbk = factor(NA))
# factor(NA) 固定；映射中的 NULL 可释放默认固定项
released <- generate_map(list(logit_log_F_y = NULL))
VB_func(Linf = 60, k = 0.2, t0 = 1/60, age = 1:15)
mat_func(L50 = 35, L95 = 40, length = seq(6, 50, 2))
```

|生物含义|ACL 参数名|ALSCL 参数名|
|---|---|---|
|生长| `log_Linf`, `log_vbk`, `t0` | `log_Linf`, `log_vbk`, `log_t0` |
|长度离散| `log_cv_len` | `log_cv_len`, `log_cv_grow` |
|过程 SD| `log_std_log_R`, `log_std_log_F`, `log_std_log_N0` | `log_sigma_log_R`, `log_sigma_log_F`, `log_sigma_log_N0` |
|F 相关| `logit_log_F_a`, `logit_log_F_y` | `logit_log_F_l`, `logit_log_F_y` |
|补充相关| `logit_log_R` | `logit_log_R` |

`log_` 参数通常用 `log(自然尺度值)`，相关系数的变换必须以对应模板为准。尤其 ACL 的 `t0` 是原尺度，而 ALSCL 的 `log_t0` 经 exp 变换，因此本版 ALSCL 的该参数化不支持负 t0。不要将 ACL 参数表直接用于 ALSCL。全部 14/15 个估计参数的名称、变换及上下界查询方法见 [初值列表参考](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#estimation_parameters)。
