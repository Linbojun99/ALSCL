# 内置数据与首个拟合

## 拟合示例

完成 [安装](getting-started.html) 后，在同一个 R 会话中执行以下代码。两次拟合使用相同调查观测，但种群结构不同。代码创建 `a`（ACL）、`b`（ALSCL）、`x` 和 `inputs`，后续教程继续使用这些对象。

```r
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
# 检查维度与来源
str(inputs)
x$provenance
# 两模型使用相同观测及条件设置
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![模型比较](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)

|对象|内容和适用方式|
|---|---|
| `YTF` |原有 26 项参数列表，无 `data.CatL`|
| `YTF_example` |23 个体长组 × 20 期；三张表、参数、模拟真值、拟合设置、来源|
| `example_data` |从原 ALSCL 保留的 9 组 × 21 年三张表；来源待核实|

`YTF_example` 复用 `YTF` 的共有参数，特别是调查对数 SD = 0.1；其余采用 flatfish 模拟预设。随机种子 4，预热 80 年、保留 20 年；2000–2019 仅为示意年份。教学单位约定为 cm、kg，数据不是论文实测样本。调用 [data-raw/YTF_example.R](https://github.com/Linbojun99/ALSCL/blob/main/data-raw/YTF_example.R) 可重建 `.rda`，脚本恢复调用者的随机状态。


`x$fit_config` 固定已知生长、部分过程误差及相关参数，让示例聚焦工作流程。区间是在这些固定参数条件下的近似区间，不包含全部生物不确定性。这种教学拟合不能证明同时估计所有参数时可识别。

## 检查拟合结果

```r
names(b)
names(b$report)
diagnose_model(x$data.CatL, b)
plot_residuals(b, type = "year")
```

`report` 包含派生种群量与调查预测，`opt` 包含优化结果，`est_std` 和 `vcov` 描述近似不确定性，`obj` 是当前会话的 TMB 目标函数。模型比较和回溯分析应保留原始设置。解释估计值前，请先阅读 [诊断指南](diagnostics.html)。
