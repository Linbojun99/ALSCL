# ACL 与 ALSCL 的评估结果有什么不同？

## 要解决的问题

你希望判断改变种群结构是否会影响评估。对同一组模拟 YTF 观测拟合 ACL 与 ALSCL，并保留示例提供的生物设置。这是条件比较，不能证明某个模型在所有情形下更好。

## 拟合同一组观测

内置数据为模拟数据。示例配置固定了已知生长及部分过程参数，因此拟合区间不包含这些假设的不确定性。

```r
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))
```

## 检查拟合是否可解释

```r
diagnose_model(x$data.CatL, a)
diagnose_model(x$data.CatL, b)
comparison <- compare_models(a, b, x$data.CatL)
comparison$summary
plot_compare_residuals(a, b, x$data.CatL)
```

比较估计轨迹前，先检查数值状态与残差结构。只有观测数据及似然定义可比较时，才比较基于似然的准则。不能仅凭目标函数值对数值不可靠的拟合排序。

![模拟 YTF：ACL 与 ALSCL 的残差比较](../assets/manual/figures/ytf/compare_residuals.png)

各面板对比残差分布以及随时间、体长的变化。年度误差线表示离散程度，不是置信区间。模型即使拟合出相似生物量趋势，也可能保留系统性观测残差。

## 比较调查观测的拟合

为两种模型选择相同年份，检查约束各次拟合的观测。

```r
plot_compare_CatL(a, b, x$data.CatL,
                  years = c(2000, 2005, 2010, 2015), ncol = 2)
```

![模拟 YTF：四个年份的调查观测及拟合中位数](../assets/manual/figures/ytf/compare_CatL.png)

点为观测，曲线为拟合中位数。讨论派生丰度差异前，先比较拟合形状，特别是观测稀疏的大体长组。上图使用仓库中的模拟拟合。

## 比较种群估计

```r
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![模拟 YTF 示例的条件种群估计](../assets/manual/figures/ytf/compare_ts.png)

观察时间位置、趋势与尺度的差异。上图为仓库已有示例图；新拟合结果由安装版本与数值设置决定。结果相近不证明参数可识别，单个示例也不能判断模型的普遍优劣。

## 检查捕捞死亡率的尺度

ACL 按年龄表示捕捞死亡率，ALSCL 按体长表示。因此面板坐标不同，不能逐单元格一一比较。

```r
plot_compare_F(a, b)
```

![模拟 YTF：年龄别与体长别捕捞死亡率](../assets/manual/figures/ytf/compare_F.png)

比较时间变化时，应考虑不同的状态变量。若使用 `plot_compare_annual_F()`，不要因函数名含 annual 就认为季度值会自动年化。

## 保存可复核的比较结果

```r
dir.create("my_model_comparison", showWarnings = FALSE)
write.csv(comparison$summary, "my_model_comparison/summary.csv", row.names = FALSE)
write.csv(comparison$fit_metrics, "my_model_comparison/fit_metrics.csv", row.names = FALSE)
p_population <- plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
ggplot2::ggsave("my_model_comparison/population.png", p_population,
               width = 10, height = 7, dpi = 300)
saveRDS(list(inputs = inputs, fit_args = x$fit_args, fit_config = x$fit_config,
             comparison = comparison, ACL_report = a$report, ALSCL_report = b$report),
        "my_model_comparison/summary.rds")
writeLines(capture.output(sessionInfo()), "my_model_comparison/session.txt")
```

|观察|能支持什么|不能证明什么|
|---|---|---|
|两次数值拟合均可接受|值得进一步比较拟合结果|全部假设正确|
|残差结构不同|模型解释观测分布的方式不同|某个模型普遍更好|
|生物量或补充量不同|评估结论依赖种群结构或相关假设|哪个是实测种群的真值|
|区间重叠|条件估计未被这些区间清楚分开|等价或管理后果相同|

处理自己的观测时，应同时替换两种模型的输入，并审查各模型的参数名。不能分别拟合不同数据，却把全部差异归因于模型结构。报告时间窗、固定生物参数、零值处理、数值诊断与敏感性结果。

## 继续追查差异

审查固定的生物假设、残差结构以及对初值的敏感性。[模型比较绘图](model-comparison.html) 提供调查拟合、死亡率与生长的对照；[回溯分析案例](case-retrospective-bias.html) 检查移除末端观测后估计如何修订。
