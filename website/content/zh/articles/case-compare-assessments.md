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

## 比较种群估计

```r
plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![模拟 YTF 示例的条件种群估计](../assets/manual/figures/ytf/compare_ts.png)

观察时间位置、趋势与尺度的差异。上图为仓库已有示例图；新拟合结果由安装版本与数值设置决定。结果相近不证明参数可识别，单个示例也不能判断模型的普遍优劣。

## 继续追查差异

审查固定的生物假设、残差结构以及对初值的敏感性。[模型比较绘图](model-comparison.html) 提供调查拟合、死亡率与生长的对照；[回溯分析案例](case-retrospective-bias.html) 检查移除末端观测后估计如何修订。
