# 移除近期数据后，估计结果会怎样变化？

## 要解决的问题

新增观测后，资源评估可能修订历史生物量或补充量。回溯分析通过移除末端观测并重新拟合同一模型，检查估计的变化。它不衡量预测准确度。

## 重新拟合截短序列

以下完整设置使用模拟 YTF 示例及其条件生物假设。需要可用的编译器，并会执行多次拟合。

```r
library(ALSCL)
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
rb <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                            list(nyear = 3, train_times = 2, silent = TRUE)))
plot_retro(rb, facet_col = 2, rho_digits = 3)
rb$rho_text
```

## 检查轨迹与 Mohn rho

![仓库提供的 ALSCL 回溯示例](../assets/manual/figures/ytf/retro_ALSCL.png)

在每次截短序列的最后一个时间点，将其估计与完整拟合比较。正 rho 表示截短拟合在这些末端时期往往偏高，负值则相反。方向相反的偏差可能抵消，所以应查看单条轨迹及每次拟合的数值状态。参考估计为零时，相对差异无定义。

## 应用于自己的评估

将输入与生物参数替换为完整数据拟合时使用的同一套配置。固定参数、观测处理和时间单位应保持一致。`nyear` 计算移除的末端时间步：三个季度步并不是三年。每次截短后应保留足够观测，详见 [函数使用说明](retrospectives.html)。

## 决定下一步检查什么

持续同向的修订提示需要检查输入变化、模型假设和参数敏感性。rho 本身不能确定成因，也不能直接作为修正系数。作出结论前，应记录完整拟合、各次截短的诊断及末端时间。
