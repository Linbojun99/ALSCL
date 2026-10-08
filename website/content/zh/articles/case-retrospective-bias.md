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

## 对照 ACL 的回溯结果

使用相同观测与截短次数，同时保留 ACL 自己的参数名与 map。

```r
ra <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl,
                          list(nyear = 3, train_times = 2, silent = TRUE)))
plot_retro(ra, facet_col = 2, rho_digits = 3)
ra$rho_text
```

![模拟 YTF：ACL 回溯轨迹](../assets/manual/figures/ytf/retro_ACL.png)

两种模型的配图均为仓库已有示例。先检查各模型内部的修订，再比较同一指标与末端时期的 rho。不同纵轴尺度可能放大或掩盖视觉差异。

## 检查返回表

```r
head(rb$results)
rb$last_points
rb$rho_text
names(rb)
```

`results` 保存轨迹，`last_points` 标出末端估计，`rho_text` 保存各指标的平均相对修订。包装函数 **不会** 返回每次截短的底层拟合对象或完整数值诊断表。

### 可选：逐次审查 ALSCL 截短拟合

以下循环额外拟合完整序列及三次截短，以保存数值状态和失败记录。它会重复上文的拟合工作，只在需要诊断审查时运行。保留标签列，并从三张表中移除相同的末端时间列。

```r
peel_checks <- lapply(0:3, function(k) {
  keep <- seq_len(ncol(inputs$data.CatL) - k)
  peeled <- lapply(inputs, function(d) d[, keep, drop = FALSE])
  tryCatch({
    fit <- do.call(run_alscl, c(peeled, x$fit_args, x$fit_config$alscl,
                               list(train_times = 2, silent = TRUE)))
    data.frame(Peel = k, Terminal = max(fit$year),
      Code = fit$convergence_code, Gradient = fit$max_abs_gradient,
      pdHess = fit$pdHess, Boundary = fit$bound_hit, Error = "")
  }, error = function(e) data.frame(Peel = k,
    Terminal = as.numeric(tail(names(peeled$data.CatL), 1)),
    Code = NA_integer_, Gradient = NA_real_, pdHess = NA,
    Boundary = NA, Error = conditionMessage(e)))
})
peel_checks <- do.call(rbind, peel_checks)
peel_checks
```

不能因为返回了曲线，就把不可靠拟合直接纳入科学解释。应报告失败并重新审查拟合假设。不存在能保证评估有效的通用梯度或 rho 阈值。

## 应用于自己的评估

将输入与生物参数替换为完整数据拟合时使用的同一套配置。固定参数、观测处理和时间单位应保持一致。`nyear` 计算移除的末端时间步：三个季度步并不是三年。每次截短后应保留足够观测，详见 [函数使用说明](retrospectives.html)。

## 导出回溯记录

```r
dir.create("my_retrospective", showWarnings = FALSE)
write.csv(rb$results, "my_retrospective/alscl_trajectories.csv", row.names = FALSE)
write.csv(rb$rho_text, "my_retrospective/alscl_rho.csv", row.names = FALSE)
write.csv(ra$rho_text, "my_retrospective/acl_rho.csv", row.names = FALSE)
p_retro <- plot_retro(rb, facet_col = 2, rho_digits = 3)
ggplot2::ggsave("my_retrospective/alscl.png", p_retro,
               width = 10, height = 7, dpi = 300)
saveRDS(list(ACL = ra, ALSCL = rb, fit_args = x$fit_args,
             fit_config = x$fit_config), "my_retrospective/results.rds")
writeLines(capture.output(sessionInfo()), "my_retrospective/session.txt")
```

|现象|解释或检查|
|---|---|
|截短后时间步太少|减少截短次数，每次拟合至少保留两个时间步|
|末端偏差正负交替|查看各次截短，平均 rho 可能掩盖抵消|
|将季度分析理解为年度分析|用 `growth_step` 将移除步数转换为时长|
|需要自定义 `zero_action`|回溯包装函数不提供该参数；使用拟合函数的显式循环，保持观测处理一致|
|短序列拟合不稳定|先检查诊断，再解释末端修订|

## 决定下一步检查什么

持续同向的修订提示需要检查输入变化、模型假设和参数敏感性。rho 本身不能确定成因，也不能直接作为修正系数。作出结论前，应记录完整拟合、各次截短的诊断及末端时间。
