# 模型比较与真值

先完成 [YTF 拟合示例](first-model.html)，得到 `a`、`b`、`x`、`inputs`，并设置 `dat <- x$data.CatL`。加载 ggplot2 和 patchwork。[回溯分析指南](retrospectives.html) 建立 `ra`、`rb`。

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## 种群时间序列比较

```r
# 绘制本图
plot_compare_ts(a, b, se = TRUE, ncol = 2)
```

[![种群时间序列比较](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)

只比较共同输出与相同时间窗；自由纵轴不能直接比较不同量的振幅。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_ts)



## 死亡率矩阵比较

```r
# 绘制本图
plot_compare_F(a, b) + scale_x_continuous(breaks = c(2000, 2010, 2019))
```

[![死亡率矩阵比较](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F.png?raw=true)

ACL 纵轴是年龄，ALSCL 纵轴是长度；颜色可辅助观察，但两行不逐格对应。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_F)



## 残差联合诊断

```r
# 绘制本图
plot_compare_residuals(a, b, dat)
```

[![残差联合诊断](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_residuals.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_residuals.png?raw=true)

包括直方图、QQ 图、年度残差和长度组箱线图。年度误差棒是标准差，不是置信区间。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_residuals)



## 生长与外部参照接口

参考数据与排版辅助函数

```r
# 模拟参考生长点，并非论文实测
set.seed(42)
ref <- data.frame(Age = rep(1:15, each = 5))
ref$Length <- VB_func(60, .2, 1/60, ref$Age) + rnorm(nrow(ref), 0, 2)
growth_display <- function(p) {
  p[[1]] <- p[[1]] + labs(subtitle = paste(strsplit(
    p[[1]]$labels$subtitle, " | ", fixed=TRUE)[[1]], collapse="\n")) +
    theme(legend.text=element_text(size=7))
  p[[2]] <- p[[2]] + scale_x_continuous(breaks=c(1,5,10,15))
  p
}
```

```r
# 绘制本图
growth_display(plot_compare_growth(a, b, age_range = c(1, 15), ref_data = ref, ref_name = "Synthetic reference",
    nls_start = list(Linf = 60, k = 0.2, t0 = 1/60)))
```

[![生长与外部参照接口](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_growth.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_growth.png?raw=true)

灰色参照曲线来自模拟 Age 与 Length；不代表论文实测生长。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_growth)



## 两模型调查拟合

```r
# 绘制本图
plot_compare_CatL(a, b, dat, years = c(2000, 2005, 2010, 2015), ncol = 2)
```

[![两模型调查拟合](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_CatL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_CatL.png?raw=true)

显示 4 个指定年份；省略 years 将显示全部年份。点为调查数据，线为拟合中位数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_CatL)



## 拟合指标比较

```r
# 绘制本图
plot_compare_metrics(a, b, dat, ncol = 2) & scale_y_continuous(expand = expansion(mult = c(0,
    0.3)))
```

[![拟合指标比较](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_metrics.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_metrics.png?raw=true)

当前数据由年龄模型生成。图示差异不能证明任一模型普遍较优；IC 还需要相同数据与似然口径。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_metrics)



## 固定调查可捕性

```r
# 绘制本图
plot_compare_selectivity(a, b)
```

[![固定调查可捕性](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_selectivity.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_selectivity.png?raw=true)

当前函数展示输入 q 的连接线，不是估计的渔业选择性。两模型使用同一输入，曲线重合。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_selectivity)



## 总体 F 比较 apical

```r
# 绘制本图
plot_compare_annual_F(a, b, method = "apical") +
  labs(title = "Apical fishing mortality")
```

[![总体 F 比较 apical](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_apical.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_apical.png?raw=true)

apical 取最大值，mean 取组间算术平均，均非丰度加权。函数名称含 annual，但不会自动把季度率年化。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_annual_F)



## 总体 F 比较 mean

```r
# 绘制本图
plot_compare_annual_F(a, b, method = "mean") +
  labs(title = "Mean fishing mortality")
```

[![总体 F 比较 mean](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_mean.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_mean.png?raw=true)

apical 取最大值，mean 取组间算术平均，均非丰度加权。函数名称含 annual，但不会自动把季度率年化。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_annual_F)



## 模拟真值与条件估计

```r
# 将生成真值与报告中的估计对齐
years <- a$year
truth <- data.frame(Year=rep(years,3),
  Quantity=rep(c("B","SSB","Rec"),each=length(years)),
  Truth=c(x$truth$TB,x$truth$SSB,x$truth$Rec))
estimates <- rbind(
  data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(a$report$B,a$report$SSB,a$report$Rec),Model="ACL"),
  data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(b$report$B,b$report$SSB,b$report$Rec),Model="ALSCL"))
```

```r
# 绘制本图
ggplot(estimates, aes(Year, Value, color = Model)) + geom_line() + scale_color_manual(values = setNames(acl_theme("compare_colors"),
    c("ACL", "ALSCL"))) + geom_line(data = truth, aes(Year, Truth), inherit.aes = FALSE, linetype = 2) +
    facet_wrap(~Quantity, scales = "free_y", ncol = 1) + theme_bw() + labs(y = NULL)
```

[![模拟真值与条件估计](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/truth.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/truth.png?raw=true)

黑虚线为生成真值；固定生物参数，且模拟与拟合的补充过程不同。
