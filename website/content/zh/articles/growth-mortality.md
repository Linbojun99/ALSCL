# 生长、死亡与残差图

先完成 [YTF 拟合示例](first-model.html)，得到 `a`、`b`、`x`、`inputs`，并设置 `dat <- x$data.CatL`。加载 ggplot2 和 patchwork。[回溯分析指南](retrospectives.html) 建立 `ra`、`rb`。

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## 生长曲线与区间

```r
# 绘制本图
plot_VB(b, age_range = c(1, 15), se = TRUE)
```

[![生长曲线与区间](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/VB.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/VB.png?raw=true)

生长参数在本例中固定，区间退化为曲线；释放参数后区间才反映估计不确定性，不是个体长度分布。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_VB)



## 年龄长度转换

```r
# 绘制本图
plot_pla(b) + scale_x_discrete(labels = as.character(1:15))
```

[![年龄长度转换](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/pla.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/pla.png?raw=true)

色值是长度组给定年龄的概率；不是增长转移矩阵 G。年龄标签为组序号。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_pla)



## ACL 捕捞死亡率 year

```r
# 绘制本图
plot_fishing_mortality(a, type = "year",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL 捕捞死亡率 year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_year.png?raw=true)

年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ACL 捕捞死亡率 age

```r
# 绘制本图
plot_fishing_mortality(a, type = "age",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL 捕捞死亡率 age](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_age.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_age.png?raw=true)

年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ACL 捕捞死亡率 length

```r
# 绘制本图
plot_fishing_mortality(a, type = "length",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL 捕捞死亡率 length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_length.png?raw=true)

当前 ACL 的 length 分支与 age 分支相同，仍是年龄曲线，不能作为长度别 F。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL 捕捞死亡率 year

```r
# 绘制本图
plot_fishing_mortality(b, type = "year",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL 捕捞死亡率 year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_year.png?raw=true)

年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL 捕捞死亡率 age

```r
# 绘制本图
plot_fishing_mortality(b, type = "age",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL 捕捞死亡率 age](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_age.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_age.png?raw=true)

年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL 捕捞死亡率 length

```r
# 绘制本图
plot_fishing_mortality(b, type = "length",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL 捕捞死亡率 length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_length.png?raw=true)

年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## 对数残差 length

```r
# 绘制本图
plot_residuals(b, type = "length", facet_ncol = 4, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015))
```

[![对数残差 length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_length.png?raw=true)

观测对数减预测对数；未经标准差标准化。平滑线用于发现结构，不是显著性检验。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_residuals)



## 对数残差 year

```r
# 绘制本图
plot_residuals(b, type = "year", facet_ncol = 4, line_size = 0.6,
  x_breaks = c(10, 30, 50))
```

[![对数残差 year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_year.png?raw=true)

观测对数减预测对数；未经标准差标准化。平滑线用于发现结构，不是显著性检验。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_residuals)



## 过程偏差 R TRUE

```r
# 绘制本图
plot_deviance(b, type = "R", log = TRUE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![过程偏差 R TRUE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_TRUE.png?raw=true)

对数过程偏差以 0 为参照。这不是似然偏差统计量。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## 过程偏差 R FALSE

```r
# 绘制本图
plot_deviance(b, type = "R", log = FALSE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![过程偏差 R FALSE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_FALSE.png?raw=true)

指数转换后以 1 为参照，区间不对称。这不是原尺度残差。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## 过程偏差 F TRUE

```r
# 绘制本图
plot_deviance(b, type = "F", log = TRUE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![过程偏差 F TRUE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_TRUE.png?raw=true)

对数过程偏差以 0 为参照。这不是似然偏差统计量。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## 过程偏差 F FALSE

```r
# 绘制本图
plot_deviance(b, type = "F", log = FALSE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![过程偏差 F FALSE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_FALSE.png?raw=true)

指数转换后以 1 为参照，区间不对称。这不是原尺度残差。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## 长度组成山脊图

```r
# 绘制本图
plot_ridges(b) # 默认 viridis 渐变
```

[![长度组成山脊图](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/ridges.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/ridges.png?raw=true)

每期分别归一化；左为观测，右为拟合。不能由此读出总丰度趋势。需当前会话中的 obj。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_ridges)
