# 调查与种群图

先完成 [YTF 拟合示例](first-model.html)，得到 `a`、`b`、`x`、`inputs`，并设置 `dat <- x$data.CatL`。加载 ggplot2 和 patchwork。[回溯分析指南](retrospectives.html) 建立 `ra`、`rb`。

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## 调查拟合 length 对数尺度

```r
# 绘制本图
plot_CatL(b, type = "length", exp_transform = FALSE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015),
  ylab = "Log survey index")
```

[![调查拟合 length 对数尺度](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_FALSE.png?raw=true)

点为观测对数，曲线为拟合的 Elog_index。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## 调查拟合 length 原尺度

```r
# 绘制本图
plot_CatL(b, type = "length", exp_transform = TRUE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015),
  ylab = "Survey index")
```

[![调查拟合 length 原尺度](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_TRUE.png?raw=true)

原尺度曲线为 exp(Elog_index)，即对数正态中位数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## 调查拟合 year 对数尺度

```r
# 绘制本图
plot_CatL(b, type = "year", exp_transform = FALSE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(10, 30, 50),
  ylab = "Log survey index")
```

[![调查拟合 year 对数尺度](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_FALSE.png?raw=true)

朱红色为观测对数，深蓝色为拟合的 Elog_index。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## 调查拟合 year 原尺度

```r
# 绘制本图
plot_CatL(b, type = "year", exp_transform = TRUE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(10, 30, 50),
  ylab = "Survey index")
```

[![调查拟合 year 原尺度](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_TRUE.png?raw=true)

原尺度曲线为 exp(Elog_index)，即对数正态中位数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## 丰度 N

```r
# 绘制本图
plot_abundance(b, type = "N", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![丰度 N](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_N.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_N.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## 丰度 NA

```r
# 绘制本图
plot_abundance(b, type = "NA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![丰度 NA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NA.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## 丰度 NL

```r
# 绘制本图
plot_abundance(b, type = "NL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![丰度 NL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NL.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## 生物量 B

```r
# 绘制本图
plot_biomass(b, type = "B", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![生物量 B](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_B.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_B.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## 生物量 BA

```r
# 绘制本图
plot_biomass(b, type = "BA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![生物量 BA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BA.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## 生物量 BL

```r
# 绘制本图
plot_biomass(b, type = "BL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![生物量 BL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BL.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## 产卵生物量 SSB

```r
# 绘制本图
plot_SSB(b, type = "SSB", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![产卵生物量 SSB](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SSB.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SSB.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## 产卵生物量 SBA

```r
# 绘制本图
plot_SSB(b, type = "SBA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![产卵生物量 SBA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBA.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## 产卵生物量 SBL

```r
# 绘制本图
plot_SSB(b, type = "SBL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![产卵生物量 SBL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBL.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 总量图显示各组之和；分面图的纵轴按组调整。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## 捕获尾数 CN

```r
# 绘制本图
plot_catch(b, type = "CN", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![捕获尾数 CN](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CN.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CN.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 捕获尾数是模型推算的渔业捕获，不是调查指数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## 捕获尾数 CNA

```r
# 绘制本图
plot_catch(b, type = "CNA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![捕获尾数 CNA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNA.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 捕获尾数是模型推算的渔业捕获，不是调查指数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## 捕获尾数 CNL

```r
# 绘制本图
plot_catch(b, type = "CNL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![捕获尾数 CNL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNL.png?raw=true)

ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。 捕获尾数是模型推算的渔业捕获，不是调查指数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## 补充量

```r
# 绘制本图
plot_recruitment(b, se = TRUE)
```

[![补充量](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/recruitment.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/recruitment.png?raw=true)

补充量对应最小年龄组；受设定的补充年龄和时间步长影响。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_recruitment)



## 亲体与补充关系

```r
# 绘制本图
plot_SSB_Rec(b, age_at_recruitment = 1)
```

[![亲体与补充关系](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/SSB_Rec.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/SSB_Rec.png?raw=true)

横轴与补充量错开 1 个观测步长。这里只画散点，不拟合资源补充函数。


[参数及默认值](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB_Rec)
