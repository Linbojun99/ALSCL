# 函数与参数参考 · Function reference

对应 ALSCL 2.0.0 的 **43 个公开函数**。签名及默认值从实际 R 函数提取；每项参数均有中文说明和英文帮助。签名中的 `c(...)` 通常列出可选值，默认取第一项。`NULL` 的含义依函数而异，未必表示停用。

All **43 exported functions** in ALSCL 2.0.0 are covered. Signatures/defaults were extracted from R. Each argument has Chinese guidance and English help. For a choice argument, `c(...)` lists options and the first is the default.

[使用指南 / Guide](USER_GUIDE.md) · [图谱 / Gallery](PLOT_GALLERY.md)

## 示例环境 / Example context

先执行使用指南中的内置数据拟合，得到 `x`、`inputs`、`a`（ACL）、`b`（ALSCL），并设 `dat <- x$data.CatL`。回溯示例生成 `ra`、`rb`；模拟依次运行 `initialize_params`、`sim_cal`、`sim_data`。独立参数函数无需拟合对象。

First run the guide's built-in fitting example to define `x`, `inputs`, `a` and `b`; set `dat <- x$data.CatL`. Retrospective examples define `ra/rb`. Run simulation constructors in order. Standalone biological functions need no fitted model.

```r
library(ALSCL)
library(ggplot2)
library(patchwork)
data(YTF_example)
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
dat <- x$data.CatL
```

## 索引 / Index

- [acl_theme](#acl_theme) — 读取全局主题 / Read global theme settings
- [acl_theme_reset](#acl_theme_reset) — 恢复默认主题 / Reset theme settings
- [acl_theme_set](#acl_theme_set) — 设置全局主题 / Set global theme settings
- [compare_models](#compare_models) — 生成两个模型的数值比较表 / Compare fitted model summaries
- [create_parameters](#create_parameters) — 建立模型初值及边界 / Construct estimation starts and bounds
- [create_parameters_alscl](#create_parameters_alscl) — ALSCL 初值及边界便捷入口 / ALSCL parameter constructor
- [diagnose_model](#diagnose_model) — 汇总诊断指标 / Summarize diagnostics
- [diagnostic_metrics](#diagnostic_metrics) — 计算预测误差和信息指标 / Calculate prediction and fit metrics
- [generate_map](#generate_map) — 合并 ACL 默认固定参数与自定义映射 / Merge ACL default and custom maps
- [initialize_params](#initialize_params) — 建立模拟生物参数 / Initialize simulation biology
- [mat_func](#mat_func) — 计算 logistic 成熟比例 / Calculate logistic maturity
- [plot_abundance](#plot_abundance) — 丰度 N / Abundance N
- [plot_biomass](#plot_biomass) — 生物量 B / Biomass B
- [plot_catch](#plot_catch) — 捕获尾数 CN / Catch numbers CN
- [plot_CatL](#plot_CatL) — 调查拟合 length 对数尺度 / Survey fit length log scale
- [plot_compare_annual_F](#plot_compare_annual_F) — 总体 F 比较 apical / Summary F comparison apical
- [plot_compare_CatL](#plot_compare_CatL) — 两模型调查拟合 / Survey fits from both models
- [plot_compare_F](#plot_compare_F) — 死亡率矩阵比较 / Fishing mortality matrices
- [plot_compare_growth](#plot_compare_growth) — 生长与外部参照接口 / Growth and reference-data interface
- [plot_compare_metrics](#plot_compare_metrics) — 拟合指标比较 / Fit metric comparison
- [plot_compare_residuals](#plot_compare_residuals) — 残差联合诊断 / Combined residual diagnostics
- [plot_compare_selectivity](#plot_compare_selectivity) — 固定调查可捕性 / Fixed survey catchability
- [plot_compare_ts](#plot_compare_ts) — 种群时间序列比较 / Population time series comparison
- [plot_deviance](#plot_deviance) — 过程偏差 R TRUE / Process deviations R TRUE
- [plot_fishing_mortality](#plot_fishing_mortality) — ACL 捕捞死亡率 year / ACL fishing mortality year
- [plot_pla](#plot_pla) — 年龄长度转换 / Age length conversion
- [plot_recruitment](#plot_recruitment) — 补充量 / Recruitment
- [plot_residuals](#plot_residuals) — 对数残差 length / Log residuals length
- [plot_retro](#plot_retro) — ACL 回溯分析 / ACL retrospective analysis
- [plot_ridges](#plot_ridges) — 长度组成山脊图 / Length composition ridges
- [plot_SSB](#plot_SSB) — 产卵生物量 SSB / Spawning biomass SSB
- [plot_SSB_Rec](#plot_SSB_Rec) — 亲体与补充关系 / Spawning biomass and recruitment
- [plot_VB](#plot_VB) — 生长曲线与区间 / Growth curve and interval
- [retro_acl](#retro_acl) — ACL 回溯入口 / ACL retrospective wrapper
- [retro_alscl](#retro_alscl) — ALSCL 回溯入口 / ALSCL retrospective wrapper
- [retro_model](#retro_model) — 统一回溯分析 / Unified retrospective analysis
- [run_acl](#run_acl) — 拟合年龄型 ACL / Fit age-based ACL
- [run_alscl](#run_alscl) — 拟合年龄与长度联合 ALSCL / Fit joint age-length ALSCL
- [sim_acl](#sim_acl) — 批量拟合 ACL 模拟重复 / Fit simulated replicates with ACL
- [sim_cal](#sim_cal) — 从参数计算生物矩阵 / Calculate biological arrays
- [sim_data](#sim_data) — 运行完整随机种群模拟 / Simulate full population dynamics
- [simulate_example_data](#simulate_example_data) — 快速年度教学数据生成 / Generate simple annual teaching data
- [VB_func](#VB_func) — 计算 VB 平均体长 / Calculate VB mean length

<a id="acl_theme"></a>
## `acl_theme()`

读取全局主题 / Read global theme settings

```r
acl_theme(what = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `what` | `NULL` | NULL 返回全部全局主题设置，字符名称返回单项。 | Optional character. Retrieve a specific element, e.g. "titles", "font_family", "compare_colors". |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list of current settings, or a single element if what is specified.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
acl_theme()
acl_theme("line_color")
```

<a id="acl_theme_reset"></a>
## `acl_theme_reset()`

恢复默认主题 / Reset theme settings

```r
acl_theme_reset()
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
acl_theme_reset()
```

<a id="acl_theme_set"></a>
## `acl_theme_set()`

设置全局主题 / Set global theme settings

```r
acl_theme_set(base_theme = NULL, font_family = NULL, title_size = NULL,
    title_hjust = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    x_expand = NULL, line_color = NULL, line_size = NULL, se_color = NULL,
    se_alpha = NULL, compare_colors = NULL, compare_linetypes = NULL,
    compare_linewidth = NULL, compare_point_size = NULL, compare_se_alpha = NULL,
    compare_legend_pos = NULL, compare_facet_ncol = NULL, compare_facet_scales = NULL,
    titles = NULL, xlab = NULL, ylab = NULL, palette = NULL,
    point_color = NULL, observed_color = NULL, smooth_color = NULL,
    hline_color = NULL, low_col = NULL, high_col = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character. Base ggplot2 theme name. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character. Font family for all text. |
| `title_size` | `NULL` | 图标题字号。 | Numeric. Plot title size in pt. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric. Title alignment (0=left, 0.5=center, 1=right). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric. Axis title size. |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric. Axis tick label size. |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric. Facet label size. |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric. Legend text size. |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks. |
| `x_expand` | `NULL` | 横轴留白设置，长度 2 的数值向量。 | Numeric vector of length 2. X-axis expansion. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character. Default line color for single-model plots. |
| `line_size` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. Default line width for single-model plots. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character. Default CI color. |
| `se_alpha` | `NULL` | 区间透明度，0–1。 | Numeric. Default CI transparency. |
| `compare_colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Character vector of 2 colors. Can be unnamed (positional) or named. The first color is always model1, second is model2. Examples: c("blue", "red"), c(ACL = "blue", ALSCL = "red"), c("M=0" = "steelblue", "M=0.5" = "tomato"). |
| `compare_linetypes` | `NULL` | 按模型顺序提供两种线型。 | Character vector of 2 linetypes. Same rules as colors. |
| `compare_linewidth` | `NULL` | 比较图全局线宽。 | Numeric. Line width for comparison plots. |
| `compare_point_size` | `NULL` | 比较图全局点大小。 | Numeric. Point size for comparison plots. |
| `compare_se_alpha` | `NULL` | 比较图区间透明度。 | Numeric. CI ribbon alpha for comparison plots. |
| `compare_legend_pos` | `NULL` | 图例位置，如 bottom、right、none。 | Character. Legend position for comparison plots. |
| `compare_facet_ncol` | `NULL` | 比较图全局分面列数。 | Integer. Default facet columns in comparison plots. |
| `compare_facet_scales` | `NULL` | 比较图全局分面尺度。 | Character. Facet scales for comparison plots. |
| `titles` | `NULL` | 以图名为键的标题列表，用于全局覆盖。 | Named list. Override specific titles. |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Named list. Override specific x-axis labels. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Named list. Override specific y-axis labels. |
| `palette` | `NULL` | ggsci 色板名；NULL 继承全局设置，部分图也保留旧色板。 | Character or NULL. ggsci palette: "npg" (default), "aaas", "nejm", "lancet", or "jama". Changing it resets palette-derived colors; explicit color arguments in the same call take precedence. No arguments leaves settings unchanged; use acl_theme_reset() to reset. |
| `point_color` | `NULL` | 点颜色。 | Character or NULL. Default point, observation, smoother and reference-line colors. |
| `observed_color` | `NULL` | 观测数据颜色；NULL 继承全局角色。 | Character or NULL. Default point, observation, smoother and reference-line colors. |
| `smooth_color` | `NULL` | 残差平滑线颜色。 | Character or NULL. Default point, observation, smoother and reference-line colors. |
| `hline_color` | `NULL` | 残差零参考线颜色。 | Character or NULL. Default point, observation, smoother and reference-line colors. |
| `low_col` | `NULL` | 概率热图的低值 / 高值颜色。 | Character or NULL. Continuous heatmap endpoints. |
| `high_col` | `NULL` | 概率热图的低值 / 高值颜色。 | Character or NULL. Continuous heatmap endpoints. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 Invisible previous settings.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
acl_theme_set(base_theme = "bw", line_color = "navy")
```

<a id="compare_models"></a>
## `compare_models()`

生成两个模型的数值比较表 / Compare fitted model summaries

```r
compare_models(model1, model2, data.CatL, model1_name = NULL, model2_name = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | A list from run_acl or run_alscl. Treated as the reference model. |
| `model2` | `必填 / required` | 第二个拟合对象。 | A list from run_acl or run_alscl. Treated as the alternative model. |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | A data frame with survey catch-at-length data (used for goodness-of-fit). |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Label for model1. Default is auto-detected from model_type. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Label for model2. Default is auto-detected from model_type. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with: summaryData frame with convergence, objective, AIC, BIC, number of parameters. growthData frame comparing VB growth parameters. fit_metricsData frame with MSE, RMSE, R-squared, MAPE, etc. for each model. correlationData frame with Pearson correlations and Ratio_Mean, the ratio of means (model2/model1), for common population quantities. model1_name,model2_nameLabels used in the comparison.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
compare_models(a, b, dat)
```

<a id="create_parameters"></a>
## `create_parameters()`

建立模型初值及边界 / Construct estimation starts and bounds

```r
create_parameters(model_type = c("acl", "alscl"), species = NULL, parameters = NULL,
    parameters.L = NULL, parameters.U = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_type` | `c("acl", "alscl")` | 模型选择：acl 或 alscl。 | Character. Either "acl" or "alscl". Determines which TMB model parameter names and defaults to use. |
| `species` | `NULL` | 物种预设 flatfish、tuna、krill；? 列出预设。 | Character or NULL. Optional species preset: "flatfish" (M7 yellowtail flounder-like, annual, Linf=60), "tuna" (M6 tuna-like, quarterly, Linf=152), or "krill" (Antarctic krill-like, Linf=60-90). When NULL (default), generic wide-range defaults are used. Each preset provides species-specific initial values and bounds that have been tested in simulation studies. See Details. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | A list of custom initial values (default NULL = use all defaults). Only the parameters you specify will be overridden; all others keep their default (or species-preset) values. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | A list of custom lower bounds (default NULL). |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | A list of custom upper bounds (default NULL). |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with three elements: parameters, parameters.L, and parameters.U. Each is a named list of parameter values with names matching the corresponding TMB model template.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
create_parameters(model_type = "acl", species = "flatfish")
```

<a id="create_parameters_alscl"></a>
## `create_parameters_alscl()`

ALSCL 初值及边界便捷入口 / ALSCL parameter constructor

```r
create_parameters_alscl(species = NULL, parameters = NULL, parameters.L = NULL,
    parameters.U = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `species` | `NULL` | 物种预设 flatfish、tuna、krill；? 列出预设。 | Character or NULL. Optional species preset: "flatfish" (M7 yellowtail flounder-like, annual, Linf=60), "tuna" (M6 tuna-like, quarterly, Linf=152), or "krill" (Antarctic krill-like, Linf=60-90). When NULL (default), generic wide-range defaults are used. Each preset provides species-specific initial values and bounds that have been tested in simulation studies. See Details. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | A list of custom initial values (default NULL = use all defaults). Only the parameters you specify will be overridden; all others keep their default (or species-preset) values. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | A list of custom lower bounds (default NULL). |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | A list of custom upper bounds (default NULL). |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with three elements: parameters, parameters.L, and parameters.U. Each is a named list of parameter values with names matching the corresponding TMB model template.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
create_parameters_alscl(species = "flatfish")
```

<a id="diagnose_model"></a>
## `diagnose_model()`

汇总诊断指标 / Summarize diagnostics

```r
diagnose_model(data.CatL, model_result)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Length labels followed by time columns. |
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A fitted model result. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A data frame of fit metrics and convergence diagnostics.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
diagnose_model(dat, b)
```

<a id="diagnostic_metrics"></a>
## `diagnostic_metrics()`

计算预测误差和信息指标 / Calculate prediction and fit metrics

```r
diagnostic_metrics(data.CatL, model_result)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Length labels followed by time columns. |
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A fitted ACL or ALSCL result. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A data frame with Metric and Value. Undefined metrics are NA.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
diagnostic_metrics(dat, b)
```

<a id="generate_map"></a>
## `generate_map()`

合并 ACL 默认固定参数与自定义映射 / Merge ACL default and custom maps

```r
generate_map(map = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | A list containing the custom values for the map elements (default is NULL). log_std_log_F: Custom value for log_std_log_F (default is NA). logit_log_F_y: Custom value for logit_log_F_y (default is NA). logit_log_F_a: Custom value for logit_log_F_a (default is NA). t0: Custom value for t0 (default is NA). Growth parameters log_vbk and log_Linf are free unless explicitly mapped. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list representing the generated map with custom or default values.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
generate_map(list(log_Linf = factor(NA), logit_log_F_y = NULL))
```

该函数提供 ACL 默认映射。ALSCL 请使用本例 `x$fit_config$alscl$map` 或模型相应命名参数。默认 `t0/log_std_log_F/logit_log_F_y/logit_log_F_a` 固定，`log_vbk/log_Linf` 默认可估计；自定义未列出的项仍保留默认映射。 / This constructor uses ACL defaults. Use model-specific names for ALSCL. Unmentioned default mappings remain in force.

<a id="initialize_params"></a>
## `initialize_params()`

建立模拟生物参数 / Initialize simulation biology

```r
initialize_params(nyear = NULL, rec.age = NULL, first.year = NULL, nage = NULL,
    M = NULL, init_Z = NULL, vbk = NULL, Linf = NULL, t0 = NULL,
    len_mid = NULL, len_border = NULL, a = NULL, b = NULL, mat_L50 = NULL,
    mat_L95 = NULL, cv_L = NULL, cv_inc = NULL, std_logR = NULL,
    std_logN0 = NULL, alpha = NULL, beta = NULL, std_SN = NULL,
    q_surv_L50 = NULL, q_surv_L95 = NULL, R_init = NULL, F_mean = NULL,
    F_ar = NULL, F_sd = NULL, R_ar = NULL, species = NULL, growth_step = NULL,
    burn_in = NULL, observation_error = c("independent", "shared_time"))
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `nyear` | `NULL` | 模拟器为总年数；回溯函数为删除的末端时间步数。 | Integer. Total simulation years including burn-in (default 100). |
| `rec.age` | `NULL` | 补充年龄，单位为年。 | Numeric. Age of recruitment (default 1; use 0.25 for quarterly). |
| `first.year` | `NULL` | 模拟返回观测的起始年份标签。 | Integer. First year of the observation window (default 2020). |
| `nage` | `NULL` | 年龄组数量，包括最大年龄加组。 | Integer. Number of age classes including plus-group (default 15). |
| `M` | `NULL` | 每个模型时间步的自然死亡瞬时率。 | Numeric. Natural mortality rate (default 0.2). |
| `init_Z` | `NULL` | 构造初始平衡年龄结构的总死亡率。 | Numeric. Initial total mortality for equilibrium (default 0.5). |
| `vbk` | `NULL` | VB 生长系数，以年为时间单位。 | Numeric. Von Bertalanffy k (default 0.2). |
| `Linf` | `NULL` | VB 渐近体长，与输入体长单位一致。 | Numeric. Asymptotic length (default 60). |
| `t0` | `NULL` | VB 理论零体长年龄，单位为年。 | Numeric. Von Bertalanffy t0 (default 1/60). |
| `len_mid` | `NULL` | 每个体长组的中心值向量，长度与数据行数相同。 | Numeric vector. Midpoints of length bins. |
| `len_border` | `NULL` | 体长边界；拟合入口需要 nlen-1 个内部边界，模拟器使用 nlen+1 个完整边界。 | Complete boundary vector of nlen+1 entries, including outer boundaries. |
| `a` | `NULL` | 长度重量关系 W=a*L^b 的系数，需匹配单位。 | Numeric. Length-weight coefficient (default exp(-12)). |
| `b` | `NULL` | 长度重量关系 W=a*L^b 的指数。 | Numeric. Length-weight exponent (default 3). |
| `mat_L50` | `NULL` | 成熟比例为 0.5 的体长。 | Numeric. Length at 50 pct maturity (default 35). |
| `mat_L95` | `NULL` | 成熟比例为 0.95 的体长，须大于 L50。 | Numeric. Length at 95 pct maturity (default 40). |
| `cv_L` | `NULL` | 给定年龄时体长的变异系数。 | Numeric. CV of length-at-age (default 0.2). |
| `cv_inc` | `NULL` | 生长增量的变异系数。 | Numeric. CV of growth increment (default 0.2). |
| `std_logR` | `NULL` | 模拟对数补充偏差的 SD 缩放。 | Numeric. SD of log-recruitment deviations (default 0.3). |
| `std_logN0` | `NULL` | 初始对数年龄组尾数偏差的 SD。 | Numeric. SD of initial log-N deviations (default 0.2). |
| `alpha` | `NULL` | Beverton–Holt 补充关系系数，量纲随丰度与生物量单位改变。 | Numeric. Beverton-Holt alpha (default 400). |
| `beta` | `NULL` | Beverton–Holt 补充关系系数，量纲随丰度与生物量单位改变。 | Numeric. Beverton-Holt beta (default 10). |
| `std_SN` | `NULL` | 调查观测对数误差的标准差。 | Numeric. Survey observation error SD (default 0.2). |
| `q_surv_L50` | `NULL` | 调查可捕性达到 50% 的体长。 | Numeric. Survey selectivity L50 (default 15). |
| `q_surv_L95` | `NULL` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Numeric. Survey selectivity L95 (default 20). |
| `R_init` | `NULL` | 初始补充尾数。 | Numeric. Initial recruitment (default 500). |
| `F_mean` | `NULL` | 每步平均捕捞死亡率。 | Numeric. Mean fishing mortality (default 0.3). |
| `F_ar` | `NULL` | 模拟 F 偏差的 AR(1) 相关参数。 | Numeric. AR1 coefficient for F deviations (default 0.75). |
| `F_sd` | `NULL` | 模拟 F 偏差标准差。 | Numeric. SD of F deviations (default 0.2). |
| `R_ar` | `NULL` | 模拟补充偏差的 AR(1) 相关参数。 | Numeric. AR1 coefficient for recruitment deviations (default 0.1). |
| `species` | `NULL` | 物种预设 flatfish、tuna、krill；? 列出预设。 | Character or NULL. Species preset name, one of "flatfish", "tuna", "krill", or "?" to list presets. NULL (default) uses generic defaults equivalent to flatfish. |
| `growth_step` | `NULL` | 每个时间步包含的年数；季度为 0.25。 | Time step in years; defaults to the species preset. |
| `burn_in` | `NULL` | 模拟预热年数，不包含在最终观测中。 | Simulation years discarded before returning observations. |
| `observation_error` | `c("independent", "shared_time")` | independent：各单元独立；shared_time：同一期共享误差。 | Independent lognormal errors per length and time cell, or shared_time for an explicitly misspecified, perfectly correlated scenario. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list containing all initialized parameters and derived variables (ages, years, etc.), ready to pass to sim_cal().

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
pa <- initialize_params(species = "flatfish", nyear = 100, burn_in = 80)
```

表中 NULL 表示从物种预设取值，英文括号中的数值是通用/flatfish 基线，并非所有物种的默认值。查看 `str(initialize_params(species="tuna"))` 得到完整季度预设。 / NULL inherits the selected preset; numeric defaults mentioned in help describe the generic/flatfish baseline.

<a id="mat_func"></a>
## `mat_func()`

计算 logistic 成熟比例 / Calculate logistic maturity

```r
mat_func(L50, L95, length)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `L50` | `必填 / required` | 成熟比例为 0.5 的体长。 | Length at which 50% of individuals are mature. |
| `L95` | `必填 / required` | 成熟比例为 0.95 的体长，须大于 L50。 | Length at which 95% of individuals are mature. |
| `length` | `必填 / required` | 计算成熟比例的体长向量。 | Length of the fish. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 Maturation probability.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
mat_func(L50 = 35, L95 = 40, length = seq(6, 50, 2))
```

<a id="plot_abundance"></a>
## `plot_abundance()`

丰度 N / Abundance N

```r
plot_abundance(model_result, line_size = 1.2, line_color = NULL, line_type = "solid",
    se = FALSE, se_color = NULL, se_alpha = 0.2, type = c("N",
        "NA", "NL"), facet_ncol = NULL, facet_scales = "free",
    return_data = FALSE)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list from run_acl or run_alscl. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Line type. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to plot confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. CI ribbon transparency. Default is 0.2. |
| `type` | `c("N", "NA", "NL")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "N" (total), "NA" (at age), or "NL" (at length). Default is "N". |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Number of columns in facet wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Scales for facet wrap. Default is "free". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Whether to return processed data. Default is FALSE. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_abundance(b, type = "N", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_abundance)。

<a id="plot_biomass"></a>
## `plot_biomass()`

生物量 B / Biomass B

```r
plot_biomass(model_result, line_size = 1.2, line_color = NULL, line_type = "solid",
    se = FALSE, se_color = NULL, se_alpha = 0.2, type = c("B",
        "BL", "BA"), facet_ncol = NULL, facet_scales = "free",
    return_data = FALSE)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list from run_acl or run_alscl. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Line type. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to plot confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. CI ribbon transparency. Default is 0.2. |
| `type` | `c("B", "BL", "BA")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "B" (total), "BL" (at length), or "BA" (at age). Default is "B". |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Integer. Columns in facet_wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Scales for facet_wrap. Default is "free". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Whether to return processed data. Default is FALSE. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_biomass(b, type = "B", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_biomass)。

<a id="plot_catch"></a>
## `plot_catch()`

捕获尾数 CN / Catch numbers CN

```r
plot_catch(model_result, line_size = 1.2, line_color = NULL, line_type = "solid",
    se = FALSE, se_color = NULL, se_alpha = 0.2, facet_ncol = NULL,
    facet_scales = "free", type = c("CN", "CNA", "CNL"), return_data = FALSE)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list from run_acl or run_alscl. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Line type. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to plot confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. CI ribbon transparency. Default is 0.2. |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Numeric. Columns in facet wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Scales for facet wrap. Default is "free". |
| `type` | `c("CN", "CNA", "CNL")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "CN" (total), "CNA" (at age), or "CNL" (at length). Default is "CN". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Whether to return processed data. Default is FALSE. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_catch(b, type = "CN", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_catch)。

<a id="plot_CatL"></a>
## `plot_CatL()`

调查拟合 length 对数尺度 / Survey fit length log scale

```r
plot_CatL(model_result, point_size = 2.5, point_color = NULL,
    point_shape = 16, line_size = 1.2, line_color = NULL, line_type = "solid",
    line_size1 = 1.8, line_color1 = NULL, line_type1 = "solid",
    line_alpha1 = 1, line_size2 = 1, line_color2 = NULL, line_type2 = "solid",
    line_alpha2 = 0.6, facet_ncol = NULL, facet_scales = "free",
    type = c("length", "year"), exp_transform = FALSE, return_data = FALSE,
    title = NULL, xlab = NULL, ylab = NULL, font_family = NULL,
    title_size = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list that contains the model output. The list should have a "report" component which contains an "Elog_index" component representing counts in length groups. |
| `point_size` | `2.5` | 点大小。 | Numeric. Specifies the size of the point in the plot. Default is 2.5. |
| `point_color` | `NULL` | 点颜色。 | Character or NULL. NULL inherits the global observed_color setting. |
| `point_shape` | `16` | 点形状编号，如 16、21。 | Numeric. Specifies the shape of the point in the plot. Default is 16 (filled circle). Common: 1=open circle, 16=filled circle, 17=triangle, 15=square. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness for type="length" plot. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Line type for type="length" plot. Default is "solid". |
| `line_size1` | `1.8` | plot_CatL 年份分面中拟合曲线的线宽/颜色/线型/透明度。 | Numeric. Line thickness for estimated (Elog_index) in type="year". Default is 1.8. |
| `line_color1` | `NULL` | plot_CatL 年份分面中拟合曲线的线宽/颜色/线型/透明度。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type1` | `"solid"` | plot_CatL 年份分面中拟合曲线的线宽/颜色/线型/透明度。 | Character. Line type for estimated line. Default is "solid". |
| `line_alpha1` | `1` | plot_CatL 年份分面中拟合曲线的线宽/颜色/线型/透明度。 | Numeric. Transparency for estimated line (0-1). Default is 1. |
| `line_size2` | `1` | plot_CatL 年份分面中观测曲线的线宽/颜色/线型/透明度。 | Numeric. Line thickness for observed (logN_at_len) in type="year". Default is 1.0. |
| `line_color2` | `NULL` | plot_CatL 年份分面中观测曲线的线宽/颜色/线型/透明度。 | Character or NULL. NULL inherits the global observed_color setting. |
| `line_type2` | `"solid"` | plot_CatL 年份分面中观测曲线的线宽/颜色/线型/透明度。 | Character. Line type for observed line. Default is "solid". |
| `line_alpha2` | `0.6` | plot_CatL 年份分面中观测曲线的线宽/颜色/线型/透明度。 | Numeric. Transparency for observed line (0-1). Default is 0.6. |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Integer. Specifies the number of columns in the facet_wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Specifies scales for facet_wrap. Default is "free". |
| `type` | `c("length", "year")` | 选择输出类型；可选值列于签名和英文说明。 | Character. It specifies whether the Elog_index is plotted across "length" or "year". Default is "length". |
| `exp_transform` | `FALSE` | 对数调查值取指数，显示原尺度；拟合曲线是中位数。 | Logical. Specifies whether to apply the exponential function to the data before plotting. If TRUE, the exponential of the data values is plotted. Default is FALSE. |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | A logical indicating whether to return the processed data alongside the plot. Default is FALSE. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot, data1 and data2.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_CatL(b, type = "length", exp_transform = FALSE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015),
  ylab = "Log survey index")
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_CatL)。

<a id="plot_compare_annual_F"></a>
## `plot_compare_annual_F()`

总体 F 比较 apical / Summary F comparison apical

```r
plot_compare_annual_F(model1, model2, model1_name = NULL, model2_name = NULL,
    method = c("apical", "mean"), colors = NULL, linetypes = NULL,
    linewidth = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `method` | `c("apical", "mean")` | apical 取组间最大 F；mean 取非加权均值，不自动年化。 | Character. How to summarize F across age/length groups each year. "apical" (default): maximum F across groups -- comparable across age-based and length-based models. "mean": arithmetic mean across all groups -- can be misleading when models differ in dimension (e.g. 20 ages vs 22 length bins with many near-zero). |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_annual_F(a, b, method = "apical") +
  labs(title = "Apical fishing mortality")
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_annual_F)。

<a id="plot_compare_CatL"></a>
## `plot_compare_CatL()`

两模型调查拟合 / Survey fits from both models

```r
plot_compare_CatL(model1, model2, data.CatL, years = NULL, model1_name = NULL,
    model2_name = NULL, colors = NULL, linetypes = NULL, linewidth = NULL,
    ncol = 4, scales = "free_y", obs_color = "grey40", obs_size = 1.5,
    obs_alpha = 0.6, obs_shape = 16, legend_pos = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Catch-at-length data frame. |
| `years` | `NULL` | 年份向量；快速模拟器要求连续年度，比较图用它筛选时期。 | Numeric vector of years to show. Default NULL = all years. Use e.g. years = seq(1991, 2011, by = 2) to select specific years. |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |
| `ncol` | `4` | 分面排列的列数。 | Integer. Facet columns. Default 4. |
| `scales` | `"free_y"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Facet scales. Default "free_y". |
| `obs_color` | `"grey40"` | 观测点颜色。 | Character. Color for observed points. Default "grey40". |
| `obs_size` | `1.5` | 观测点大小。 | Numeric. Size of observed points. Default 1.5. |
| `obs_alpha` | `0.6` | 观测点透明度。 | Numeric. Transparency of observed points. Default 0.6. |
| `obs_shape` | `16` | 观测点形状。 | Integer. Shape of observed points. Default 16. |
| `legend_pos` | `NULL` | 图例位置，如 bottom、right、none。 | Character. Legend position. NULL = use global theme. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_CatL(a, b, dat, years = c(2000, 2005, 2010, 2015), ncol = 2)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_CatL)。

<a id="plot_compare_F"></a>
## `plot_compare_F()`

死亡率矩阵比较 / Fishing mortality matrices

```r
plot_compare_F(model1, model2, model1_name = NULL, model2_name = NULL,
    palette = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `palette` | `NULL` | ggsci 色板名；NULL 继承全局设置，部分图也保留旧色板。 | Character or NULL. NULL uses the global sequential colors; a ggsci palette name or a viridis option overrides them. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_F(a, b) + scale_x_continuous(breaks = c(2000, 2010, 2019))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_F)。

<a id="plot_compare_growth"></a>
## `plot_compare_growth()`

生长与外部参照接口 / Growth and reference-data interface

```r
plot_compare_growth(model1, model2, age_range = c(1, 20), model1_name = NULL,
    model2_name = NULL, colors = NULL, linetypes = NULL, linewidth = NULL,
    ref_data = NULL, ref_name = "Observed VB", ref_color = "grey30",
    ref_linetype = "dotdash", show_points = TRUE, point_size = NULL,
    point_alpha = 0.4, nls_start = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `age_range` | `c(1, 20)` | 生长图的年龄上下限，单位为年。 | Numeric vector of length 2. Default c(1, 20). |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |
| `ref_data` | `NULL` | 可选 Age、Length 列的参考个体生长数据。 | Optional data.frame with columns Age and Length containing observed age-length measurements (e.g. from otolith readings). If provided, a VB curve is fitted via nls() and overlaid as a reference, with the raw data as scatter points. |
| `ref_name` | `"Observed VB"` | 参考数据曲线的图例名称。 | Character. Label for the reference curve. Default "Observed VB". |
| `ref_color` | `"grey30"` | 参考曲线颜色。 | Character. Color for reference curve/points. Default "grey30". |
| `ref_linetype` | `"dotdash"` | 参考曲线线型。 | Character. Linetype for reference curve. Default "dotdash". |
| `show_points` | `TRUE` | 是否显示参考数据点。 | Logical. Show raw data points from ref_data. Default TRUE. |
| `point_size` | `NULL` | 点大小。 | Numeric. Size of data points. NULL = use global theme. |
| `point_alpha` | `0.4` | 点透明度，0–1。 | Numeric. Transparency of data points. Default 0.4. |
| `nls_start` | `NULL` | 参考数据 VB 非线性回归初值：Linf、k、t0。 | Named list. Starting values for nls. Default attempts auto-detection from data range. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
growth_display(plot_compare_growth(a, b, age_range = c(1, 15), ref_data = ref, ref_name = "Synthetic reference",
    nls_start = list(Linf = 60, k = 0.2, t0 = 1/60)))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_growth)。

如使用 `ref_data`，请先执行图谱中的模拟参考数据代码；示例的 `growth_display` 是脚本排版辅助函数，不属于 API。 / Define the synthetic reference data and layout helper from the gallery before using them.

<a id="plot_compare_metrics"></a>
## `plot_compare_metrics()`

拟合指标比较 / Fit metric comparison

```r
plot_compare_metrics(model1, model2, data.CatL, model1_name = NULL, model2_name = NULL,
    colors = NULL, metrics = NULL, ncol = NULL, legend_pos = NULL,
    label_size = 3.2, bar_width = 0.6, bar_alpha = 0.85)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Catch-at-length data frame. |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `metrics` | `NULL` | 选择 error、R2、MAPE、IC、npar 面板。 | Character vector selecting which panels to show. Options: "error" (RMSE, MAE), "R2", "MAPE", "IC" (AIC, BIC, NLL), "npar". Default NULL = all five panels. E.g. metrics = c("error", "R2", "IC"). |
| `ncol` | `NULL` | 分面排列的列数。 | Integer. Columns for panel layout. Default: auto. |
| `legend_pos` | `NULL` | 图例位置，如 bottom、right、none。 | Character. Legend position. NULL = use global theme. |
| `label_size` | `3.2` | 柱图数值标签字号。 | Numeric. Text label size on bars. Default 3.2. |
| `bar_width` | `0.6` | 柱子宽度。 | Numeric. Bar width. Default 0.6. |
| `bar_alpha` | `0.85` | 柱子透明度。 | Numeric. Bar transparency. Default 0.85. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_metrics(a, b, dat, ncol = 2) & scale_y_continuous(expand = expansion(mult = c(0,
    0.3)))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_metrics)。

<a id="plot_compare_residuals"></a>
## `plot_compare_residuals()`

残差联合诊断 / Combined residual diagnostics

```r
plot_compare_residuals(model1, model2, data.CatL, model1_name = NULL, model2_name = NULL,
    colors = NULL, linetypes = NULL, linewidth = NULL, point_size = NULL,
    legend_pos = NULL, bins = 35, qq_alpha = 0.6, errorbar_width = 0.3)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Catch-at-length data frame. |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |
| `point_size` | `NULL` | 点大小。 | Numeric. Point size. NULL = use global theme. |
| `legend_pos` | `NULL` | 图例位置，如 bottom、right、none。 | Character. Legend position. NULL = use global theme (default "bottom"). |
| `bins` | `35` | 残差直方图分箱数量。 | Integer. Histogram bins. Default 35. |
| `qq_alpha` | `0.6` | QQ 图点的透明度。 | Numeric. QQ point transparency. Default 0.6. |
| `errorbar_width` | `0.3` | 年度残差 SD 误差棒端帽宽度。 | Numeric. Error bar cap width. Default 0.3. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object (combined 4-panel via patchwork).

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_residuals(a, b, dat)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_residuals)。

<a id="plot_compare_selectivity"></a>
## `plot_compare_selectivity()`

固定调查可捕性 / Fixed survey catchability

```r
plot_compare_selectivity(model1, model2, model1_name = NULL, model2_name = NULL,
    colors = NULL, linetypes = NULL, linewidth = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First fit; its fixed survey catchability is plotted. |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_selectivity(a, b)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_selectivity)。

<a id="plot_compare_ts"></a>
## `plot_compare_ts()`

种群时间序列比较 / Population time series comparison

```r
plot_compare_ts(model1, model2, quantities = c("SSB", "R", "B", "N",
    "CN", "CB"), se = FALSE, model1_name = NULL, model2_name = NULL,
    colors = NULL, linetypes = NULL, linewidth = NULL, ncol = NULL,
    scales = NULL, return_data = FALSE)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model1` | `必填 / required` | 第一个拟合对象。 | First model result (reference). |
| `model2` | `必填 / required` | 第二个拟合对象。 | Second model result (alternative). |
| `quantities` | `c("SSB", "R", "B", "N", "CN", "CB")` | 选定输出：SSB、R/Rec、B、N、CN、CB、F。 | Character vector. Quantities to compare. Accepts: "SSB", "R" (or "Rec"), "B", "N", "CN", "CB", "F". Default c("SSB","R","B","N","CN","CB"). |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Show confidence intervals. Default FALSE. |
| `model1_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `model2_name` | `NULL` | 图例中的模型名称；NULL 自动判断。 | Character. Labels. Auto-detected if NULL. |
| `colors` | `NULL` | 按两个模型顺序提供两种颜色；名称可选。 | Named character vector of 2 colors. NULL = use global theme. |
| `linetypes` | `NULL` | 按模型顺序提供两种线型。 | Named character vector of 2 linetypes. NULL = use global theme. |
| `linewidth` | `NULL` | 线宽；参数名按函数区别使用。 | Numeric. NULL = use global theme. |
| `ncol` | `NULL` | 分面排列的列数。 | Integer. Facet columns. NULL = use global theme. |
| `scales` | `NULL` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Facet scales. NULL = use global theme. |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Return data alongside plot. Default FALSE. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_compare_ts(a, b, se = TRUE, ncol = 2)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_compare_ts)。

<a id="plot_deviance"></a>
## `plot_deviance()`

过程偏差 R TRUE / Process deviations R TRUE

```r
plot_deviance(model_result, se = TRUE, point_size = 3, point_color = "white",
    point_shape = 21, line_size = 1, line_color = NULL, line_type = "solid",
    se_color = NULL, se_width = 0.5, facet_ncol = NULL, facet_scales = "free",
    log = TRUE, type = c("R", "F"))
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list from run_acl or run_alscl. |
| `se` | `TRUE` | 是否绘制可用标准误对应的近似区间。 | Logical, whether to draw the standard error bar. |
| `point_size` | `3` | 点大小。 | Numeric. The size of the point. Default is 3. |
| `point_color` | `"white"` | 点颜色。 | Character. The color of the point. Default is "white". |
| `point_shape` | `21` | 点形状编号，如 16、21。 | Numeric. The shape of the point. Default is 21. |
| `line_size` | `1` | 线宽；参数名按函数区别使用。 | Numeric, the line size. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character, the line type. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_width` | `0.5` | 误差棒端帽宽度。 | Numeric, the width of the standard error. |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Number of columns in facet wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Scales for facet wrap. Default is "free". |
| `log` | `TRUE` | TRUE 为对数过程偏差；FALSE 取指数，参考值为 1。 | Logical, whether to keep log scale (TRUE) or apply exp transform (FALSE). |
| `type` | `c("R", "F")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "R" for recruitment deviation, "F" for fishing mortality deviation. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot2 object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_deviance(b, type = "R", log = TRUE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_deviance)。

<a id="plot_fishing_mortality"></a>
## `plot_fishing_mortality()`

ACL 捕捞死亡率 year / ACL fishing mortality year

```r
plot_fishing_mortality(model_result, line_size = 1, line_color = NULL, line_type = "solid",
    facet_ncol = NULL, facet_scales = "free", se = FALSE, se_color = NULL,
    se_alpha = 0.2, se_type = "ribbon", type = c("year", "age",
        "length"), return_data = FALSE, x_breaks = NULL, title = NULL,
    xlab = NULL, ylab = NULL, font_family = NULL, title_size = NULL,
    base_theme = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list obtained from run_acl or run_alscl. |
| `line_size` | `1` | 线宽；参数名按函数区别使用。 | Numeric. The thickness of the line. Default is 1. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. The type of the line. Default is "solid". |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Integer. The number of columns in facet_wrap. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Scales for facet_wrap. Default is "free". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to plot confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. The transparency of the ribbon. Default is 0.2. |
| `se_type` | `"ribbon"` | 区间样式：ribbon 或 errorbar。 | Character. "ribbon" or "errorbar". Default is "ribbon". |
| `type` | `c("year", "age", "length")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "year" (faceted by group), "age" (faceted by year), or "length" (ALSCL: faceted by year showing F vs length). Default is "year". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Whether to return data alongside the plot. Default is FALSE. |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector. Custom x-axis breaks. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character. Custom plot title. |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character. Custom x-axis label. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character. Custom y-axis label. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character. Font family for text elements. |
| `title_size` | `NULL` | 图标题字号。 | Numeric. Title text size. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | A ggplot2 theme object. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_fishing_mortality(a, type = "year",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_fishing_mortality)。

ACL 的 `type="length"` 当前与 `"age"` 是兼容别名，不能解释为长度别 F。 / For ACL, length remains an age-plot compatibility alias.

<a id="plot_pla"></a>
## `plot_pla()`

年龄长度转换 / Age length conversion

```r
plot_pla(model_result, low_col = NULL, high_col = NULL, title = NULL,
    xlab = NULL, ylab = NULL, font_family = NULL, title_size = NULL,
    axis_title_size = NULL, axis_text_size = NULL, strip_text_size = NULL,
    legend_text_size = NULL, x_breaks = NULL, base_theme = NULL,
    title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list that contains model output. The model output should include a PLA report. |
| `low_col` | `NULL` | 概率热图的低值 / 高值颜色。 | Character or NULL. NULL inherits the global low_col setting. |
| `high_col` | `NULL` | 概率热图的低值 / 高值颜色。 | Character or NULL. NULL inherits the global high_col setting. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object showing the heatmap of the probability length at age.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_pla(b) + scale_x_discrete(labels = as.character(1:15))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_pla)。

<a id="plot_recruitment"></a>
## `plot_recruitment()`

补充量 / Recruitment

```r
plot_recruitment(model_result, line_size = 1.5, line_color = NULL, line_type = "solid",
    se = FALSE, se_color = NULL, se_alpha = 0.2, se_type = c("ribbon",
        "errorbar"), return_data = FALSE, title = NULL, xlab = NULL,
    ylab = NULL, font_family = NULL, title_size = NULL, axis_title_size = NULL,
    axis_text_size = NULL, strip_text_size = NULL, legend_text_size = NULL,
    x_breaks = NULL, base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list that contains the model output. The list should have a "report" component which contains a "Rec" component representing Recruitment. |
| `line_size` | `1.5` | 线宽；参数名按函数区别使用。 | Numeric. Specifies the thickness of the line in the plot. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Specifies the type of the line in the plot. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Determines whether to calculate and plot the standard error as confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. The transparency of the confidence interval ribbon. Default is 0.2. |
| `se_type` | `c("ribbon", "errorbar")` | 区间样式：ribbon 或 errorbar。 | Character. Type of CI display: "ribbon" (shaded area) or "errorbar" (error bars). Default is "ribbon". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | A logical indicating whether to return the processed data alongside the plot. Default is FALSE. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_recruitment(b, se = TRUE)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_recruitment)。

<a id="plot_residuals"></a>
## `plot_residuals()`

对数残差 length / Log residuals length

```r
plot_residuals(model_result, f = 0.4, line_color = NULL, smooth_color = NULL,
    hline_color = NULL, line_size = 1, facet_scales = "free",
    facet_ncol = NULL, type = c("length", "year"), resid_cap = NULL,
    return_data = FALSE, title = NULL, xlab = NULL, ylab = NULL,
    font_family = NULL, title_size = NULL, axis_title_size = NULL,
    axis_text_size = NULL, strip_text_size = NULL, legend_text_size = NULL,
    x_breaks = NULL, base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list obtained from the run_acl function. This list should contain a "report" component that includes a "resid_index" component, representing the residuals for each index. |
| `f` | `0.4` | 残差平滑参数，越大越平滑。 | Numeric. The smoother span for the loess smooth line in the plot. This gives the proportion of points in the plot which influence the smooth at each value. Larger values result in more smoothing. Default is 0.4. Sparse facets use a linear trend; otherwise the span is enlarged when needed for at least five neighbors. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `smooth_color` | `NULL` | 残差平滑线颜色。 | Character or NULL. NULL inherits the global smooth_color setting. |
| `hline_color` | `NULL` | 残差零参考线颜色。 | Character or NULL. NULL inherits the global hline_color setting. |
| `line_size` | `1` | 线宽；参数名按函数区别使用。 | Numeric. The size of the lines in the plot. Default is 1. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. The "scales" argument for the facet_wrap function in ggplot2. Default is "free". |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Numeric, optional. Number of columns in facet wrap. Default is NULL. |
| `type` | `c("length", "year")` | 选择输出类型；可选值列于签名和英文说明。 | Character. Specifies whether the residuals are calculated over "length" or "year". Default is "length". |
| `resid_cap` | `NULL` | 可选残差显示裁剪阈值；使用时需说明。 | Numeric or NULL. Symmetrically cap residuals at +/- this value. Useful for species with many length bins where tail bins (e.g. >120 for tuna) have very few observations, producing extreme log-residuals that distort the plot. Default is NULL (no capping). A value of 1 is a reasonable starting point. |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | A logical indicating whether to return the processed data alongside the plot. Default is FALSE. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_residuals(b, type = "length", facet_ncol = 4, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015))
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_residuals)。

<a id="plot_retro"></a>
## `plot_retro()`

ACL 回溯分析 / ACL retrospective analysis

```r
plot_retro(retro_result, rho_digits = 4, rho_position = "top_right",
    rho_size = 3.5, line_size = 1.2, point_size = 3, point_shape = 21,
    facet_scales = "free", facet_col = NULL, facet_row = NULL,
    title = NULL, xlab = NULL, ylab = NULL, font_family = NULL,
    title_size = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    base_theme = NULL, title_hjust = NULL, palette = NULL, colors = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `retro_result` | `必填 / required` | retro_model/retro_acl/retro_alscl 的回溯结果。 | A list returned by retro_model(), retro_acl(), or retro_alscl(). |
| `rho_digits` | `4` | Mohn rho 显示的小数位数。 | Integer. Decimal places for Mohn's rho. Default is 4. |
| `rho_position` | `"top_right"` | Mohn rho 标注位置。 | Character. Position of rho text. Default is "top_right". |
| `rho_size` | `3.5` | Mohn rho 标注字号。 | Numeric. Font size of rho text. Default is 3.5. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `point_size` | `3` | 点大小。 | Numeric. Point size. Default is 3. |
| `point_shape` | `21` | 点形状编号，如 16、21。 | Numeric. Point shape. Default is 21. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Facet scales. Default is "free". |
| `facet_col` | `NULL` | 分面排列的列数。 | Integer or NULL. Number of facet columns. |
| `facet_row` | `NULL` | 分面排列的行数。 | Integer or NULL. Number of facet rows. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Font family. |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Title size in pt. |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment. |
| `palette` | `NULL` | ggsci 色板名；NULL 继承全局设置，部分图也保留旧色板。 | Character or NULL. ggsci palette name; NULL inherits acl_theme(). Native colors are interpolated when there are more peels than colors. |
| `colors` | `NULL` | 按截止期降序给色，或按截止期命名；需覆盖每条曲线。 | Character vector or NULL. Explicit period colors, ordered by descending terminal period or named by terminal period. Overrides palette. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_retro(rb, facet_col = 2, rho_digits = 3)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_retro)。

<a id="plot_ridges"></a>
## `plot_ridges()`

长度组成山脊图 / Length composition ridges

```r
plot_ridges(model_result, ridges_alpha = 0.8, ridges_scale = NULL,
    palette = "viridis", title = NULL, xlab = NULL, ylab = NULL, font_family = NULL,
    title_size = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list containing the model result, which includes a data frame of Elog_index and logN_at_len, len_mid, and year. |
| `ridges_alpha` | `0.8` | 山脊图透明度。 | Numeric. Specifies the transparency of the ridgeline plots. Default is 0.8. |
| `ridges_scale` | `NULL` | 山脊图高度/重叠缩放。 | Numeric or NULL. Multiplier for ridge height. When NULL (default), auto-scales so the tallest peak fills ~80\% of inter-year spacing. Increase for taller ridges, decrease for flatter. Useful when many length bins make proportions small (e.g. tuna with 22 bins). |
| `palette` | `"viridis"` | 按时间顺序渐变；NULL 同样使用 viridis，不继承全局分类色板。可选 cividis、plasma、ocean 或自定义渐变端点。 | Sequential gradient across ordered years; NULL also uses viridis independently of the global categorical palette. Accepts other named gradients or a color vector; explicit journal palettes remain available. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A grid arranged ggplot object.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_ridges(b)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_ridges)。

<a id="plot_SSB"></a>
## `plot_SSB()`

产卵生物量 SSB / Spawning biomass SSB

```r
plot_SSB(model_result, line_size = 1.2, line_color = NULL, line_type = "solid",
    se = FALSE, se_color = NULL, se_alpha = 0.2, type = c("SSB",
        "SBL", "SBA"), facet_ncol = NULL, facet_scales = "free",
    return_data = FALSE)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list from run_acl or run_alscl. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. Line type. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to plot confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. CI ribbon transparency. Default is 0.2. |
| `type` | `c("SSB", "SBL", "SBA")` | 选择输出类型；可选值列于签名和英文说明。 | Character. "SSB" (total), "SBL" (at length), or "SBA" (at age). Default is "SSB". |
| `facet_ncol` | `NULL` | 分面排列的列数。 | Numeric. Columns in facet wrap. Default is NULL. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Scales for facet wrap. Default is "free". |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | Logical. Whether to return processed data. Default is FALSE. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_SSB(b, type = "SSB", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_SSB)。

<a id="plot_SSB_Rec"></a>
## `plot_SSB_Rec()`

亲体与补充关系 / Spawning biomass and recruitment

```r
plot_SSB_Rec(model_result, age_at_recruitment = 1, point_size = 2,
    point_color = NULL, point_shape = 16, return_data = FALSE,
    title = NULL, xlab = NULL, ylab = NULL, font_family = NULL,
    title_size = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list containing the model results. It includes a component "report" that comprises "SSB" and "Rec". |
| `age_at_recruitment` | `1` | 亲体补充图的错位步数；单位是观测步，不自动转换为年。 | Numeric. The age at which individuals are recruited. Used to align SSB with Rec. Default is 1. |
| `point_size` | `2` | 点大小。 | Numeric. The size of the points in the scatter plot. Default is 2. |
| `point_color` | `NULL` | 点颜色。 | Character or NULL. NULL inherits the global point_color setting. |
| `point_shape` | `16` | 点形状编号，如 16、21。 | Numeric. The shape of the points in the scatter plot, as an integer value (see ?points in base R for more info). Default is 16 (filled circle). |
| `return_data` | `FALSE` | TRUE 同时返回绘图对象及整理后的数据；结构见返回说明。 | A logical indicating whether to return the processed data alongside the plot. Default is FALSE. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot; return_data=TRUE returns a list with plot and data.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_SSB_Rec(b, age_at_recruitment = 1)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_SSB_Rec)。

<a id="plot_VB"></a>
## `plot_VB()`

生长曲线与区间 / Growth curve and interval

```r
plot_VB(model_result, age_range = c(1, 25), line_size = 1.2,
    line_color = NULL, line_type = "solid", se = FALSE, se_color = NULL,
    se_alpha = 0.2, se_type = c("ribbon", "errorbar"), text_color = "black",
    text_size = 5, title = NULL, xlab = NULL, ylab = NULL, font_family = NULL,
    title_size = NULL, axis_title_size = NULL, axis_text_size = NULL,
    strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL,
    base_theme = NULL, title_hjust = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `model_result` | `必填 / required` | 单次模型拟合结果。 | A list that contains model output. The list should have a "report" component which contains "Linf", "vbk" and "t0" components. |
| `age_range` | `c(1, 25)` | 生长图的年龄上下限，单位为年。 | Numeric vector of length 2, defining the range of ages to consider. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. The thickness of the line in the plot. Default is 1.2. |
| `line_color` | `NULL` | 单模型线条颜色。 | Character or NULL. NULL inherits the global line_color setting. |
| `line_type` | `"solid"` | 线型，如 solid、dashed。 | Character. The type of the line in the plot. Default is "solid". |
| `se` | `FALSE` | 是否绘制可用标准误对应的近似区间。 | Logical. Whether to calculate and plot standard error as confidence intervals. Default is FALSE. |
| `se_color` | `NULL` | 区间阴影或误差棒的颜色。 | Character or NULL. NULL inherits the global se_color setting. |
| `se_alpha` | `0.2` | 区间透明度，0–1。 | Numeric. The transparency of the confidence interval ribbon. Default is 0.2. |
| `se_type` | `c("ribbon", "errorbar")` | 区间样式：ribbon 或 errorbar。 | Character. Type of CI display: "ribbon" (shaded area) or "errorbar" (error bars). Default is "ribbon". |
| `text_color` | `"black"` | 生长参数标注文字颜色。 | Character. The color of the text in the plot. Default is "black". |
| `text_size` | `5` | 生长参数标注文字大小。 | Numeric. The thickness of the text in the plot. Default is 5. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. If NULL, uses global theme setting. See acl_theme_set(). |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. If NULL, uses global theme setting. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. If NULL, uses global theme setting. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Custom font family. If NULL, uses global theme setting (default "sans"). |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14). |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12). |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10). |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10). |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10). |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks (e.g. seq(1, 20, by = 2)). NULL = auto. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting. |

**返回 / Returns:** 图对象；如支持 return_data，可同时取回整理数据。 A ggplot object representing the plot.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
plot_VB(b, age_range = c(1, 15), se = TRUE)
```

[查看实际 YTF 图和解读 / Rendered YTF examples](PLOT_GALLERY.md#plot_VB)。

<a id="retro_acl"></a>
## `retro_acl()`

ACL 回溯入口 / ACL retrospective wrapper

```r
retro_acl(nyear = 5, data.CatL, data.wgt, data.mat, rec.age,
    nage, M, sel_L50, sel_L95, parameters = NULL, parameters.L = NULL,
    parameters.U = NULL, map = NULL, len_mid = NULL, len_border = NULL,
    plot = FALSE, line_size = 1.2, point_size = 3, point_shape = 21,
    facet_scales = "free", facet_col = NULL, facet_row = NULL,
    train_times = 1, title = NULL, xlab = NULL, ylab = NULL,
    font_family = NULL, title_size = NULL, axis_title_size = NULL,
    axis_text_size = NULL, strip_text_size = NULL, legend_text_size = NULL,
    x_breaks = NULL, base_theme = NULL, title_hjust = NULL, rho_digits = 4,
    rho_position = "top_right", rho_size = 3.5, ncores = 1, ...)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `nyear` | `5` | 模拟器为总年数；回溯函数为删除的末端时间步数。 | Number of terminal observation steps peeled, not calendar years for quarterly data. |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | The Catch-at-Length data. |
| `data.wgt` | `必填 / required` | 平均单尾体重表，与调查表完全对齐。 | The weight data. |
| `data.mat` | `必填 / required` | 成熟比例表，0–1，与调查表完全对齐。 | The maturity data. |
| `rec.age` | `必填 / required` | 补充年龄，单位为年。 | The recruitment age. |
| `nage` | `必填 / required` | 年龄组数量，包括最大年龄加组。 | The number of age classes. |
| `M` | `必填 / required` | 每个模型时间步的自然死亡瞬时率。 | The natural mortality rate. |
| `sel_L50` | `必填 / required` | 调查可捕性达到 50% 的体长。 | Length at 50 percent selectivity. |
| `sel_L95` | `必填 / required` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Length at 95 percent selectivity. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | Optional custom starting parameter list. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | Optional lower bounds for parameters. |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | Optional upper bounds for parameters. |
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | Optional TMB map argument for fixing parameters. |
| `len_mid` | `NULL` | 每个体长组的中心值向量，长度与数据行数相同。 | Optional numeric vector of length bin midpoints. |
| `len_border` | `NULL` | 体长边界；拟合入口需要 nlen-1 个内部边界，模拟器使用 nlen+1 个完整边界。 | Optional numeric vector of length bin borders. |
| `plot` | `FALSE` | TRUE 在回溯返回对象中附带图。 | Logical. If TRUE, include plot in return value. Default is FALSE. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `point_size` | `3` | 点大小。 | Numeric. Point size for endpoints. Default is 3. |
| `point_shape` | `21` | 点形状编号，如 16、21。 | Numeric. Point shape for endpoints. Default is 21. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Facet scales argument. Default is "free". |
| `facet_col` | `NULL` | 分面排列的列数。 | Integer or NULL. Number of facet columns. |
| `facet_row` | `NULL` | 分面排列的行数。 | Integer or NULL. Number of facet rows. |
| `train_times` | `1` | 每个初始点连续优化的次数。 | Number of successive optimization passes per start; not independent random starts. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Font family. |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Title size in pt. |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment. |
| `rho_digits` | `4` | Mohn rho 显示的小数位数。 | Integer. Decimal places for Mohn's rho display. Default is 4. |
| `rho_position` | `"top_right"` | Mohn rho 标注位置。 | Character. Position of rho text. Default is "top_right". |
| `rho_size` | `3.5` | Mohn rho 标注字号。 | Numeric. Font size of rho text. Default is 3.5. |
| `ncores` | `1` | 并行数量；具体含义依函数而异，见英文说明。 | Integer. Number of CPU cores for parallel peel computation. Default is 1. |
| `...` | `必填 / required` | 转发参数；仅使用目标函数接受的命名参数。 | Additional arguments passed to retro_model(). |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list. See retro_model().

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
ra <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl, list(nyear = 3)))
```

<a id="retro_alscl"></a>
## `retro_alscl()`

ALSCL 回溯入口 / ALSCL retrospective wrapper

```r
retro_alscl(nyear = 5, data.CatL, data.wgt, data.mat, rec.age,
    nage, M, sel_L50, sel_L95, growth_step = 1, ...)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `nyear` | `5` | 模拟器为总年数；回溯函数为删除的末端时间步数。 | Number of terminal observation steps peeled, not calendar years for quarterly data. |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | The Catch-at-Length data. |
| `data.wgt` | `必填 / required` | 平均单尾体重表，与调查表完全对齐。 | The weight data. |
| `data.mat` | `必填 / required` | 成熟比例表，0–1，与调查表完全对齐。 | The maturity data. |
| `rec.age` | `必填 / required` | 补充年龄，单位为年。 | The recruitment age. |
| `nage` | `必填 / required` | 年龄组数量，包括最大年龄加组。 | The number of age classes. |
| `M` | `必填 / required` | 每个模型时间步的自然死亡瞬时率。 | The natural mortality rate. |
| `sel_L50` | `必填 / required` | 调查可捕性达到 50% 的体长。 | Length at 50 percent selectivity. |
| `sel_L95` | `必填 / required` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Length at 95 percent selectivity. |
| `growth_step` | `1` | 每个时间步包含的年数；季度为 0.25。 | Numeric. Growth transition time step. Default is 1. |
| `...` | `必填 / required` | 转发参数；仅使用目标函数接受的命名参数。 | Additional arguments passed to retro_model(). |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list. See retro_model().

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
rb <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl, list(nyear = 3)))
```

<a id="retro_model"></a>
## `retro_model()`

统一回溯分析 / Unified retrospective analysis

```r
retro_model(nyear = 5, data.CatL, data.wgt, data.mat, rec.age,
    nage, M, sel_L50, sel_L95, model_type = c("acl", "alscl"),
    growth_step = NULL, parameters = NULL, parameters.L = NULL,
    parameters.U = NULL, map = NULL, len_mid = NULL, len_border = NULL,
    len_lower = NULL, len_upper = NULL, train_times = 1, ncores = 1,
    silent = FALSE, plot = FALSE, line_size = 1.2, point_size = 3,
    point_shape = 21, facet_scales = "free", facet_col = NULL,
    facet_row = NULL, title = NULL, xlab = NULL, ylab = NULL,
    font_family = NULL, title_size = NULL, axis_title_size = NULL,
    axis_text_size = NULL, strip_text_size = NULL, legend_text_size = NULL,
    x_breaks = NULL, base_theme = NULL, title_hjust = NULL, rho_digits = 4,
    rho_position = "top_right", rho_size = 3.5)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `nyear` | `5` | 模拟器为总年数；回溯函数为删除的末端时间步数。 | Number of terminal observation steps peeled, not calendar years for quarterly data. |
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | The Catch-at-Length data. |
| `data.wgt` | `必填 / required` | 平均单尾体重表，与调查表完全对齐。 | The weight data. |
| `data.mat` | `必填 / required` | 成熟比例表，0–1，与调查表完全对齐。 | The maturity data. |
| `rec.age` | `必填 / required` | 补充年龄，单位为年。 | The recruitment age. |
| `nage` | `必填 / required` | 年龄组数量，包括最大年龄加组。 | The number of age classes. |
| `M` | `必填 / required` | 每个模型时间步的自然死亡瞬时率。 | The natural mortality rate. |
| `sel_L50` | `必填 / required` | 调查可捕性达到 50% 的体长。 | The length at 50 percent selectivity. |
| `sel_L95` | `必填 / required` | 调查可捕性达到 95% 的体长，必须大于 L50。 | The length at 95 percent selectivity. |
| `model_type` | `c("acl", "alscl")` | 模型选择：acl 或 alscl。 | Character. Model to use: "acl" (default) or "alscl". |
| `growth_step` | `NULL` | 每个时间步包含的年数；季度为 0.25。 | Time step in years. NULL uses recruitment age when below 1, otherwise 1. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | Optional custom starting parameter list. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | Optional lower bounds for parameters. |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | Optional upper bounds for parameters. |
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | Optional TMB map argument for fixing parameters. |
| `len_mid` | `NULL` | 每个体长组的中心值向量，长度与数据行数相同。 | Optional numeric vector of length bin midpoints. |
| `len_border` | `NULL` | 体长边界；拟合入口需要 nlen-1 个内部边界，模拟器使用 nlen+1 个完整边界。 | Optional numeric vector of length bin borders. |
| `len_lower` | `NULL` | ALSCL 每组下界/上界向量，须与分组一致。 | Optional numeric vector of lower bin boundaries (ALSCL only). |
| `len_upper` | `NULL` | ALSCL 每组下界/上界向量，须与分组一致。 | Optional numeric vector of upper bin boundaries (ALSCL only). |
| `train_times` | `1` | 每个初始点连续优化的次数。 | Number of successive optimization passes per start; not independent random starts. |
| `ncores` | `1` | 并行数量；具体含义依函数而异，见英文说明。 | Integer. Number of CPU cores for parallel peel computation. Default is 1. |
| `silent` | `FALSE` | TRUE 减少拟合进度输出。 | Logical. If TRUE, suppress all progress messages. Default is FALSE. |
| `plot` | `FALSE` | TRUE 在回溯返回对象中附带图。 | Logical. If TRUE, include plot in return value. Default is FALSE. |
| `line_size` | `1.2` | 线宽；参数名按函数区别使用。 | Numeric. Line thickness. Default is 1.2. |
| `point_size` | `3` | 点大小。 | Numeric. Point size for endpoints. Default is 3. |
| `point_shape` | `21` | 点形状编号，如 16、21。 | Numeric. Point shape for endpoints. Default is 21. |
| `facet_scales` | `"free"` | fixed、free、free_x 或 free_y；自由尺度有助观察组内变化。 | Character. Facet scales argument. Default is "free". |
| `facet_col` | `NULL` | 分面排列的列数。 | Integer or NULL. Number of facet columns. |
| `facet_row` | `NULL` | 分面排列的行数。 | Integer or NULL. Number of facet rows. |
| `title` | `NULL` | 自定义图标题；NULL 继承默认。 | Character or NULL. Custom plot title. |
| `xlab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom x-axis label. |
| `ylab` | `NULL` | 坐标轴标签；主题设置函数中使用命名列表。 | Character or NULL. Custom y-axis label. |
| `font_family` | `NULL` | 文本字体；请使用本机可用字体。 | Character or NULL. Font family. |
| `title_size` | `NULL` | 图标题字号。 | Numeric or NULL. Title size in pt. |
| `axis_title_size` | `NULL` | 坐标轴标题字号。 | Numeric or NULL. Axis title size in pt. |
| `axis_text_size` | `NULL` | 坐标刻度文字字号。 | Numeric or NULL. Axis tick label size in pt. |
| `strip_text_size` | `NULL` | 分面标签字号。 | Numeric or NULL. Facet label size in pt. |
| `legend_text_size` | `NULL` | 图例文字字号。 | Numeric or NULL. Legend text size in pt. |
| `x_breaks` | `NULL` | 横轴刻度向量；NULL 自动。 | Numeric vector or NULL. Custom x-axis breaks. |
| `base_theme` | `NULL` | 基础主题名称，如 bw、minimal 或 theme_bw。 | Character or NULL. Base ggplot2 theme name. |
| `title_hjust` | `NULL` | 标题水平对齐：0 左、0.5 中、1 右。 | Numeric or NULL. Title horizontal alignment. |
| `rho_digits` | `4` | Mohn rho 显示的小数位数。 | Integer. Decimal places for Mohn's rho display. Default is 4. |
| `rho_position` | `"top_right"` | Mohn rho 标注位置。 | Character. Position of rho text. Default is "top_right". |
| `rho_size` | `3.5` | Mohn rho 标注字号。 | Numeric. Font size of rho text. Default is 3.5. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list containing results, rho_text, last_points, model_type, and optionally plot.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
rb <- do.call(retro_model, c(inputs, x$fit_args, x$fit_config$alscl,
                            list(model_type = "alscl", nyear = 3)))
```

<a id="run_acl"></a>
## `run_acl()`

拟合年龄型 ACL / Fit age-based ACL

```r
run_acl(data.CatL, data.wgt, data.mat, rec.age, nage, M, sel_L50,
    sel_L95, parameters = NULL, parameters.L = NULL, parameters.U = NULL,
    map = NULL, len_mid = NULL, len_border = NULL, output = FALSE,
    train_times = 1, ncores = 1, silent = FALSE, growth_step = NULL,
    zero_action = c("missing", "error"), control = list(), nstarts = ncores)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Data frame or matrix with length labels in its first column, followed by at least two numerically labelled time columns. |
| `data.wgt` | `必填 / required` | 平均单尾体重表，与调查表完全对齐。 | Weight-at-length, with identical labels, columns and dimensions. |
| `data.mat` | `必填 / required` | 成熟比例表，0–1，与调查表完全对齐。 | Maturity-at-length between 0 and 1, with the same layout. |
| `rec.age` | `必填 / required` | 补充年龄，单位为年。 | Recruitment age in years, greater than t0. |
| `nage` | `必填 / required` | 年龄组数量，包括最大年龄加组。 | Number of age classes including the plus group (at least two). |
| `M` | `必填 / required` | 每个模型时间步的自然死亡瞬时率。 | Nonnegative natural mortality per model time step. |
| `sel_L50` | `必填 / required` | 调查可捕性达到 50% 的体长。 | Lengths at 50 and 95 percent survey catchability. |
| `sel_L95` | `必填 / required` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Lengths at 50 and 95 percent survey catchability. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | Named list of initial fixed-effect values. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | Named lower and upper bounds. Bounds are aligned with free parameters after applying map. Initial free values are clamped to bounds. |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | Named lower and upper bounds. Bounds are aligned with free parameters after applying map. Initial free values are clamped to bounds. |
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | Named TMB factor mappings. factor(NA) fixes a parameter at its supplied value; NULL releases a default fixed parameter. Other defaults remain in force. |
| `len_mid` | `NULL` | 每个体长组的中心值向量，长度与数据行数相同。 | Finite length-bin midpoints; overrides automatic parsing. |
| `len_border` | `NULL` | 体长边界；拟合入口需要 nlen-1 个内部边界，模拟器使用 nlen+1 个完整边界。 | Interior boundaries, one fewer than the number of length bins. |
| `output` | `FALSE` | 是否导出拟合诊断与图；TRUE 使用 output 目录。 | Save diagnostic tables and plots below output/ when TRUE. |
| `train_times` | `1` | 每个初始点连续优化的次数。 | Number of successive optimization passes per start; not independent random starts. |
| `ncores` | `1` | 独立初始点拟合的最大并行进程数，实际不超过 nstarts；不是单个 TMB 拟合的线程数。 | Maximum socket workers across independent starts, capped at nstarts; not threads within one start. |
| `nstarts` | `ncores` | 初始点总数；默认等于 ncores 以兼容旧用法。固定此值可公平比较不同并行数；后续起点仅扰动自由参数。 | Total independent starts; defaults to ncores for compatibility. Fix this value for equal-work comparisons. Later starts deterministically jitter free parameters only. |
| `silent` | `FALSE` | TRUE 减少拟合进度输出。 | Suppress fitting progress messages when TRUE. |
| `growth_step` | `NULL` | 每个时间步包含的年数；季度为 0.25。 | Years per model step. ACL NULL uses rec.age when below 1, otherwise 1; ALSCL defaults to 1. Supply explicitly for quarterly data. |
| `zero_action` | `c("missing", "error")` | missing 排除调查零值；error 报错。NA 始终为缺失。 | Exclude zero observations as missing (historical behavior), or reject them with "error". NA is always missing; negative/infinite values fail. |
| `control` | `list()` | 传入 nlminb 的命名控制列表，如 iter.max、eval.max。 | A named list of nlminb control settings. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with report, opt, obj, est_std, vcov, pdHess, gradient, max_abs_gradient (also final_outer_mgc), convergence_code, bound_hit, year, age, length-bin metadata, growth_step, elapsed and multi-start diagnostics. start_diagnostics contains each start's PID, start/end timestamps, elapsed and CPU seconds, objective and convergence code. nstarts and workers record work and concurrency. The lowest finite objective is selected; check convergence separately. The DLL remains loaded so the returned obj can be evaluated in the same session.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl))
```

<a id="run_alscl"></a>
## `run_alscl()`

拟合年龄与长度联合 ALSCL / Fit joint age-length ALSCL

```r
run_alscl(data.CatL, data.wgt, data.mat, rec.age, nage, M, sel_L50,
    sel_L95, growth_step = 1, parameters = NULL, parameters.L = NULL,
    parameters.U = NULL, map = NULL, len_mid = NULL, len_border = NULL,
    len_lower = NULL, len_upper = NULL, output = FALSE, train_times = 1,
    ncores = 1, silent = FALSE, zero_action = c("missing", "error"),
    control = list(), nstarts = ncores)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `data.CatL` | `必填 / required` | 调查数量指数表：首列体长标签，其余列为时间。 | Data frame or matrix with length labels in its first column, followed by at least two numerically labelled time columns. |
| `data.wgt` | `必填 / required` | 平均单尾体重表，与调查表完全对齐。 | Weight-at-length, with identical labels, columns and dimensions. |
| `data.mat` | `必填 / required` | 成熟比例表，0–1，与调查表完全对齐。 | Maturity-at-length between 0 and 1, with the same layout. |
| `rec.age` | `必填 / required` | 补充年龄，单位为年。 | Recruitment age in years, greater than t0. |
| `nage` | `必填 / required` | 年龄组数量，包括最大年龄加组。 | Number of age classes including the plus group (at least two). |
| `M` | `必填 / required` | 每个模型时间步的自然死亡瞬时率。 | Nonnegative natural mortality per model time step. |
| `sel_L50` | `必填 / required` | 调查可捕性达到 50% 的体长。 | Lengths at 50 and 95 percent survey catchability. |
| `sel_L95` | `必填 / required` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Lengths at 50 and 95 percent survey catchability. |
| `growth_step` | `1` | 每个时间步包含的年数；季度为 0.25。 | Years per model step. ACL NULL uses rec.age when below 1, otherwise 1; ALSCL defaults to 1. Supply explicitly for quarterly data. |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | Named list of initial fixed-effect values. |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | Named lower and upper bounds. Bounds are aligned with free parameters after applying map. Initial free values are clamped to bounds. |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | Named lower and upper bounds. Bounds are aligned with free parameters after applying map. Initial free values are clamped to bounds. |
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | Named TMB factor mappings. factor(NA) fixes a parameter at its supplied value; NULL releases a default fixed parameter. Other defaults remain in force. |
| `len_mid` | `NULL` | 每个体长组的中心值向量，长度与数据行数相同。 | Finite length-bin midpoints; overrides automatic parsing. |
| `len_border` | `NULL` | 体长边界；拟合入口需要 nlen-1 个内部边界，模拟器使用 nlen+1 个完整边界。 | Interior boundaries, one fewer than the number of length bins. |
| `len_lower` | `NULL` | ALSCL 每组下界/上界向量，须与分组一致。 | Optional bounds consistent with len_border and len_mid. Infinite outer boundaries are allowed; the model treats the outer bins as tails. |
| `len_upper` | `NULL` | ALSCL 每组下界/上界向量，须与分组一致。 | Optional bounds consistent with len_border and len_mid. Infinite outer boundaries are allowed; the model treats the outer bins as tails. |
| `output` | `FALSE` | 是否导出拟合诊断与图；TRUE 使用 output 目录。 | Save diagnostic tables and plots below output/ when TRUE. |
| `train_times` | `1` | 每个初始点连续优化的次数。 | Number of successive optimization passes per start; not independent random starts. |
| `ncores` | `1` | 独立初始点拟合的最大并行进程数，实际不超过 nstarts；不是单个 TMB 拟合的线程数。 | Maximum socket workers across independent starts, capped at nstarts; not threads within one start. |
| `nstarts` | `ncores` | 初始点总数；默认等于 ncores 以兼容旧用法。固定此值可公平比较不同并行数；后续起点仅扰动自由参数。 | Total independent starts; defaults to ncores for compatibility. Fix this value for equal-work comparisons. Later starts deterministically jitter free parameters only. |
| `silent` | `FALSE` | TRUE 减少拟合进度输出。 | Suppress fitting progress messages when TRUE. |
| `zero_action` | `c("missing", "error")` | missing 排除调查零值；error 报错。NA 始终为缺失。 | Exclude zero observations as missing (historical behavior), or reject them with "error". NA is always missing; negative/infinite values fail. |
| `control` | `list()` | 传入 nlminb 的命名控制列表，如 iter.max、eval.max。 | A named list of nlminb control settings. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with report, opt, obj, est_std, vcov, pdHess, gradient, max_abs_gradient (also final_outer_mgc), convergence_code, bound_hit, year, age, length-bin metadata, growth_step, elapsed and multi-start diagnostics. start_diagnostics contains each start's PID, start/end timestamps, elapsed and CPU seconds, objective and convergence code. nstarts and workers record work and concurrency. The lowest finite objective is selected; check convergence separately. The DLL remains loaded so the returned obj can be evaluated in the same session.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl))
```

<a id="sim_acl"></a>
## `sim_acl()`

批量拟合 ACL 模拟重复 / Fit simulated replicates with ACL

```r
sim_acl(iter_range = 4:100, sim_data_path = ".", output_dir = ".",
    parameters = NULL, parameters.L = NULL, parameters.U = NULL,
    map = NULL, M = 0.2, ncores = 1, train_times = 1, control = list())
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `iter_range` | `4:100` | 重复编号；模拟时也用作随机种子。 | A numeric vector specifying the range of iterations to run the stock assessment model for (default is 4:100). |
| `sim_data_path` | `"."` | 包含 sim_rep 文件的目录。 | A character string specifying the path to the folder containing the simulation data files (default is the current working directory). |
| `output_dir` | `"."` | 输出目录；必要时自动建立。 | Character, the directory where simulation results will be saved.(default is the current working directory). |
| `parameters` | `NULL` | 命名参数初值列表，采用优化器尺度。 | A list containing the custom initial values for the parameters (default is NULL). |
| `parameters.L` | `NULL` | 命名下界列表，采用优化器尺度。 | A list containing the custom lower bounds for the parameters (default is NULL). |
| `parameters.U` | `NULL` | 命名上界列表，采用优化器尺度。 | A list containing the custom upper bounds for the parameters (default is NULL). |
| `map` | `NULL` | TMB 映射；factor(NA) 固定，NULL 可释放默认固定项。 | A list containing the custom values for the map elements (default is NULL). |
| `M` | `0.2` | 每个模型时间步的自然死亡瞬时率。 | Numeric, natural mortality (default: 0.2) |
| `ncores` | `1` | 并行数量；具体含义依函数而异，见英文说明。 | Integer. Number of CPU cores to use. 1 = sequential (default). Uses socket workers on all platforms. Use parallel::detectCores() to see available cores. |
| `train_times` | `1` | 每个初始点连续优化的次数。 | Number of successive optimization passes per start; not independent random starts. |
| `control` | `list()` | 传入 nlminb 的命名控制列表，如 iter.max、eval.max。 | Named nlminb control settings. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list containing the results of the stock assessment model.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
do.call(sim_acl, c(list(iter_range = 4:5,
  sim_data_path = "simulation_examples/flatfish",
  output_dir = "simulation_examples/fit", M = pa$M), x$fit_config$acl))
```

<a id="sim_cal"></a>
## `sim_cal()`

从参数计算生物矩阵 / Calculate biological arrays

```r
sim_cal(params)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `params` | `必填 / required` | initialize_params() 返回的完整模拟参数列表。 | A list of parameters from initialize_params(). Must contain model_type ("age_based" or "length_based") to determine which transition matrix to build. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list containing: model_typeCharacter: "age_based" or "length_based" LatALength-at-age vector W_at_lenWeight-at-length vector matMaturity-at-length vector q_survSurvey catchability-at-length vector selSelectivity vector (at-age for age_based, at-length for length_based) plaAge-length transition matrix (always built for initial conditions) GijGrowth transition matrix (length_based only) nlen, len_lower, len_upperLength bin dimensions

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
bio <- sim_cal(pa)
```

<a id="sim_data"></a>
## `sim_data()`

运行完整随机种群模拟 / Simulate full population dynamics

```r
sim_data(bio_vars, params, sim_year = NULL, output_dir = ".",
    iter_range = 4:100, return_iter = NULL)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `bio_vars` | `必填 / required` | sim_cal() 返回的生物矩阵及生成模型类型。 | A list from sim_cal() containing biological variables. |
| `params` | `必填 / required` | initialize_params() 返回的完整模拟参数列表。 | A list from initialize_params() containing parameters. |
| `sim_year` | `NULL` | 模拟总年数；NULL 使用 params$nyear。 | Total simulation years; NULL uses params$nyear. |
| `output_dir` | `"."` | 输出目录；必要时自动建立。 | Directory for saving results (default: current directory). |
| `iter_range` | `4:100` | 重复编号；模拟时也用作随机种子。 | Integer vector, which iterations (seeds) to run (default: 4:100). |
| `return_iter` | `NULL` | 选择在内存中返回的模拟重复编号。 | Integer or NULL, which iteration to return in memory. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list containing simulated fishery data for the observation window (post-burn-in), including SN_at_len, N_at_len, N_at_age, SSB, etc.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
sm <- sim_data(bio, pa, iter_range = 4:5, return_iter = 4,
               output_dir = "simulation_examples/flatfish")
```

<a id="simulate_example_data"></a>
## `simulate_example_data()`

快速年度教学数据生成 / Generate simple annual teaching data

```r
simulate_example_data(years = 2000:2025, bin_breaks = c(0, 20, 25, 30, 35,
    40, 45, 50, 55, 60), Linf = 55, vbk = 0.45, t0 = -0.1, M = 0.8,
    nage = 7, L50_sel = 28, L95_sel = 36, L50_mat = 33, L95_mat = 40,
    wgt_a = 2.2e-06, wgt_b = 3.2, mean_F = 0.3, cv_catch = 0.3,
    rec_sigma = 0.5, seed = NULL, save_csv = FALSE, output_dir = ".")
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `years` | `2000:2025` | 年份向量；快速模拟器要求连续年度，比较图用它筛选时期。 | Numeric vector of years. Default 2000:2025. |
| `bin_breaks` | `c(0, 20, 25, 30, 35, 40, 45, 50, 55, 60)` | 递增完整边界向量，相邻边界定义一个长度组。 | Numeric vector of bin boundary points. Defines n+1 boundaries for n bins. Default c(0, 20, 25, 30, 35, 40, 45, 50, 55, 60) giving 9 bins: 0-20, 20-25, 25-30, ..., 55-60. Works with any bin width (1mm, 2mm, 3mm, 5mm, etc.) and unequal widths. |
| `Linf` | `55` | VB 渐近体长，与输入体长单位一致。 | Numeric. VB asymptotic length. Default 55. |
| `vbk` | `0.45` | VB 生长系数，以年为时间单位。 | Numeric. VB growth coefficient. Default 0.45. |
| `t0` | `-0.1` | VB 理论零体长年龄，单位为年。 | Numeric. VB theoretical age at length 0. Default -0.1. |
| `M` | `0.8` | 每个模型时间步的自然死亡瞬时率。 | Numeric. Natural mortality. Default 0.8. |
| `nage` | `7` | 年龄组数量，包括最大年龄加组。 | Integer. Maximum age (plus group). Default 7. |
| `L50_sel` | `28` | 调查可捕性达到 50% 的体长。 | Numeric. Length at 50 percent selectivity. Default 28. |
| `L95_sel` | `36` | 调查可捕性达到 95% 的体长，必须大于 L50。 | Numeric. Length at 95 percent selectivity. Default 36. |
| `L50_mat` | `33` | 成熟比例为 0.5 的体长。 | Numeric. Length at 50 percent maturity. Default 33. |
| `L95_mat` | `40` | 成熟比例为 0.95 的体长，须大于 L50。 | Numeric. Length at 95 percent maturity. Default 40. |
| `wgt_a` | `2.2e-06` | 长度重量关系 W=a*L^b 的系数，需匹配单位。 | Numeric. Weight-length parameter a in W = a * L^b. Default 2.2e-6. |
| `wgt_b` | `3.2` | 长度重量关系 W=a*L^b 的指数。 | Numeric. Weight-length parameter b in W = a * L^b. Default 3.2. |
| `mean_F` | `0.3` | 每步平均捕捞死亡率。 | Numeric. Mean fishing mortality. Default 0.3. |
| `cv_catch` | `0.3` | 简易模拟调查误差的 CV，而非对数 SD。 | Numeric. CV of observation noise on catches. Default 0.3. |
| `rec_sigma` | `0.5` | 简易模拟对数补充扰动 SD。 | Numeric. SD of log-recruitment variation. Default 0.5. |
| `seed` | `NULL` | 随机种子；设定后可重复生成。 | Integer or NULL. Random seed for reproducibility. Default NULL (different result each run). Set e.g. seed = 42 for reproducible output. |
| `save_csv` | `FALSE` | TRUE 将生成的三张输入表写为 CSV。 | Logical. Save CSV files to working directory. Default FALSE. |
| `output_dir` | `"."` | 输出目录；必要时自动建立。 | Character. Directory for CSV output. Default ".". |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 A list with components: data.CatLCatch-at-length data.frame (LengthBin x Years) data.wgtWeight-at-length data.frame (same structure) data.matMaturity-at-length data.frame (same structure) true_paramsList of true parameter values used in simulation

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
simulate_example_data(years = 2000:2019, seed = 42,
                      save_csv = TRUE, output_dir = "quick_simulation")
```

<a id="VB_func"></a>
## `VB_func()`

计算 VB 平均体长 / Calculate VB mean length

```r
VB_func(Linf, k, t0, age)
```

| 参数 / Argument | 默认值 / Default | 中文说明 | English |
|---|---|---|---|
| `Linf` | `必填 / required` | VB 渐近体长，与输入体长单位一致。 | The asymptotic length. |
| `k` | `必填 / required` | VB 生长系数，以年为时间单位。 | The growth coefficient. |
| `t0` | `必填 / required` | VB 理论零体长年龄，单位为年。 | The theoretical age at zero length. |
| `age` | `必填 / required` | 用于计算生长的年龄向量，单位为年。 | Age of the fish. |

**返回 / Returns:** 返回结构及字段如下；可用 str() 检查。 Length at age.

**示例 / Example:**

```r
# 调用示例；变量按上文创建 / Use objects defined above
VB_func(Linf = 60, k = 0.2, t0 = 1/60, age = 1:15)
```

<a id="estimation_parameters"></a>
## 初值列表内部可以改什么 / Inside estimation parameter lists

`parameters`、`parameters.L`、`parameters.U` 的键名如下。表中变换指从优化参数到自然尺度；初值和界限都要在优化尺度填写。`create_parameters()` 返回的三个命名列表就是传给拟合入口的同名参数。

The keys below apply to starting values and bounds. Transformations map optimizer values to the natural scale; enter starts and bounds on the optimizer scale. The constructor's three named components match the corresponding fit arguments.

| ACL 键名 / Key | ALSCL 键名 / Key | 变换 / Transform | 含义 / Meaning |
|---|---|---|---|
| `log_init_Z` | `log_init_Z` | exp | 初始平衡总死亡率 / Initial equilibrium total mortality |
| `log_std_log_N0` | `log_sigma_log_N0` | exp | 初始年龄尾数偏差 SD / Initial log-number deviation SD |
| `mean_log_R` | `mean_log_R` | exp | 对数补充过程的位置参数 / Location of the log recruitment process |
| `log_std_log_R` | `log_sigma_log_R` | exp | 对数补充过程 SD 缩放 / Log recruitment process SD scale |
| `logit_log_R` | `logit_log_R` | plogis | 补充 AR(1) 相关，限定 (0,1) / Recruitment correlation restricted to (0,1) |
| `mean_log_F` | `mean_log_F` | exp | 对数 F 过程的位置参数 / Location of the log F process |
| `log_std_log_F` | `log_sigma_log_F` | exp | F 偏差 SD 缩放 / F deviation SD scale |
| `logit_log_F_y` | `logit_log_F_y` | plogis | F 的时间相关 / Temporal F correlation |
| `logit_log_F_a` | `logit_log_F_l` | plogis | F 的年龄/体长相关 / Age/length F correlation |
| `log_vbk` | `log_vbk` | exp | VB k，单位为每年 / Annual VB k |
| `log_Linf` | `log_Linf` | exp | 渐近体长 / Asymptotic length |
| `t0` | `log_t0` | ACL 原值；ALSCL exp / Identity vs exp | VB t0；ALSCL 本版限定正值 / ALSCL restricts t0 to positive values |
| `log_cv_len` | `log_cv_len` | exp | 年龄体长 CV / Length-at-age CV |
| — | `log_cv_grow` | exp | 增长增量 CV / Growth-increment CV |
| `log_std_index` | `log_sigma_index` | exp | 对数调查观测 SD / Log survey observation SD |

`exp(mean_log_R)` 和 `exp(mean_log_F)` 是对数位置的指数转换；有随机偏差时不能直接称作无条件算术平均数。两个模板当前都不支持负的 AR(1) 相关。SD 在相关过程中的具体尺度见指南原理部分。

Exponentiating a log location does not give the unconditional arithmetic mean when random deviations are present. Both templates currently restrict AR(1) correlations to positive values. See the theory section for SD conventions in correlated processes.

```r
# 查看某物种的全部初值和上下界 / Inspect all starts and bounds for a preset
p <- create_parameters(model_type = "alscl", species = "flatfish")
keys <- names(p$parameters)
parameter_table <- data.frame(
  Parameter = keys,
  Start = unlist(p$parameters[keys]),
  Lower = unlist(p$parameters.L[keys]),
  Upper = unlist(p$parameters.U[keys]))
parameter_table
# 改一个初值及上界：此处并未自动固定该参数 / Change a start and upper bound
p$parameters$log_sigma_index <- log(0.15)
p$parameters.U$log_sigma_index <- log(0.5)
# 要固定该 SD，还需显式 map；固定值取自 parameters / Explicit map fixes the value
fixed_index <- list(log_sigma_index = factor(NA))
```

拟合会按数据维度建立 `dev_log_R`、`dev_log_F`、`dev_log_N0` 随机效应。它们不是本表的普通标量初值；不要把生成器的 `std_logR`、`F_mean` 等直接作为估计列表键名。

Fits construct random effects `dev_log_R`, `dev_log_F` and `dev_log_N0` from data dimensions. These are not ordinary scalar entries above. Generator names such as `std_logR` and `F_mean` are not estimator keys.
