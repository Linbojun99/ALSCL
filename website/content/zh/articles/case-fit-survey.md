# 如何拟合自己的实测调查数据？

## 要解决的问题

你已有按体长和时间记录的调查数量，以及平均个体重量和成熟比例，希望用 ALSCL 拟合。目标是建立可追溯的资源评估流程。本案例不提供实测数据，也不替你的研究种群指定生物假设。

## 准备三张表

复制 [Excel 模板](../assets/manual/data/ALSCL_Data_Entry_Example.xlsx)，将模拟观测替换为自己的数据。保留三张工作表：`data.CatL`、`data.wgt` 和 `data.mat`。第一列为 `LengthBin`，其余表头为递增数值时间；体长组、时间及顺序必须一致。

记录物种、调查设计、指数标准化方式、体长与重量单位、体长边界和零值含义。重量填写平均个体重量，成熟比例取值为 [0,1]。[数据准备](data-preparation.html) 说明缺失时间点与不等宽体长组的处理。

## 工作簿应是什么样子

以下截图来自配套的 **模拟 YTF 工作簿**，只说明填写格式。将数值替换为自己的观测，并另存原模板。

![调查数量输入表](../assets/manual/figures/excel_CatL.png)

填写口径一致的数量指数。空白/NA 表示缺失，不能将缺失观测改填为零。

![平均个体重量输入表](../assets/manual/figures/excel_wgt.png)

每个单元格填写平均个体重量，不是总捕获重量。记录重量单位及所用体长—体重关系。

![成熟比例输入表](../assets/manual/figures/excel_mat.png)

填写零到一之间的比例。三张表的体长组、时间及顺序必须完全一致。

## 导入与验证

先安装一次 `readxl`，再从下载的仓库根目录运行，以便载入指南辅助函数。将工作簿路径改成自己的文件。

```r
library(ALSCL)
source("scripts/real_data_workflow.R", encoding = "UTF-8")
my_inputs <- read_survey_excel("my_survey.xlsx")
step_years <- 1
validate_survey_tables(my_inputs, growth_step = step_years, zero_action = "error")
```

只有在时间序列为规则季度间隔且生物参数单位一致时，才使用 `step_years <- 0.25`。验证默认在遇到零值时停止。软件包默认排除零观测的行为不一定适合真实零捕获；包内未实现零膨胀观测模型。

## 先完整运行配套工作簿

这是配置完整的可选练习，使用模板内模拟 YTF 数据及 **与其匹配的生物设置**，便于在替换观测前检查安装、导入与编译。它并不是对你的研究种群进行评估。

```r
library(ALSCL)
source("scripts/real_data_workflow.R", encoding = "UTF-8")
data("YTF_example")
demo_inputs <- read_survey_excel("docs/data/ALSCL_Data_Entry_Example.xlsx")
demo_config <- c(YTF_example$fit_args, YTF_example$fit_config$alscl,
                 list(zero_action = "error", train_times = 2, silent = TRUE))
demo_fit <- fit_survey_tables(demo_inputs, demo_config, model_type = "alscl")
diagnose_model(demo_inputs$data.CatL, demo_fit)
plot_CatL(demo_fit, type = "year", exp_transform = FALSE, facet_ncol = 4)
```

![模拟 YTF 示意：各时期的观测与拟合调查指数](../assets/manual/figures/ytf/CatL_year_FALSE.png)

在图示的对数尺度上比较观测与拟合分布。某一体长范围内持续低估或高估，提示值得检查的结构性问题。上图对应模板练习，并不是你自己的工作簿结果。

## 指定评估假设

下面是需要填写的配置模板。拟合前必须根据研究种群的证据替换所有 `NA`。检查语句会主动阻止不完整配置继续运行。

```r
config <- list(
  rec.age = NA_real_, nage = NA_integer_, M = NA_real_,
  sel_L50 = NA_real_, sel_L95 = NA_real_,
  growth_step = step_years, zero_action = "error",
  train_times = 2, silent = TRUE)
required_biology <- c("rec.age", "nage", "M", "sel_L50", "sel_L95")
stopifnot(all(is.finite(unlist(config[required_biology]))))
```

这只是配置起点，并非完整生物设定。还应审查实际体长边界、生长假设、初值、边界与固定参数，将相应的命名拟合参数加入 `config`。参见 [拟合控制](fitting-controls.html)、[估计参数](estimation-parameters.html) 与 [run_alscl()](../reference/run_alscl.html)。不能因为表格维度相似就直接套用 YTF 的生长或可捕性设置。死亡率须与模型时间步一致。

### 拟合实测观测前审查这些输入

|输入|需要说明的依据|
|---|---|
|`rec.age`、`nage`|补充年龄与模型年龄组数量|
|`M`|每个模型时间步的自然死亡率、来源与不确定性|
|`sel_L50`、`sel_L95`|固定的调查可捕性位置，不是估计的渔业选择性|
|`growth_step`|每步年数，例如 1 或 0.25，须与时间表头匹配|
|`len_mid`、`len_border`|实际体长中点与内部边界；不等宽组需明确含义|
|`parameters`、`parameters.L`、`parameters.U`|优化尺度上的初值与上下界|
|`map`|固定了哪些参数及其理由；区间以固定设置为条件|

ALSCL 的生长初值使用 `log_Linf`、`log_vbk` 与 `log_t0`。此参数化限制 `t0` 为正，不能直接对负的外部估计取对数。应审查模型是否适用，而不是悄悄改符号。提供初值不等于固定参数，固定还需相应的 map 项，例如 `factor(NA)`。

`train_times` 是连续优化轮数，不是独立随机起点。不同初值需另行检查，参见 [拟合控制](fitting-controls.html)。

## 拟合并检查

补全并审查 `config` 后进行拟合，在解释种群估计之前检查数值诊断。

```r
my_fit <- fit_survey_tables(my_inputs, config, model_type = "alscl")
checks <- diagnose_model(my_inputs$data.CatL, my_fit)
checks
plot_residuals(my_fit, type = "year")
plot_biomass(my_fit)
```

分别查看优化状态、梯度、Hessian、参数边界与残差结构。拟合完成不代表全部参数可识别。如果诊断或敏感性检查较差，应重新审查假设；[诊断说明](diagnostics.html) 介绍相关工具。

## 解释诊断图与种群结果图

以下两图使用仓库 YTF 拟合说明读图方法。对 `my_fit` 运行命令，才能得到自己观测的对应图。

![模拟 YTF 示意：逐年残差结构](../assets/manual/figures/ytf/residuals_year.png)

观察残差是否持续同向以及离散程度是否变化。这些是原始对数残差；发现结构时应检查指数标准化、可捕性假设和过程结构，而不是自动删去观测。

![模拟 YTF 示意：生物量轨迹与条件不确定性](../assets/manual/figures/ytf/plot_biomass_B.png)

结合输入重量单位与固定假设解释轨迹。图示区间不包含全部外部生物输入的不确定性；看起来精确不能替代敏感性分析。

### 导出自己的结果

```r
dir.create("my_assessment", showWarnings = FALSE)
write.csv(checks, "my_assessment/diagnostics.csv", row.names = FALSE)
p_index <- plot_CatL(my_fit, type = "year", facet_ncol = 4)
p_residual <- plot_residuals(my_fit, type = "year")
p_biomass <- plot_biomass(my_fit)
ggplot2::ggsave("my_assessment/survey_fit.png", p_index, width = 10, height = 7, dpi = 300)
ggplot2::ggsave("my_assessment/residuals.png", p_residual, width = 10, height = 6, dpi = 300)
ggplot2::ggsave("my_assessment/biomass.png", p_biomass, width = 8, height = 5, dpi = 300)
```

### 遇到问题时

|现象|下一步处理|
|---|---|
|工作表或表头错误|恢复准确的工作表名、`LengthBin` 与数值时间表头|
|三张表不一致|统一体长组与时间列及其顺序|
|零观测被拒绝|先核实零值含义，不要替换成任意小常数|
|编译错误|按操作系统安装说明检查工具链与可写 TMB 缓存|
|Hessian 非正定或梯度较大|解释不确定性前检查初值、边界、固定假设与可识别性|
|不同 M 或 q 下结果差异大|报告有依据的外部值的敏感性，不要只选择最方便的结果|

## 保存评估记录

```r
saveRDS(list(inputs = my_inputs, config = config,
             diagnostics = checks, report = my_fit$report),
        "my_assessment_summary.rds")
writeLines(capture.output(sessionInfo()), "my_assessment_session.txt")
```

同时保留输入来源与软件包提交版本。保存的汇总结果可用于报告；需要当前 TMB 目标函数的操作，必须在当前 R 会话中重新拟合。处理实测数据前，可先运行配置完整的 [模拟 YTF 示例](first-model.html)。
