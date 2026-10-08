# 如何拟合自己的实测调查数据？

## 要解决的问题

你已有按体长和时间记录的调查数量，以及平均个体重量和成熟比例，希望用 ALSCL 拟合。目标是建立可追溯的资源评估流程。本案例不提供实测数据，也不替你的研究种群指定生物假设。

## 准备三张表

复制 [Excel 模板](../assets/manual/data/ALSCL_Data_Entry_Example.xlsx)，将模拟观测替换为自己的数据。保留三张工作表：`data.CatL`、`data.wgt` 和 `data.mat`。第一列为 `LengthBin`，其余表头为递增数值时间；体长组、时间及顺序必须一致。

记录物种、调查设计、指数标准化方式、体长与重量单位、体长边界和零值含义。重量填写平均个体重量，成熟比例取值为 [0,1]。[数据准备](data-preparation.html) 说明缺失时间点与不等宽体长组的处理。

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

## 保存评估记录

```r
saveRDS(list(inputs = my_inputs, config = config,
             diagnostics = checks, report = my_fit$report),
        "my_assessment_summary.rds")
writeLines(capture.output(sessionInfo()), "my_assessment_session.txt")
```

同时保留输入来源与软件包提交版本。保存的汇总结果可用于报告；需要当前 TMB 目标函数的操作，必须在当前 R 会话中重新拟合。处理实测数据前，可先运行配置完整的 [模拟 YTF 示例](first-model.html)。
