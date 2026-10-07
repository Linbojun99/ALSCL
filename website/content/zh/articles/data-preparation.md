# 数据、Excel 与实测导入

下载 [ALSCL_Data_Entry_Example.xlsx](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx)。工作簿的 `Instructions` 是双语说明页；其余三张表已填入与 `YTF_example` 相同的数值。请另存副本再替换数据。


|位置|填写规则|
|---|---|
|工作表名|必须为 `data.CatL`、`data.wgt`、`data.mat`|
| A1 |`LengthBin`，不要在表头上方加标题|
| A2:A… |每行一个体长组；示例为中心值 6、8、…、50|
| B1、C1、… |递增数字年份，如 2000、2001；季度如 2000、2000.25|
|B2 等数值单元|对应这一体长组、这一期的值|
|三表顺序|体长组、时间、行列顺序必须完全相同|

### 调查数量指数 `data.CatL`

![调查数据填表示例](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_CatL.png?raw=true)

填数量或已标准化为可比尺度的数量指数；若使用 CPUE，各期努力量、单位和标准化方式必须一致。正数可为小数。未知调查单元留空或填 `NA`；不要用 0、破折号、`<1` 或“未测”代替缺失值。


### 平均个体重量 `data.wgt`

![平均体重填表示例](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_wgt.png?raw=true)

填每组每期的平均单尾重量，不是该组总重。选择同一重量单位并记录在说明页。没有逐年数据时，可以在有依据的条件下重复同一长度重量关系，但要在分析报告说明。不得留缺失值。


### 成熟比例 `data.mat`

![成熟比例填表示例](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_mat.png?raw=true)

填写 0–1，例如 60% 填 `0.6`，不要填 `60`。0 和 1 在成熟表中是合法值；不要将“调查零值”的规则误用于成熟表。三张表的数据区不要插入单位行、合计、合并单元格或解释文字。


Excel 使用科学记数显示小数，例如 `5.1361E-02` 等于 `0.051361`；显示精度不会改变保存的底层数值。

**体长边界与缺测期 / Bin boundaries and missing periods.** 本示例中心值是 6:50，每隔 2 cm；拟合时显式传入 22 个内部边界 7:49，首尾尾部按模型约定处理。自己的不等宽组应提供真实边界。区间标签如 `0-20`、`20-25` 要连续且不重叠；原标签若为 `0-20`、`21-25`，不能未经核实就消除间隙。缺整个年度时保留对应列并在调查表填 NA，重量和成熟表仍需完整；不要删掉这一列而假装时间间隔没变。


## 导入自己的实测数据

[real_data_workflow.R](https://github.com/Linbojun99/ALSCL/blob/main/scripts/real_data_workflow.R) 提供指南辅助函数，需 `source()`，不属于包的导出 API。Excel 读取采用 [readxl::read_excel](https://readxl.tidyverse.org/reference/read_excel.html)，保留原始年份表头并显式转换数值。


```r
library(ALSCL)
data("YTF_example")
install.packages("readxl") # 只需安装一次
source("scripts/real_data_workflow.R", encoding = "UTF-8")
# 先用附带工作簿练习
my_inputs <- read_survey_excel("docs/data/ALSCL_Data_Entry_Example.xlsx")
validate_survey_tables(my_inputs, growth_step = 1, zero_action = "error")
# CSV 也可以，目录内需有三张同名 CSV
csv_inputs <- read_survey_csv("docs/data")
# 教学工作簿使用其对应的生物设定
config <- c(YTF_example$fit_args, YTF_example$fit_config$alscl,
            list(zero_action = "error", train_times = 2, silent = TRUE))
# 真正重新拟合；需要编译环境
my_fit <- fit_survey_tables(my_inputs, config, model_type = "alscl")
```

换成自己的文件时，替换路径，并依据自己的物种设置 `rec.age`、`nage`、`M`、`sel_L50`、`sel_L95`、`growth_step`、边界及拟合初值/固定参数。不能因为输入表尺寸相同，就沿用 YTF 的生物参数或固定参数表。生长、长度重量关系、成熟信息和调查可捕性应有独立依据。


`validate_survey_tables()` 默认遇到调查零值报错，要求先确认其含义。包的拟合入口默认 `zero_action="missing"`，把零值从对数似然中排除。如果零是真实的无捕获结果，排除它会改变推断；本包没有提供零膨胀观测模型，不应随意加小常数伪装为正数。


### 1 原有 `example_data`

```r
data("example_data")
names(example_data)
head(example_data$data.CatL)
# 查看标签与零值，暂不盲目拟合
example_data$data.CatL[[1]]
sum(as.matrix(example_data$data.CatL[-1]) == 0, na.rm = TRUE)
```

这三张表从原 ALSCL 仓库原样保留，原说明称匿名调查数据，列为 2001–2021，9 个长度组。物种、单位、原始出处和零值含义没有得到核实，而且 `0-20`、`21-25` 等标签存在间隔。当前长度检查会拒绝这些不连续区间。读者应先补齐元数据和真实分组边界，不能把它标注成本文论文的已核实实测数据，也不能擅自改标签后给出正式评估。
