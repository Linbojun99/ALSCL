# 年度与季度案例

两案例重新使用当前包拟合。它们用于学习年度/季度工作流，不能替代原论文数据复现；单个随机重复也不能给两个模型排名。所有拟合图的通用操作仍以 `YTF_example` 为主线。


|设置| 黄尾鲽 Flatfish | 大眼金枪鱼 Tuna |
|---|---|---|
|生成动力学| age_based | length_based |
|时间步| 1 | 0.25 |
| 总年数 / 预热年数 · Total / burn-in | 100 / 80 | 100 / 95 |
|保留步数| 20 | 20 |
| rec.age / nage | 1 / 15 | 0.25 / 20 |
| Linf / annual k | 60 / 0.20 | 152 / 0.38 |
| t0, years | 1/60 | 1/152 |
| M / mean F, per step | 0.20 / 0.30 | 0.20 / 0.20 |
|成熟 L50/L95| 35 / 40 | 100 / 120 |
|调查 q L50/L95| 15 / 20 | 30 / 50 |
|长度中点| 6, 8, …, 50 | 12.5, 17.5, …, 122.5 |
|调查 log SD| 0.20 | 0.10 |

这里 flatfish 的调查 SD = 0.20，主线 `YTF_example` 则沿用 YTF 的 0.10，不能混称同一数据。tuna 的 M = 0.2 为每季度率。长度与体重单位必须与幂函数系数匹配；示例年份从 2000 开始仅用于坐标标签。


## 读取和重建

```r
# 读取已有案例
source("scripts/real_data_workflow.R", encoding="UTF-8")
flatfish_inputs <- read_survey_csv("docs/data/cases/flatfish")
tuna_inputs <- read_survey_csv("docs/data/cases/tuna")
validate_survey_tables(flatfish_inputs, growth_step=1)
validate_survey_tables(tuna_inputs, growth_step=.25)
pa <- readRDS("docs/data/cases/tuna/truth.rds")$parameters
str(pa)
# 重新产生数据
source("scripts/simulation_workflow.R", encoding="UTF-8")
simulations <- run_simulation_examples(out="my_simulations", fit_batch=FALSE)
# 重跑两个案例的拟合和七张图
source("scripts/case_studies.R", encoding="UTF-8")
case_status <- run_case_studies(out="docs")
```

最后一行读取 `out/data/cases` 的已保存真值并将输出写入同一目录下。要换目录，先复制案例输入及真值。完整拟合初值、map 和生物参数构造均在 [case_studies.R](https://github.com/Linbojun99/ALSCL/blob/main/scripts/case_studies.R)，数据下载见 [两个案例目录](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/cases/README.md)。


## 年度黄尾鲽

![模拟调查](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_data.png?raw=true)

颜色为 log10 调查数量指数。沿时间和体长延伸的条带可反映队列结构；它们仍同时受调查可捕性和观测误差影响。


![生物设定](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_biology.png?raw=true)

生长曲线决定年龄与平均长度的联系；成熟概率和调查 q 有不同的 L50/L95，不可相互替代。


![模拟真值与条件估计](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/flatfish_truth.png?raw=true)

黑虚线为生成真值，彩色线是固定部分生物参数的条件估计。比较每种量的趋势和偏移；同一时期曲线接近并不证明所有参数可识别。


## 季度大眼金枪鱼

![模拟调查](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_data.png?raw=true)

20 个保留步长覆盖 5 年，非 20 年。物种预设同时改变生成动力学、长度组、年龄网格及生物参数。


![生物设定](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_biology.png?raw=true)

![模拟真值与条件估计](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/cases/tuna_truth.png?raw=true)

tuna 由长度型过程生成；拟合仍要明确 `growth_step=.25`。案例固定了真实生长、CV、部分过程离散及相关参数，估计的区间不能覆盖这些假设本身的误差。


## 如何扩展为模拟实验

改变 `iter_range` 生成多个随机重复；系统改变观测误差、M、q、长度分箱或数据长度。保留所有失败拟合和诊断，再汇总偏差、RMSE、区间覆盖率及回溯表现。仅平均“成功重复”会形成选择偏差。


[本次四个拟合的数值状态](https://github.com/Linbojun99/ALSCL/blob/main/docs/results/case_convergence.csv)。
