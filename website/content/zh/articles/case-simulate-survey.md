# 如何模拟调查数据？

## 要解决的问题

你需要可控的调查观测来测试分析流程，或研究模型表现。本案例区分“快速检查输入表”和“保留生成真值的模拟实验”两条路径。两条路径均不提供实测渔业数据。

## 准备会话

安装 ALSCL，下载或克隆仓库，并从仓库根目录运行 R。下面的辅助脚本随仓库提供，不属于已安装软件包的公开函数。

```r
library(ALSCL)
source("scripts/simulation_workflow.R", encoding = "UTF-8")
```

## 路径一：生成年度输入表

用这条路径测试导入、验证与绘图。若要保留以前的文件，请使用新的输出目录。

```r
quick <- simulate_example_data(
  years = 2000:2019, seed = 42, bin_breaks = seq(5, 51, 2),
  Linf = 60, vbk = 0.2, t0 = 1/60, M = 0.2, nage = 15,
  L50_sel = 15, L95_sel = 20, L50_mat = 35, L95_mat = 40,
  wgt_a = exp(-12), wgt_b = 3, mean_F = 0.3,
  cv_catch = 0.2, rec_sigma = 0.3,
  save_csv = TRUE, output_dir = "my_annual_simulation")
str(quick)
```

检查返回的表和导出的 CSV 文件。这个简化年度模拟器与完整操作模型不同；把年份标签改成季度，并不会产生季度动态。

## 路径二：保留真值与模拟重复

当问题涉及偏差、不确定性或操作模型与估计模型的差异时，使用这条路径。

```r
pa <- initialize_params(species = "flatfish", observation_error = "independent")
bio <- sim_cal(pa)
sm <- sim_data(bio, pa, iter_range = 4:5, return_iter = 4,
               output_dir = "my_full_simulation")
inputs <- simulation_to_tables(sm, pa)
saveRDS(list(parameters = pa, truth = sm, inputs = inputs),
        "my_full_simulation/truth_and_inputs.rds")
str(inputs)
```

`iter_range` 指定重复编号与随机种子；`return_iter` 选择返回内存的重复。`simulation_to_tables()` 转置调查数量，并根据 `growth_step` 构造时间表头。保存的 `sim_rep*` 文件用 `load()` 读取，不是 `readRDS()`。

## 应检查与报告什么

拟合前确认三表对齐、时间单位、保留时长、生成生物参数与观测误差。两个重复只能演示文件处理，不能确定模型表现。正式模拟研究应有目的地改变目标假设，并在汇总偏差、RMSE 与区间覆盖率时同时保留拟合失败记录。

继续阅读 [模拟函数说明](simulation.html)、[年度与季度设置](case-studies.html)，或 [数据表拟合与生物参数设置](case-fit-survey.html)。
