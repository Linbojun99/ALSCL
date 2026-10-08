# 如何模拟调查数据？

## 要解决的问题

你需要可控的调查观测来测试分析流程，或研究模型表现。本案例区分“快速检查输入表”和“保留生成真值的模拟实验”两条路径。两条路径均不提供实测渔业数据。

## 先明确模拟实验的问题

|问题|路径|应保留的内容|
|---|---|---|
|导入或绘图能否正常工作？|快速年度表|种子、设置及三张输入表|
|估计与生成真值有何差异？|完整操作模型|参数、生物数组、真值、观测及全部拟合结果|
|年度与季度设置有何差异？|显式指定时间步的完整预设|时间单位及所有改变的生物假设|

本案例配图来自仓库的完整 flatfish 预设、第 4 次重复，对应下文完整流程，不是快速输入表模拟器的结果。flatfish 在 80 年预热后保留 20 个年度观测；tuna 在 95 年预热后保留 20 个季度观测。

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

## 验证并展示模拟观测

转换器返回三张数据框，第一列为体长组标签。将其交给估计模型前先做验证。热图使用完整预设中的实际体长中点。

```r
source("scripts/real_data_workflow.R", encoding = "UTF-8")
validate_survey_tables(inputs, growth_step = pa$growth_step, zero_action = "error")
library(ggplot2)
years <- as.numeric(names(inputs$data.CatL)[-1])
observations <- data.frame(
  Year = rep(years, each = length(pa$len_mid)),
  Length = rep(pa$len_mid, length(years)),
  Index = as.vector(as.matrix(inputs$data.CatL[-1])))
p_data <- ggplot(observations, aes(Year, Length, fill = log10(Index))) +
  geom_tile() + scale_fill_gradient(low = "white", high = "steelblue") +
  labs(fill = "log10 index") + theme_bw()
p_data
ggsave("my_full_simulation/survey.png", p_data, width = 10, height = 4, dpi = 160)
```

![完整 flatfish 模拟：体长与年份上的调查数量](../assets/manual/figures/cases/flatfish_data.png)

**读图：**颜色为调查指数的 log10。体长组间移动的条带反映队列推进，同时受可捕性和噪声影响。这是观测分布，不是总丰度；调查可捕性变化可使其改变，而种群规模不一定按同样幅度改变。

![生成模型的生长曲线、成熟比例与调查可捕性](../assets/manual/figures/cases/flatfish_biology.png)

上图将年龄与平均体长联系起来；下图的成熟比例与调查可捕性承担不同作用，不能混用各自的 L50/L95。这些生物值属于 flatfish 预设，不应直接作为未知种群的假设。

## 对同一次重复拟合两种模型

作为条件估计演示，按照仓库案例脚本固定已知生长及部分过程参数。ACL 与 ALSCL 使用不同的参数名与变换。下面代码使用刚生成的 `pa` 与 `inputs`。

```r
pA <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk), t0 = pa$t0,
           log_cv_len = log(pa$cv_L), log_std_log_F = log(pa$F_sd),
           log_std_log_N0 = log(pa$std_logN0), logit_log_R = qlogis(pa$R_ar))
pB <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk), log_t0 = log(pa$t0),
           log_cv_len = log(pa$cv_L), log_cv_grow = log(pa$cv_inc),
           log_sigma_log_F = log(pa$F_sd), log_sigma_log_N0 = log(pa$std_logN0),
           logit_log_R = qlogis(pa$R_ar))
mA <- lapply(pA, function(value) factor(NA))
mB <- lapply(pB, function(value) factor(NA))
mB$logit_log_F_l <- mB$logit_log_F_y <- factor(NA)
common <- c(inputs, list(
  rec.age = pa$rec.age, nage = pa$nage, M = pa$M,
  sel_L50 = pa$q_surv_L50, sel_L95 = pa$q_surv_L95,
  len_mid = pa$len_mid, len_border = pa$len_border[-c(1, length(pa$len_border))],
  growth_step = pa$growth_step, train_times = 2, ncores = 1, silent = TRUE))
a <- do.call(run_acl, c(common, list(parameters = pA, map = mA)))
b <- do.call(run_alscl, c(common, list(parameters = pB, map = mB)))
checks <- list(ACL = diagnose_model(inputs$data.CatL, a),
               ALSCL = diagnose_model(inputs$data.CatL, b))
```

## 将估计与生成真值比较

模拟器的总生物量名为 `TB`，拟合报告中为 `B`。在计算差值或比值前，先对齐相同时间段和指标。

```r
truth <- data.frame(Year = rep(years, 3),
  Quantity = rep(c("B", "SSB", "Rec"), each = length(years)),
  Truth = c(sm$TB, sm$SSB, sm$Rec))
estimates <- rbind(
  data.frame(Year = rep(years, 3), Quantity = truth$Quantity,
    Value = c(a$report$B, a$report$SSB, a$report$Rec), Model = "ACL"),
  data.frame(Year = rep(years, 3), Quantity = truth$Quantity,
    Value = c(b$report$B, b$report$SSB, b$report$Rec), Model = "ALSCL"))
p_truth <- ggplot(estimates, aes(Year, Value, color = Model)) + geom_line() +
  geom_line(data = truth, aes(Year, Truth), inherit.aes = FALSE, linetype = 2) +
  facet_wrap(~ Quantity, scales = "free_y", ncol = 1) + theme_bw()
p_truth
ggsave("my_full_simulation/truth_comparison.png", p_truth,
       width = 10, height = 6, dpi = 160)
saveRDS(list(parameters = pa, truth = sm, inputs = inputs, checks = checks,
             ACL_report = a$report, ALSCL_report = b$report),
        "my_full_simulation/assessment_summary.rds")
```

![flatfish 第 4 次重复：生成真值与 ACL/ALSCL 条件估计](../assets/manual/figures/cases/flatfish_truth.png)

虚线为生成真值，彩色线为仓库案例中的条件估计。在同一面板内比较方向与幅度；各指标的纵轴自由缩放，不能跨面板比较视觉振幅。该例采用年龄结构生成模型，且只展示一次重复，不能据此判断模型普遍优劣。

## 扩展到季度或更多重复

季度模拟从 `initialize_params(species = "tuna")` 开始，重新生成生物数组、观测与匹配的拟合设置，不要仅修改年度表头。多次模拟应保留种子、真值与失败记录；`sim_acl()` 只自动拟合 ACL，ALSCL 需对各次输入循环调用 `run_alscl()`。

|现象|优先检查|
|---|---|
|时间表头与 `growth_step` 不一致|转换时使用生成该重复的同一参数列表|
|轨迹不合理|单位、预热时长、生长、补充年龄与调查可捕性|
|拟合失败|保留错误和种子，先检查数值诊断再汇总成功拟合|
|无法读取文件|`sim_rep*` 用 `load()`；`saveRDS()` 保存的文件用 `readRDS()`|

## 应检查与报告什么

拟合前确认三表对齐、时间单位、保留时长、生成生物参数与观测误差。两个重复只能演示文件处理，不能确定模型表现。正式模拟研究应有目的地改变目标假设，并在汇总偏差、RMSE 与区间覆盖率时同时保留拟合失败记录。

继续阅读 [模拟函数说明](simulation.html)、[年度与季度设置](case-studies.html)，或 [数据表拟合与生物参数设置](case-fit-survey.html)。
