> 用于渔业资源评估的年龄与体长结构统计捕获体长模型

**ALSCL** 是基于 [Template Model Builder（TMB）](https://github.com/kaskr/adcomp) 的 R 包，利用独立于渔业的调查体长数据拟合种群模型。包内提供年龄结构模型 **ACL** 和年龄—体长联合结构模型 **ALSCL**，以及模拟、诊断、模型比较和绘图工具。

可以从 [完整拟合示例](articles/first-model.html) 开始，从 [首页目录](#contents) 查阅基本功能，阅读 [专题案例](articles/index.html)，或通过 [函数参考](reference/index.html) 查找具体用法。模型框架基于 [Zhang & Cadigan（2022）](https://doi.org/10.1111/faf.12673)。

<h2 id="contents">目录</h2>

- [安装](#installation)
- [两种模型结构](#models)
- [基本使用](#basic-use)
- [检查与解释拟合](#check-fit)
<!-- ARTICLE_CONTENTS -->
- [获取帮助](#help)
- [引用](#citation)

<h2 id="installation">安装</h2>

从 GitHub 安装当前开发版本：

```r
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
```

拟合时会编译包内的 TMB 模型，因此需要与 R 匹配的 C++ 工具链。各操作系统的说明见 [安装与编译](articles/getting-started.html)。

<h2 id="models">两种模型结构</h2>

| | ACL | ALSCL |
|---|---|---|
| 种群状态 | 年龄 × 时间的尾数 | 体长 × 年龄 × 时间的尾数 |
| 生长表达 | 年龄—体长概率矩阵 | 生长转移矩阵 |
| 捕捞死亡率 | 按年龄估计 | 按体长估计；年龄别 F 由模型推导 |
| 拟合函数 | `run_acl()` | `run_alscl()` |

两种模型都使用三张对齐的输入表：调查体长别数量、体长别平均个体重量和体长别成熟比例。自然死亡率及调查可捕性设置由使用者提供。这些假设影响估计的尺度与解释，详见 [模型原理](articles/model-theory.html)。

<h2 id="basic-use">基本使用</h2>

内置 `YTF_example` 包含黄尾鲽风格的模拟数据、生物参数设置和已知模拟真值。它是教学数据，不是论文的实测样本。

```r
data("YTF_example")
x <- YTF_example
inputs <- x[c("data.CatL", "data.wgt", "data.mat")]

# 使用相同观测与生物设置拟合两个模型
a <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,
                       list(train_times = 2, silent = TRUE)))
b <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
                         list(train_times = 2, silent = TRUE)))

plot_compare_ts(a, b, quantities = c("B", "SSB", "Rec"), se = TRUE)
```

![模拟 YTF 示例中 ACL 与 ALSCL 的生物量、产卵生物量和补充量估计](assets/manual/figures/ytf/compare_ts.png)

本示例固定了部分生物与过程参数，区间以这些设置为条件。[完整示例](articles/first-model.html) 解释输入、拟合结果与数据来源。

<h2 id="check-fit">检查与解释拟合</h2>

```r
diagnose_model(x$data.CatL, b)
plot_residuals(b, type = "year")
comparison <- compare_models(a, b, x$data.CatL)
comparison$summary
```

分别检查优化状态、最大绝对梯度、Hessian 及参数是否触及边界，再检查残差结构和对固定假设的敏感性。优化器显示成功，不能单独证明推断可靠。详见 [诊断](articles/diagnostics.html) 与 [模型比较](articles/model-comparison.html)。

<h2 id="learn-more">基本功能</h2>

<!-- ARTICLE_GUIDE -->

<h2 id="help">获取帮助</h2>

通过 [GitHub Issues](https://github.com/Linbojun99/ALSCL/issues) 报告问题或提出功能需求。建议附上最小可复现示例、包版本、输入维度及 `sessionInfo()` 输出。常见数据和编译问题见 [故障处理](articles/reproducibility.html)。

<h2 id="citation">引用</h2>

Zhang, F. & Cadigan, N. G. (2022). An age- and length-structured statistical catch-at-length model for hard-to-age fisheries stocks. *Fish and Fisheries*, **23**(5), 1121–1135. [doi:10.1111/faf.12673](https://doi.org/10.1111/faf.12673)。

<h3 id="krill-paper">相关应用：南极磷虾</h3>

Dong, S., Zhang, F. & Zhu, G. (2025). Length-dependent growth and mortality within each cohort cannot be ignored in stock assessment: a case study of Antarctic krill *Euphausia superba*. *Marine Ecology Progress Series*, **769**, 1–22. [doi:10.3354/meps14923](https://doi.org/10.3354/meps14923).

该研究将 ACL 与 ALSCL 应用于南极磷虾，考察同一队列内体长相关的生长和死亡率差异如何影响资源评估结果。

报告分析时，也应记录 ALSCL 版本、代码提交号、数据来源与拟合假设。
