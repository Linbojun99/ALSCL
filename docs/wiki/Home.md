# ALSCL 完整使用手册 · Complete user manual

<img src="https://github.com/Linbojun99/ALSCL/blob/main/ALSCLlogo.png?raw=true" alt="ALSCL logo" align="right" height="140" />

**调查体长数据的种群评估 · Survey catch-at-length stock assessment**

**ALSCL 2.0.0 · 简体中文与英文 · 15 章 / chapters**

本手册介绍模型原理、输入数据、模拟、拟合、诊断、绘图及回溯分析，提供双语 R 示例、Excel 填表示例、45 张 YTF 图和 7 张案例图。普通图使用 ggsci 配色；山脊图使用按年份排列的 viridis 渐变。

This manual covers model principles, data, simulation, fitting, diagnostics, plotting and retrospectives, with bilingual R examples, input worksheets, 45 YTF figures and seven case-study figures. Plots use ggsci colors; ridges use an ordered viridis gradient.

[开始使用 · Getting started](https://github.com/Linbojun99/ALSCL/wiki/01-Start) · [仓库主页 · Repository](https://github.com/Linbojun99/ALSCL/tree/main) · [单页书稿 · Single-page book](https://github.com/Linbojun99/ALSCL/blob/main/docs/BOOK.md) · [Excel 示例 · Workbook](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx)

## 目录 · Contents

| 章 / Chapter | 内容 / Contents |
|---|---|
| 01 | [开始使用 · Getting started](https://github.com/Linbojun99/ALSCL/wiki/01-Start) |
| 02 | [模型原理与时间尺度 · Model principles](https://github.com/Linbojun99/ALSCL/wiki/02-Theory) |
| 03 | [数据、Excel 与实测导入 · Data and observations](https://github.com/Linbojun99/ALSCL/wiki/03-Data-and-Excel) |
| 04 | [内置数据与首个拟合 · Built-in YTF](https://github.com/Linbojun99/ALSCL/wiki/04-Built-in-YTF) |
| 05 | [模拟与批量实验 · Simulation](https://github.com/Linbojun99/ALSCL/wiki/05-Simulation) |
| 06 | [初值、边界与拟合 · Fitting controls](https://github.com/Linbojun99/ALSCL/wiki/06-Fitting) |
| 07 | [诊断与结果解读 · Diagnostics](https://github.com/Linbojun99/ALSCL/wiki/07-Diagnostics) |
| 08 | [调查与种群图 · Survey and population plots](https://github.com/Linbojun99/ALSCL/wiki/08-Population-Plots) |
| 09 | [生长、死亡与残差图 · Growth, mortality and residuals](https://github.com/Linbojun99/ALSCL/wiki/09-Growth-and-Mortality) |
| 10 | [模型比较与真值 · Model comparisons](https://github.com/Linbojun99/ALSCL/wiki/10-Model-Comparisons) |
| 11 | [回溯分析 · Retrospectives](https://github.com/Linbojun99/ALSCL/wiki/11-Retrospectives) |
| 12 | [年度与季度案例 · Annual and quarterly cases](https://github.com/Linbojun99/ALSCL/wiki/12-Case-Studies) |
| 13 | [ggsci 配色与导出 · Colors and export](https://github.com/Linbojun99/ALSCL/wiki/13-Colors-and-Export) |
| 14 | [全部函数与参数 · Complete function reference](https://github.com/Linbojun99/ALSCL/wiki/14-Function-Reference) |
| 15 | [进阶、复现与故障处理 · Reproducibility](https://github.com/Linbojun99/ALSCL/wiki/15-Reproducibility) |

## 相关链接 · Links

**张帆教授 · Prof. Fan Zhang**

- [个人主页 · Faculty homepage](https://hyxy.shou.edu.cn/2021/0721/c18717a291861/page.htm)
- [GitHub · fzhang-shou](https://github.com/fzhang-shou)
- 邮箱 · Email: [f-zhang@shou.edu.cn](mailto:f-zhang@shou.edu.cn)

**董思宋 · Sisong Dong**

- [GitHub · dongworks97](https://github.com/dongworks97)
- [南极磷虾评估论文 · Antarctic krill assessment (Dong, Zhang & Zhu, 2025)](https://doi.org/10.3354/meps14923)

理论来源 / Reference: Zhang & Cadigan (2022), *Fish and Fisheries* 23,1121–1135. [Paper and Appendix S1](https://doi.org/10.1111/faf.12673).
