# 开始使用

安装 ALSCL 2.0.0 后即可使用内置数据、拟合模型和绘图。以下代码在 R 中运行。


## 从 GitHub 安装

```r
# 安装 ALSCL
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
packageVersion("ALSCL")
data("YTF_example")
str(YTF_example[c("data.CatL", "data.wgt", "data.mat")])
```

首次拟合需要与 R 匹配的 C++ 工具链。安装程序自动处理依赖；读取 Excel 示例另需 `readxl`。包含 `source("scripts/...")` 的示例请从下载的仓库根目录执行。


[仓库主页](https://github.com/Linbojun99/ALSCL/blob/main/README.md) · [Excel 示例](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx)

## 准备编译工具

| 系统 | 要求 |
|---|---|
| Windows | 安装与 R 版本对应的 [Rtools](https://cran.r-project.org/bin/windows/Rtools/)。 |
| macOS | 安装 Xcode Command Line Tools，使用所安装 R 版本推荐的编译器设置。 |
| Linux | 安装 C++ 编译器、make 及 R 发行版所需的开发工具。 |

拟合函数会将包内 C++ 模板编译到可写的会话缓存，无需手动编译，也无需向已安装包目录写入文件。首次拟合可能比之后的拟合耗时更长。

## 包内示例与仓库脚本

加载 `library(ALSCL)` 后即可使用 `run_acl()`、`plot_SSB()` 等包函数。调用 `source("scripts/...")` 的示例还需要下载仓库，并将仓库根目录设为 R 工作目录，以正确读取脚本、Excel 模板与 CSV 文件。

## 开始第一个分析

继续阅读 [YTF 完整示例](first-model.html)，拟合两个模型并检查输出。分析自己的数据前，请先阅读 [数据准备](data-preparation.html)。函数参考页提供每个公开函数的完整参数和示例准备代码。
