# Getting started

Install ALSCL 2.0.0 to use the bundled data, fit models and create plots. Run the following in R.

## Install from GitHub

```r
# Install ALSCL
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
packageVersion("ALSCL")
data("YTF_example")
str(YTF_example[c("data.CatL", "data.wgt", "data.mat")])
```


The first fit needs an R-compatible C++ toolchain. Dependencies are installed automatically; Excel examples additionally use readxl. Run examples containing `source("scripts/...")` from the downloaded repository root.

[Repository](https://github.com/Linbojun99/ALSCL/blob/main/README.md) · [Workbook](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx)

## Prepare the compiler

| Platform | Requirement |
|---|---|
| Windows | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version. |
| macOS | Install the Xcode Command Line Tools; use the R distribution's recommended compiler configuration. |
| Linux | Install a C++ compiler, make and the development tools required by your R distribution. |

The fitting functions compile the bundled C++ templates into a writable session cache. You do not need to compile the models manually or write into the installed package directory. The first fit may take longer than later fits.

## Package examples and repository scripts

Functions such as `run_acl()` and `plot_SSB()` are available after `library(ALSCL)`. Examples that call `source("scripts/...")` additionally require a downloaded copy of the repository. Open that repository as your R working directory so the paths to scripts, input worksheets and example CSV files resolve correctly.

## Your first analysis

Continue to the [worked YTF example](first-model.html) to fit both models and inspect their outputs. For your own data, follow [data preparation](data-preparation.html) before fitting. The function reference includes a complete argument list and example setup for each public function.
