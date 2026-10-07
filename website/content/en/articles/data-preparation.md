# Data and observations

Download the [Excel example workbook](https://github.com/Linbojun99/ALSCL/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx). `Instructions` contains bilingual directions, and the other three sheets contain the exact YTF example. Save a copy before replacing values.

|Location|Entry rule|
|---|---|
| Sheet names | Exactly `data.CatL`, `data.wgt`, `data.mat` |
| A1 | `LengthBin`; no title rows above the header |
| A2:A… | One bin per row; example midpoints are 6, 8, …, 50 |
| B1, C1, … | Increasing numeric times: 2000, 2001 for annual data; 2000, 2000.25 for quarters |
|Numeric cells|Value for that bin and period|
|Alignment|Exact same bins, periods and order|

### Survey number index: `data.CatL`

![Survey data worksheet](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_CatL.png?raw=true)


Enter counts or a consistently standardized number index. For CPUE, maintain comparable effort, units and standardization. Decimals are valid. Use blank/`NA` for unknown observations; do not use zero, dashes, `<1` or text labels as missing-value substitutes.

### Mean individual weight: `data.wgt`

![Mean-weight worksheet](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_wgt.png?raw=true)


Enter mean weight per fish, not the total weight in the bin. Use one unit and document it separately. Repeating a defensible length-weight relationship across periods is possible but must be reported. Missing weights are not accepted.

### Maturity proportion: `data.mat`

![Maturity worksheet](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/excel_mat.png?raw=true)


Enter proportions in [0,1], e.g. `0.6` for 60%. Both 0 and 1 are valid maturity values. Do not insert units, totals, merged cells or prose inside the three data grids.

Scientific notation avoids displaying tiny positive values as zero; it does not round the stored values.


The example uses midpoints 6:50 by 2 and 22 explicit internal boundaries 7:49. Supply the actual boundaries for unequal bins. Interval labels must be contiguous and non-overlapping; do not silently close gaps in legacy labels. Retain an entirely missing period as an NA survey column while supplying weight and maturity, so the temporal grid remains regular.

## Importing your observations


These are guide helpers, sourced from the repository, not package exports. The Excel reader preserves original time headers and converts numeric cells explicitly.

```r
library(ALSCL)
data("YTF_example")
install.packages("readxl") # Install once
source("scripts/real_data_workflow.R", encoding = "UTF-8")
# Start with the supplied workbook
my_inputs <- read_survey_excel("docs/data/ALSCL_Data_Entry_Example.xlsx")
validate_survey_tables(my_inputs, growth_step = 1, zero_action = "error")
# Three named CSV files also work
csv_inputs <- read_survey_csv("docs/data")
# Use the workbook's matching biology
config <- c(YTF_example$fit_args, YTF_example$fit_config$alscl,
            list(zero_action = "error", train_times = 2, silent = TRUE))
# This refits and requires a compiler
my_fit <- fit_survey_tables(my_inputs, config, model_type = "alscl")
```


For your own workbook, replace the path and supply species-specific biology, boundaries, starts and fixed parameters. Equal table dimensions do not justify reusing YTF biology. Growth, weight, maturity and survey catchability need independent support.


The guide validator rejects survey zeros by default. Package fits retain the historical default `zero_action="missing"`, excluding zeros from the log likelihood. Excluding genuine zero catches changes inference. This package does not implement a zero-inflated observation model; arbitrary small constants are not a principled substitute.

### Legacy data inspection

```r
data("example_data")
names(example_data)
head(example_data$data.CatL)
# Inspect bins and zeros before fitting
example_data$data.CatL[[1]]
sum(as.matrix(example_data$data.CatL[-1]) == 0, na.rm = TRUE)
```


The legacy tables are preserved unchanged from ALSCL. Their old help calls them anonymized survey data, but species, units, original source and zero meanings are unverified. Their gapped interval labels are rejected by the current bin parser. Recover metadata and actual boundaries before assessment; these are not verified empirical data from the cited paper.
