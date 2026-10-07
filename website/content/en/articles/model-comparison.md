# Model comparisons

First complete the [worked YTF example](first-model.html) to create `a`, `b`, `x` and `inputs`. Define `dat <- x$data.CatL` and load ggplot2 and patchwork. The [retrospective guide](retrospectives.html) creates `ra` and `rb`.

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## Population time series comparison

```r
# Draw this figure
plot_compare_ts(a, b, se = TRUE, ncol = 2)
```

[![Population time series comparison](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_ts.png?raw=true)


Compare common quantities over the same time window. Free y axes preclude amplitude comparisons across quantities.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_ts)



## Fishing mortality matrices

```r
# Draw this figure
plot_compare_F(a, b) + scale_x_continuous(breaks = c(2000, 2010, 2019))
```

[![Fishing mortality matrices](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F.png?raw=true)


ACL uses age and ALSCL length on the vertical axis; cells do not correspond one to one.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_F)



## Combined residual diagnostics

```r
# Draw this figure
plot_compare_residuals(a, b, dat)
```

[![Combined residual diagnostics](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_residuals.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_residuals.png?raw=true)


Includes histograms, QQ plots, annual residuals and length-bin boxplots. Annual error bars represent SD, not confidence intervals.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_residuals)



## Growth and reference-data interface

Reference data and layout helper:

```r
# Synthetic reference points, not paper observations
set.seed(42)
ref <- data.frame(Age = rep(1:15, each = 5))
ref$Length <- VB_func(60, .2, 1/60, ref$Age) + rnorm(nrow(ref), 0, 2)
growth_display <- function(p) {
  p[[1]] <- p[[1]] + labs(subtitle = paste(strsplit(
    p[[1]]$labels$subtitle, " | ", fixed=TRUE)[[1]], collapse="\n")) +
    theme(legend.text=element_text(size=7))
  p[[2]] <- p[[2]] + scale_x_continuous(breaks=c(1,5,10,15))
  p
}
```

```r
# Draw this figure
growth_display(plot_compare_growth(a, b, age_range = c(1, 15), ref_data = ref, ref_name = "Synthetic reference",
    nls_start = list(Linf = 60, k = 0.2, t0 = 1/60)))
```

[![Growth and reference-data interface](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_growth.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_growth.png?raw=true)


The reference curve uses synthetic Age and Length data, not empirical measurements from the paper.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_growth)



## Survey fits from both models

```r
# Draw this figure
plot_compare_CatL(a, b, dat, years = c(2000, 2005, 2010, 2015), ncol = 2)
```

[![Survey fits from both models](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_CatL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_CatL.png?raw=true)


Four selected years are shown. Omit years to show all periods. Points are observations and curves fitted medians.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_CatL)



## Fit metric comparison

```r
# Draw this figure
plot_compare_metrics(a, b, dat, ncol = 2) & scale_y_continuous(expand = expansion(mult = c(0,
    0.3)))
```

[![Fit metric comparison](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_metrics.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_metrics.png?raw=true)


The generating model is age based. This comparison cannot establish universal superiority; information criteria also require comparable data and likelihoods.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_metrics)



## Fixed survey catchability

```r
# Draw this figure
plot_compare_selectivity(a, b)
```

[![Fixed survey catchability](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_selectivity.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_selectivity.png?raw=true)


This function connects the input survey q; it does not estimate fishery selectivity. Identical inputs produce overlapping curves.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_selectivity)



## Summary F comparison apical

```r
# Draw this figure
plot_compare_annual_F(a, b, method = "apical") +
  labs(title = "Apical fishing mortality")
```

[![Summary F comparison apical](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_apical.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_apical.png?raw=true)


Apical is the maximum; mean is an unweighted group average. Despite its name, the function does not annualize quarterly rates.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_annual_F)



## Summary F comparison mean

```r
# Draw this figure
plot_compare_annual_F(a, b, method = "mean") +
  labs(title = "Mean fishing mortality")
```

[![Summary F comparison mean](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_mean.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/compare_F_mean.png?raw=true)


Apical is the maximum; mean is an unweighted group average. Despite its name, the function does not annualize quarterly rates.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_compare_annual_F)



## Simulation truth and conditional estimates

```r
# Align truth and fitted reports
years <- a$year
truth <- data.frame(Year=rep(years,3),
  Quantity=rep(c("B","SSB","Rec"),each=length(years)),
  Truth=c(x$truth$TB,x$truth$SSB,x$truth$Rec))
estimates <- rbind(
  data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(a$report$B,a$report$SSB,a$report$Rec),Model="ACL"),
  data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(b$report$B,b$report$SSB,b$report$Rec),Model="ALSCL"))
```

```r
# Draw this figure
ggplot(estimates, aes(Year, Value, color = Model)) + geom_line() + scale_color_manual(values = setNames(acl_theme("compare_colors"),
    c("ACL", "ALSCL"))) + geom_line(data = truth, aes(Year, Truth), inherit.aes = FALSE, linetype = 2) +
    facet_wrap(~Quantity, scales = "free_y", ncol = 1) + theme_bw() + labs(y = NULL)
```

[![Simulation truth and conditional estimates](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/truth.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/truth.png?raw=true)


Black dashed curves are generating truth. Biology is fixed; recruitment processes differ.
