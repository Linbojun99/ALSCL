# Growth, mortality and residuals

First complete the [worked YTF example](first-model.html) to create `a`, `b`, `x` and `inputs`. Define `dat <- x$data.CatL` and load ggplot2 and patchwork. The [retrospective guide](retrospectives.html) creates `ra` and `rb`.

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## Growth curve and interval

```r
# Draw this figure
plot_VB(b, age_range = c(1, 15), se = TRUE)
```

[![Growth curve and interval](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/VB.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/VB.png?raw=true)


Growth is fixed in this example, so its interval collapses. With estimated growth, intervals reflect parameter uncertainty, not individual length variation.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_VB)



## Age length conversion

```r
# Draw this figure
plot_pla(b) + scale_x_discrete(labels = as.character(1:15))
```

[![Age length conversion](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/pla.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/pla.png?raw=true)


Cells are length-bin probabilities conditional on age, not the growth transition matrix G. Age labels are class indices.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_pla)



## ACL fishing mortality year

```r
# Draw this figure
plot_fishing_mortality(a, type = "year",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL fishing mortality year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_year.png?raw=true)


F is an instantaneous rate per model step. Native F dimensions differ between the models.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ACL fishing mortality age

```r
# Draw this figure
plot_fishing_mortality(a, type = "age",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL fishing mortality age](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_age.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_age.png?raw=true)


F is an instantaneous rate per model step. Native F dimensions differ between the models.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ACL fishing mortality length

```r
# Draw this figure
plot_fishing_mortality(a, type = "length",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ACL fishing mortality length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ACL_length.png?raw=true)


In this version, ACL length is an alias of age; the output is still age-specific F.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL fishing mortality year

```r
# Draw this figure
plot_fishing_mortality(b, type = "year",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL fishing mortality year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_year.png?raw=true)


F is an instantaneous rate per model step. Native F dimensions differ between the models.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL fishing mortality age

```r
# Draw this figure
plot_fishing_mortality(b, type = "age",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL fishing mortality age](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_age.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_age.png?raw=true)


F is an instantaneous rate per model step. Native F dimensions differ between the models.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## ALSCL fishing mortality length

```r
# Draw this figure
plot_fishing_mortality(b, type = "length",
  se = TRUE, facet_ncol = 4, line_size = 0.6)
```

[![ALSCL fishing mortality length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/F_ALSCL_length.png?raw=true)


F is an instantaneous rate per model step. Native F dimensions differ between the models.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_fishing_mortality)



## Log residuals length

```r
# Draw this figure
plot_residuals(b, type = "length", facet_ncol = 4, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015))
```

[![Log residuals length](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_length.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_length.png?raw=true)


Observed minus fitted log index, without standardization. Smooth curves reveal structure, not significance.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_residuals)



## Log residuals year

```r
# Draw this figure
plot_residuals(b, type = "year", facet_ncol = 4, line_size = 0.6,
  x_breaks = c(10, 30, 50))
```

[![Log residuals year](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_year.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/residuals_year.png?raw=true)


Observed minus fitted log index, without standardization. Smooth curves reveal structure, not significance.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_residuals)



## Process deviations R TRUE

```r
# Draw this figure
plot_deviance(b, type = "R", log = TRUE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![Process deviations R TRUE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_TRUE.png?raw=true)


Log process deviations use zero as reference; these are not likelihood deviance statistics.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## Process deviations R FALSE

```r
# Draw this figure
plot_deviance(b, type = "R", log = FALSE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![Process deviations R FALSE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_R_FALSE.png?raw=true)


Exponentiated deviations use one as reference and asymmetric intervals; these are not original-scale residuals.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## Process deviations F TRUE

```r
# Draw this figure
plot_deviance(b, type = "F", log = TRUE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![Process deviations F TRUE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_TRUE.png?raw=true)


Log process deviations use zero as reference; these are not likelihood deviance statistics.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## Process deviations F FALSE

```r
# Draw this figure
plot_deviance(b, type = "F", log = FALSE, se = TRUE,
  facet_ncol = 3, point_size = 1.2, line_size = 0.5) +
  scale_x_continuous(breaks = c(2000, 2010, 2019), expand = expansion(mult = 0.04))
```

[![Process deviations F FALSE](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/deviation_F_FALSE.png?raw=true)


Exponentiated deviations use one as reference and asymmetric intervals; these are not original-scale residuals.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_deviance)



## Length composition ridges

```r
# Draw this figure
plot_ridges(b) # Default viridis gradient
```

[![Length composition ridges](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/ridges.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/ridges.png?raw=true)


Each period is normalized separately: observations left, fit right. This does not show total abundance trends. A live obj is required.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_ridges)
